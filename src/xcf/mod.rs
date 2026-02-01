/*
This module provides the VCF/BCF-specific I/O functions.
*/
use crate::io::{consolidate_list, get_gt_index, gt2gcount, GenoData};
use anyhow::Result;
use counter::Counter;
use itertools::MultiUnzip;
use rust_htslib::bcf::{
    record::{Genotype, GenotypeAllele},
    IndexedReader, Read, Reader,
};
use std::{fmt::Display, path::Path};

// Data structure
// Define multiple readers for the indexed and unindexed XCF file
pub enum XcfReader {
    Indexed(IndexedReader),
    Readthrough(Reader),
}

/// Read an XCF (VCF/BCF) file, either indexed or unindexed.
///
/// # Examples
///
/// ```ignore
/// let reader = xpclrs::xcf::read_xcf("in.bcf", true).unwrap();
/// ```
pub fn read_xcf<P: AsRef<Path> + Display>(path: P, has_index: bool) -> Result<XcfReader> {
    let xcf_reader: XcfReader = if has_index {
        XcfReader::Indexed(
            IndexedReader::from_path(path).expect("Cannot load indexed BCF/VCF file"),
        )
    } else {
        XcfReader::Readthrough(Reader::from_path(path).expect("Cannot load BCF/VCF file"))
    };
    Ok(xcf_reader)
}

/// Process an indexed XCF file into `GenoData`.
///
/// # Examples
///
/// ```ignore
/// let data = xpclrs::xcf::indexed_xcf("in.bcf".to_string(), &s1, &s2, "1", 0, None, (None, None, None, 1)).unwrap();
/// ```
pub fn indexed_xcf(
    xcf_fn: String,
    s1: &[String],
    s2: &[String],
    chrom: &str,
    start: u64,
    end: Option<u64>,
    (phased, rrate, gdistkey, n_threads): (Option<bool>, Option<f64>, Option<String>, usize),
) -> Result<GenoData> {
    log::info!("Indexed reader.");
    // Prepare the indexed reader
    let mut reader = IndexedReader::from_path(xcf_fn).expect("Cannot load indexed BCF/VCF file");
    reader
        .set_threads(n_threads)
        .expect("Failed to set threads");
    let rrate = rrate.unwrap_or(1e-8);
    // Resolve options once to avoid per-record branching.
    let phased = phased.unwrap_or(false);

    // Load the XCF file
    let xcf_header = reader.header().clone();
    log::info!("Samples in VCF: {}", xcf_header.sample_count());

    // Load sample lists as an u8 array
    let s1 = consolidate_list(&xcf_header.samples(), s1).expect("Failed to subset sampleA");
    let s2 = consolidate_list(&xcf_header.samples(), s2).expect("Failed to subset sampleB");

    // Fetch the indices of each sample in each list
    let i1 = get_gt_index(&xcf_header.samples(), &s1).expect("Failed to get indeces of sampleA");
    let i2 = get_gt_index(&xcf_header.samples(), &s2).expect("Failed to get indeces of sampleB");

    // Print number of samples
    log::info!("Samples A: {}", i1.len());
    log::info!("Samples B: {}", i2.len());

    // Dies if no samples are retained
    if s1.is_empty() || s2.is_empty() {
        eprintln!("No samples found in the lists.");
        std::process::exit(1);
    }

    // Start loading the genotypes here
    // First, find the sequence index
    let rid = reader
        .header()
        .name2rid(chrom.as_bytes())
        .expect("RID not found");
    log::info!("Chromosome {chrom} (ID: {rid})");

    // Jump to target position in place
    let _ = reader.fetch(rid, start, end);
    // Load the records, defining the counters of how many sites we skip
    let mut multiallelic = 0;
    let mut monom_gt2 = 0;
    let mut miss_gt1 = 0;
    let mut miss_gt2 = 0;
    let mut pass = 0;
    let mut skipped = 0;
    let mut tot = 0;

    let (positions, gt1_data, gt2_data, gd_data): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) = if phased {
        reader
            .records()
            .filter_map(|r| {
                let record = r.ok()?;
                tot += 1;
                let genotypes = record.genotypes().expect("Cannot fetch the genotypes");
                let gt1_g = i1
                    .iter()
                    .map(|i| genotypes.get(*i))
                    .collect::<Vec<Genotype>>();
                let gt2_g = i2
                    .iter()
                    .map(|i| genotypes.get(*i))
                    .collect::<Vec<Genotype>>();

                // Count alleles from both populations in a single pass
                let mut alleles1: Counter<u32> = Counter::new();
                let mut alleles2: Counter<u32> = Counter::new();

                for g in &gt1_g {
                    for a in g.iter().filter_map(|a| a.index()) {
                        alleles1[&a] += 1;
                    }
                }
                for g in &gt2_g {
                    for a in g.iter().filter_map(|a| a.index()) {
                        alleles2[&a] += 1;
                    }
                }

                // Union both sets
                let all_alleles: Counter<u32> = alleles1
                    .iter()
                    .chain(alleles2.iter())
                    .map(|(&k, &v)| (k, v))
                    .collect();

                // Define genetic position
                let gd = match &gdistkey {
                    Some(key) => record
                        .info(key.as_bytes())
                        .float()
                        .ok()
                        .flatten()
                        .expect("Missing info field for genetic position")[0]
                        as f64,
                    None => record.pos() as f64 * rrate,
                };

                // Perform filtering and counting
                if all_alleles.len() > 2 {
                    skipped += 1;
                    multiallelic += 1;
                    None
                } else if alleles1.is_empty() || alleles2.is_empty() {
                    skipped += 1;
                    if alleles1.is_empty() {
                        miss_gt1 += 1;
                    };
                    if alleles2.is_empty() {
                        miss_gt2 += 1;
                    };
                    None
                } else if alleles2.len() == 1 || alleles2.values().min().copied()? == 1 {
                    skipped += 1;
                    monom_gt2 += 1;
                    None
                } else {
                    pass += 1;
                    // Define reference allele as the minimum allele index (consistent with methods)
                    let ref_ix = *all_alleles
                        .keys()
                        .min()
                        .expect("Can't compute reference allele index");
                    // Encode genotypes to compact i8 counts relative to ref allele
                    let mut gt1 = Vec::with_capacity(gt1_g.len());
                    for gt in gt1_g {
                        gt1.push(gt2gcount(gt, ref_ix));
                    }
                    // If phased, store haplotypes in gt2; else store dosages
                    let mut gt2 = Vec::with_capacity(gt2_g.len() * 2);
                    for gt in &gt2_g {
                        for a in gt.iter() {
                            gt2.push(match a {
                                GenotypeAllele::PhasedMissing => -9_i8,
                                GenotypeAllele::UnphasedMissing => -9_i8,
                                _ => a.index().unwrap() as i8,
                            });
                        }
                    }
                    Some((record.pos() as usize, gt1, gt2, gd))
                }
            })
            .multiunzip()
    } else {
        reader
            .records()
            .filter_map(|r| {
                let record = r.ok()?;
                tot += 1;
                let genotypes = record.genotypes().expect("Cannot fetch the genotypes");
                let gt1_g = i1
                    .iter()
                    .map(|i| genotypes.get(*i))
                    .collect::<Vec<Genotype>>();
                let gt2_g = i2
                    .iter()
                    .map(|i| genotypes.get(*i))
                    .collect::<Vec<Genotype>>();

                // Count alleles from both populations in a single pass
                let mut alleles1: Counter<u32> = Counter::new();
                let mut alleles2: Counter<u32> = Counter::new();

                for g in &gt1_g {
                    for a in g.iter().filter_map(|a| a.index()) {
                        alleles1[&a] += 1;
                    }
                }
                for g in &gt2_g {
                    for a in g.iter().filter_map(|a| a.index()) {
                        alleles2[&a] += 1;
                    }
                }

                // Union both sets
                let all_alleles: Counter<u32> = alleles1
                    .iter()
                    .chain(alleles2.iter())
                    .map(|(&k, &v)| (k, v))
                    .collect();

                // Define genetic position
                let gd = match &gdistkey {
                    Some(key) => record
                        .info(key.as_bytes())
                        .float()
                        .ok()
                        .flatten()
                        .expect("Missing info field for genetic position")[0]
                        as f64,
                    None => record.pos() as f64 * rrate,
                };

                // Perform filtering and counting
                if all_alleles.len() > 2 {
                    skipped += 1;
                    multiallelic += 1;
                    None
                } else if alleles1.is_empty() || alleles2.is_empty() {
                    skipped += 1;
                    if alleles1.is_empty() {
                        miss_gt1 += 1;
                    };
                    if alleles2.is_empty() {
                        miss_gt2 += 1;
                    };
                    None
                } else if alleles2.len() == 1 || alleles2.values().min().copied()? == 1 {
                    skipped += 1;
                    monom_gt2 += 1;
                    None
                } else {
                    pass += 1;
                    // Define reference allele as the minimum allele index (consistent with methods)
                    let ref_ix = *all_alleles
                        .keys()
                        .min()
                        .expect("Can't compute reference allele index");
                    // Encode genotypes to compact i8 counts relative to ref allele
                    let mut gt1 = Vec::with_capacity(gt1_g.len());
                    for gt in gt1_g {
                        gt1.push(gt2gcount(gt, ref_ix));
                    }
                    // If phased, store haplotypes in gt2; else store dosages
                    let mut gt2 = Vec::with_capacity(gt2_g.len());
                    for gt in gt2_g {
                        gt2.push(gt2gcount(gt, ref_ix));
                    }
                    Some((record.pos() as usize, gt1, gt2, gd))
                }
            })
            .multiunzip()
    };

    // Assess everything looks good
    if gt1_data.len() != gt2_data.len() {
        panic!("Inconsistent data")
    };
    if positions.len() != gt2_data.len() {
        panic!("Inconsistent data")
    };
    if positions.len() as i32 != pass {
        panic!("Inconsistent data")
    };
    if tot != (pass + skipped) {
        panic!("Inconsistent counts")
    }
    // Print some info
    log::info!("Processed {tot} variants");
    log::info!("Loaded {pass} variants");
    log::info!("Skipped {skipped} variants because:");
    log::info!(" - {multiallelic} multiallelic");
    log::info!(" - {monom_gt2} monomorphic/singleton in pop B");
    log::info!(" - {miss_gt1} all-missing in pop A");
    log::info!(" - {miss_gt2} all-missing in pop B");
    Ok(GenoData {
        positions,
        gt1: gt1_data,
        gt2: gt2_data,
        gdistances: gd_data,
    })
}

/// Process an unindexed XCF file into `GenoData`.
///
/// # Examples
///
/// ```ignore
/// let data = xpclrs::xcf::readthrough_xcf("in.bcf".to_string(), &s1, &s2, "1", 0, None, (None, None, None, 1)).unwrap();
/// ```
pub fn readthrough_xcf(
    xcf_fn: String,
    s1: &[String],
    s2: &[String],
    chrom: &str,
    start: u64,
    end: Option<u64>,
    (phased, rrate, gdistkey, n_threads): (Option<bool>, Option<f64>, Option<String>, usize),
) -> Result<GenoData> {
    log::info!("Streamed reader.");
    log::info!("This is substantially slower than the indexed one.");
    log::info!("Consider generating an index for your BCF/VCF file.");
    // Prepare the indexed reader
    let mut reader = Reader::from_path(xcf_fn).expect("Cannot load unindexed/streamed BCF/VCF file");
    reader
        .set_threads(n_threads)
        .expect("Failed to set threads");
    let end = end.unwrap_or(999999999);
    let rrate = rrate.unwrap_or(1e-8);
    // Resolve options once to avoid per-record branching.
    let phased = phased.unwrap_or(false);

    // Load the XCF file
    let xcf_header = reader.header().clone();
    log::info!("Samples in VCF: {}", xcf_header.sample_count());

    // Load sample lists as an u8 array
    let s1 = consolidate_list(&xcf_header.samples(), s1).expect("Failed to subset sampleA");
    let s2 = consolidate_list(&xcf_header.samples(), s2).expect("Failed to subset sampleB");

    // Fetch the indices of each sample in each list
    let i1 = get_gt_index(&xcf_header.samples(), &s1).expect("Failed to get indeces of sampleA");
    let i2 = get_gt_index(&xcf_header.samples(), &s2).expect("Failed to get indeces of sampleB");

    // Print number of samples
    log::info!("Samples A: {}", s1.len());
    log::info!("Samples B: {}", s2.len());

    // Dies if no samples are retained
    if s1.is_empty() || s2.is_empty() {
        eprintln!("No samples found in the lists.");
        std::process::exit(1);
    }

    // Start loading the genotypes here
    // First, find the sequence index
    let rid = reader
        .header()
        .name2rid(chrom.as_bytes())
        .unwrap_or_else(|_| panic!("Chromosome ID not found {chrom}"));
    log::info!("Chromosome {chrom} (ID: {rid})");

    // Load the records, defining the counters of how many sites we skip
    let mut multiallelic = 0;
    let mut monom_gt2 = 0;
    let mut miss_gt1 = 0;
    let mut miss_gt2 = 0;
    let mut pass = 0;
    let mut skipped = 0;
    let mut tot = 0;
    let (positions, gt1_data, gt2_data, gd_data): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) = if phased {
        reader
            .records()
            .filter(|r| {
                let record = r.as_ref().unwrap();
                let pos = record.pos() as u64;
                record.rid().unwrap() == rid && (pos >= start && pos < end)
            })
            .filter_map(|r| {
                let record = r.ok()?;
                tot += 1;
                let genotypes = record.genotypes().expect("Cannot fetch the genotypes");
                let gt1_g = i1
                    .iter()
                    .map(|i| genotypes.get(*i))
                    .collect::<Vec<Genotype>>();
                let gt2_g = i2
                    .iter()
                    .map(|i| genotypes.get(*i))
                    .collect::<Vec<Genotype>>();

                // Count alleles from both populations in a single pass
                let mut alleles1: Counter<u32> = Counter::new();
                let mut alleles2: Counter<u32> = Counter::new();

                for g in &gt1_g {
                    for a in g.iter().filter_map(|a| a.index()) {
                        alleles1[&a] += 1;
                    }
                }
                for g in &gt2_g {
                    for a in g.iter().filter_map(|a| a.index()) {
                        alleles2[&a] += 1;
                    }
                }

                // Union both sets
                let all_alleles: Counter<u32> = alleles1
                    .iter()
                    .chain(alleles2.iter())
                    .map(|(&k, &v)| (k, v))
                    .collect();

                // Define genetic position
                let gd = match &gdistkey {
                    Some(key) => record
                        .info(key.as_bytes())
                        .float()
                        .ok()
                        .flatten()
                        .expect("Missing info field for genetic position")[0]
                        as f64,
                    None => record.pos() as f64 * rrate,
                };

                // Perform filtering and counting
                if all_alleles.len() > 2 {
                    skipped += 1;
                    multiallelic += 1;
                    None
                } else if alleles1.is_empty() || alleles2.is_empty() {
                    skipped += 1;
                    if alleles1.is_empty() {
                        miss_gt1 += 1;
                    };
                    if alleles2.is_empty() {
                        miss_gt2 += 1;
                    };
                    None
                } else if alleles2.len() == 1 || alleles2.values().min().copied()? == 1 {
                    skipped += 1;
                    monom_gt2 += 1;
                    None
                } else {
                    pass += 1;
                    // Define reference allele as the minimum allele index
                    let ref_ix = *all_alleles
                        .keys()
                        .min()
                        .expect("Can't compute reference allele index");
                    // Encode genotypes to compact i8 counts relative to ref allele
                    let mut gt1 = Vec::with_capacity(gt1_g.len());
                    for gt in gt1_g {
                        gt1.push(gt2gcount(gt, ref_ix));
                    }
                    // If phased, store haplotypes in gt2; else store dosages
                    let mut gt2 = Vec::with_capacity(gt2_g.len() * 2);
                    for gt in &gt2_g {
                        for a in gt.iter() {
                            gt2.push(match a {
                                GenotypeAllele::PhasedMissing => -9_i8,
                                GenotypeAllele::UnphasedMissing => -9_i8,
                                _ => a.index().unwrap() as i8,
                            });
                        }
                    }
                    Some((record.pos() as usize, gt1, gt2, gd))
                }
            })
            .multiunzip()
    } else {
        reader
            .records()
            .filter(|r| {
                let record = r.as_ref().unwrap();
                let pos = record.pos() as u64;
                record.rid().unwrap() == rid && (pos >= start && pos < end)
            })
            .filter_map(|r| {
                let record = r.ok()?;
                tot += 1;
                let genotypes = record.genotypes().expect("Cannot fetch the genotypes");
                let gt1_g = i1
                    .iter()
                    .map(|i| genotypes.get(*i))
                    .collect::<Vec<Genotype>>();
                let gt2_g = i2
                    .iter()
                    .map(|i| genotypes.get(*i))
                    .collect::<Vec<Genotype>>();

                // Count alleles from both populations in a single pass
                let mut alleles1: Counter<u32> = Counter::new();
                let mut alleles2: Counter<u32> = Counter::new();

                for g in &gt1_g {
                    for a in g.iter().filter_map(|a| a.index()) {
                        alleles1[&a] += 1;
                    }
                }
                for g in &gt2_g {
                    for a in g.iter().filter_map(|a| a.index()) {
                        alleles2[&a] += 1;
                    }
                }

                // Union both sets
                let all_alleles: Counter<u32> = alleles1
                    .iter()
                    .chain(alleles2.iter())
                    .map(|(&k, &v)| (k, v))
                    .collect();

                // Define genetic position
                let gd = match &gdistkey {
                    Some(key) => record
                        .info(key.as_bytes())
                        .float()
                        .ok()
                        .flatten()
                        .expect("Missing info field for genetic position")[0]
                        as f64,
                    None => record.pos() as f64 * rrate,
                };

                // Perform filtering and counting
                if all_alleles.len() > 2 {
                    skipped += 1;
                    multiallelic += 1;
                    None
                } else if alleles1.is_empty() || alleles2.is_empty() {
                    skipped += 1;
                    if alleles1.is_empty() {
                        miss_gt1 += 1;
                    };
                    if alleles2.is_empty() {
                        miss_gt2 += 1;
                    };
                    None
                } else if alleles2.len() == 1 || alleles2.values().min().copied()? == 1 {
                    skipped += 1;
                    monom_gt2 += 1;
                    None
                } else {
                    pass += 1;
                    // Define reference allele as the minimum allele index
                    let ref_ix = *all_alleles
                        .keys()
                        .min()
                        .expect("Can't compute reference allele index");
                    // Encode genotypes to compact i8 counts relative to ref allele
                    let mut gt1 = Vec::with_capacity(gt1_g.len());
                    for gt in gt1_g {
                        gt1.push(gt2gcount(gt, ref_ix));
                    }
                    // If phased, store haplotypes in gt2; else store dosages
                    let mut gt2 = Vec::with_capacity(gt2_g.len());
                    for gt in gt2_g {
                        gt2.push(gt2gcount(gt, ref_ix));
                    }
                    Some((record.pos() as usize, gt1, gt2, gd))
                }
            })
            .multiunzip()
    };

    // Assess everything looks good
    if gt1_data.len() != gt2_data.len() {
        panic!("Inconsistent data")
    };
    if positions.len() != gt2_data.len() {
        panic!("Inconsistent data")
    };
    if positions.len() as i32 != pass {
        panic!("Inconsistent data")
    };
    if tot != (pass + skipped) {
        panic!("Inconsistent counts")
    }
    // Print some info

    log::info!("Processed {tot} variants");
    log::info!("Loaded {pass} variants");
    log::info!("Skipped {skipped} variants because:");
    log::info!(" - {multiallelic} multiallelic");
    log::info!(" - {monom_gt2} monomorphic/singleton in pop B");
    log::info!(" - {miss_gt1} all-missing in pop A");
    log::info!(" - {miss_gt2} all-missing in pop B");
    Ok(GenoData {
        positions,
        gt1: gt1_data,
        gt2: gt2_data,
        gdistances: gd_data,
    })
}
