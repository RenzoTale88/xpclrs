/*
This module provides the I/O functions (e.g. VCF/BCF readers, text readers etc).
*/
use anyhow::Result;
use counter::Counter;
use itertools::MultiUnzip;
use rayon::prelude::*;
use rust_htslib::bcf::{self, record::{Genotype, GenotypeAllele}, IndexedReader, Read, Reader};
use std::{
    collections::HashSet,
    fmt::Display,
    io::{BufRead, BufReader},
    path::Path,
};

// Define multople readers for the indexed and unindexed XCF file
pub enum XcfReader {
    Indexed(bcf::IndexedReader),
    Readthrough(bcf::Reader),
}

/// Same as smakcr, but single threaded for now
pub fn read_xcf<P: AsRef<Path> + Display>(path: P, has_index: bool) -> Result<XcfReader> {
    let xcf_reader: XcfReader = if has_index {
        XcfReader::Indexed(
            bcf::IndexedReader::from_path(path).expect("Cannot load indexed BCF/VCF file"),
        )
    } else {
        XcfReader::Readthrough(bcf::Reader::from_path(path).expect("Cannot load BCF/VCF file"))
    };
    Ok(xcf_reader)
}

/// Same as smakcr, but single threaded for now
pub fn read_file<P: AsRef<Path> + Display>(
    path: P,
) -> Result<impl Iterator<Item = Result<String>>> {
    let sniffed_reader: std::result::Result<
        (Box<dyn std::io::Read>, niffler::Format),
        niffler::Error,
    > = niffler::from_path(&path);
    // if the file has fewer than 5 bytes, `niffler` can't sniff the compression format and will
    // return a `FileTooShort` error; this could be due to
    // * an empty file
    // * a file containing only a single FASTA record with the ID consisting only of a single
    //   character and the sequence being empty
    // * a file containing a single sequence with a one-character ID and one-character sequence and
    //   missing newline character at the end
    // we don't want to fail at this stage in these cases and thus handle the `FileTooShort` error
    // separately
    let reader = match sniffed_reader {
        Ok(rdr) => Ok(rdr.0),
        Err(e) => Err(e),
    }?;
    let buffer = BufReader::new(reader);
    let records = buffer
        .lines()
        .map(|item: std::result::Result<String, std::io::Error>| {
            item.map_err(|e| anyhow::anyhow!("Error reading line: {}", e))
        });
    Ok(records)
}

// Consolidate the list
fn consolidate_list(full_list: &Vec<&[u8]>, subset: &[String]) -> Result<Vec<String>> {
    let filtered = subset
        .iter()
        .filter(|&s| full_list.contains(&s.as_bytes()))
        .cloned()
        .collect::<Vec<String>>();
    Ok(filtered)
}

// Function to get the index of each sample in the lists
fn get_gt_index(full_list: &Vec<&[u8]>, subset: &[String]) -> Result<Vec<usize>> {
    let subset_set: HashSet<_> = subset.iter().collect();
    let indices: Vec<usize> = full_list.iter()
        .enumerate()
        .filter_map(|(idx, &val)| {
            if subset_set.contains(&String::from_utf8(val.to_owned()).unwrap()) { Some(idx) } else { None }
        })
        .collect();
    Ok(indices)
}


// Process indexed XCF file
fn indexed_xcf(
    xcf_fn: String,
    s1: &[String],
    s2: &[String],
    chrom: String,
    start: u64,
    end: Option<u64>,
    _gdistkey: Option<String>,
) -> Result<(Vec<i64>, Vec<Vec<Genotype>>, Vec<Vec<Genotype>>)> {
    log::info!("Indexed reader.");
    // Prepare the indexed reader
    let mut reader = IndexedReader::from_path(xcf_fn).expect("Cannot load indexed BCF/VCF file");

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
    let mut monom_gt1 = 0;
    let mut monom_gt2 = 0;
    let mut miss_gt1 = 0;
    let mut miss_gt2 = 0;
    let mut pass = 0;
    let mut skipped = 0;
    let mut tot = 0;
    let (positions, gt1_data, gt2_data): (Vec<_>, Vec<_>, Vec<_>) = reader.records().filter_map(|r| {
        let record = r.ok()?;
        tot += 1;
        let genotypes = record.genotypes().expect("Cannot fetch the genotypes");
        let gt1 = i1.iter().map(|i| genotypes.get(*i)).collect::<Vec<Genotype>>();
        let gt2 = i2.iter().map(|i| genotypes.get(*i)).collect::<Vec<Genotype>>();
        // Check that the site is biallelic
        let alleles1: Counter<u32> = gt1.iter()
            .flat_map(|g| g.iter().filter_map(|&a: &GenotypeAllele| a.index()))
            .collect::<Counter<u32>>();

        let alleles2: Counter<u32> = gt2.iter()
            .flat_map(|g| g.iter().filter_map(|a| a.index()))
            .collect::<Counter<u32>>();

        // Union both sets
        let mut all_alleles: Counter<_> = alleles1.clone();
        all_alleles.extend(&alleles2);

        // Perform filtering and counting
        if all_alleles.len() > 2 {
            skipped += 1;
            multiallelic += 1;
            None
        } else if alleles1.len() == 0 || alleles2.len() == 0 {
            skipped += 1;
            if alleles1.len() == 0 {
                miss_gt1 += 1;
            };
            if alleles2.len() == 0 {
                miss_gt2 += 1;
            };
            None
        } else {
            if alleles1.len() == 1 || alleles1.values().min().copied()? == 1 || alleles2.len() == 1 || alleles2.values().min().copied()? == 1 {
                skipped += 1;
                if alleles1.len() == 1 || alleles1.values().min().copied()? == 1 {
                    monom_gt1 += 1;
                };
                if alleles2.len() == 1 || alleles2.values().min().copied()? == 1 {
                    monom_gt2 += 1;
                }
                None
            } else {
                pass += 1;
                Some((record.pos(), gt1, gt2))
            }
        }
    }).multiunzip();

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
    println!("{tot} = {pass} {skipped}");
    if tot != (pass + skipped){
        panic!("Inconsistent counts")
    }
    // Print some info
    log::info!("Processed {} variants", tot);
    log::info!("Loaded {} variants", pass);
    log::info!("Skipped {} variants because:", skipped);
    log::info!(" - {} multiallelic", multiallelic);
    log::info!(" - {} monomorphic in pop A", monom_gt1);
    log::info!(" - {} monomorphic in pop B", monom_gt2);
    log::info!(" - {} all-missing in pop A", miss_gt1);
    log::info!(" - {} all-missing in pop B", miss_gt2);
    Ok((positions, gt1_data, gt2_data))
}

// Process unindexed XCF
fn readthrough_xcf(
    xcf_fn: String,
    s1: &[String],
    s2: &[String],
    chrom: String,
    start: u64,
    end: Option<u64>,
    _gdistkey: Option<String>,
) -> Result<
        (Vec<i64>, Vec<Vec<Genotype>>, Vec<Vec<Genotype>>)
> {
    log::info!("Streamed reader.");
    log::info!("This is substantially slower than the indexed one.");
    log::info!("Consider generating an index for your BCF/VCF file.");
    // Prepare the indexed reader
    let mut reader = Reader::from_path(xcf_fn).expect("Cannot load indexed BCF/VCF file");
    let end = end.unwrap_or(999999999);

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
    let mut monom_gt1 = 0;
    let mut monom_gt2 = 0;
    let mut miss_gt1 = 0;
    let mut miss_gt2 = 0;
    let mut pass = 0;
    let mut skipped = 0;
    let mut tot = 0;
    let (positions, gt1_data, gt2_data): (Vec<_>, Vec<_>, Vec<_>) = reader
        .records()
        .filter(|r| {
                let record = r.as_ref().unwrap();
                let pos = record.pos() as u64;
                record.rid().unwrap() == rid && (pos >= start && pos < end )
            }
        )
        .filter_map(|r| {
            let record = r.ok()?;
            tot += 1;
            let genotypes = record.genotypes().expect("Cannot fetch the genotypes");
            let gt1 = i1.iter().map(|i| genotypes.get(*i)).collect::<Vec<Genotype>>();
            let gt2 = i2.iter().map(|i| genotypes.get(*i)).collect::<Vec<Genotype>>();

            // Check that the site is biallelic
            let alleles1: Counter<u32> = gt1.iter()
                .flat_map(|g| g.iter().filter_map(|&a: &GenotypeAllele| a.index()))
                .collect::<Counter<u32>>();

            let alleles2: Counter<u32> = gt2.iter()
                .flat_map(|g| g.iter().filter_map(|a| a.index()))
                .collect::<Counter<u32>>();

            // Union both sets
            let mut all_alleles: Counter<_> = alleles1.clone();
            all_alleles.extend(&alleles2);

            // Perform filtering and counting
            if all_alleles.len() > 2 {
                skipped += 1;
                multiallelic += 1;
                None
            } else if alleles1.len() == 0 || alleles2.len() == 0 {
                skipped += 1;
                if alleles1.len() == 0 {
                    miss_gt1 += 1;
                };
                if alleles2.len() == 0 {
                    miss_gt2 += 1;
                };
                None
            } else {
                if alleles1.len() == 1 || alleles1.values().min().copied()? == 1 || alleles2.len() == 1 || alleles2.values().min().copied()? == 1 {
                    skipped += 1;
                    if alleles1.len() == 1 || alleles1.values().min().copied()? == 1 {
                        monom_gt1 += 1;
                    };
                    if alleles2.len() == 1 || alleles2.values().min().copied()? == 1 {
                        monom_gt2 += 1;
                    }
                    None
                } else {
                    pass += 1;
                    Some((record.pos(), gt1, gt2))
                }
            }
        }).multiunzip();

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
    println!("{tot} = {pass} {skipped}");
    if tot != (pass + skipped){
        panic!("Inconsistent counts")
    }
    // Print some info
    log::info!("Processed {} variants", tot);
    log::info!("Loaded {} variants", pass);
    log::info!("Skipped {} variants because:", skipped);
    log::info!(" - {} multiallelic", multiallelic);
    log::info!(" - {} monomorphic in pop A", monom_gt1);
    log::info!(" - {} monomorphic in pop B", monom_gt2);
    log::info!(" - {} all-missing in pop A", miss_gt1);
    log::info!(" - {} all-missing in pop B", miss_gt2);
    Ok((positions, gt1_data, gt2_data))
}

// Load the genotypes for the given samples
pub fn process_xcf(
    xcf_fn: String,
    s1: &[String],
    s2: &[String],
    chrom: String,
    start: Option<u64>,
    end: Option<u64>,
    _gdistkey: Option<String>,
) -> Result<()> {
    // Process the data depending on the presence of the index
    let start = start.unwrap_or(0);
    // Prepare the input VCF
    let tbi_path = format!("{xcf_fn}.tbi");
    let csi_path = format!("{xcf_fn}.csi");
    let has_index = Path::exists(Path::new(&tbi_path)) || Path::exists(Path::new(&csi_path));
    let _ = match has_index {
        true => indexed_xcf(xcf_fn, s1, s2, chrom, start, end, None),
        false => readthrough_xcf(xcf_fn, s1, s2, chrom, start, end, None),
    };
    Ok(())
}
