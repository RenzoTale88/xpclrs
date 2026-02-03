/*
This module provides the binary plink-specific I/O functions.
*/
use crate::io::{consolidate_list, get_gt_index, GenoData};
use itertools::izip;
use itertools::MultiUnzip;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

pub struct PlinkData {
    pub genotypes: Vec<Vec<u8>>, // Matrix of genotypes (0,1,2)
    pub snp_ids: Vec<Vec<u8>>,   // SNP IDs from .bim file
    pub n_samples: usize,        // Number of samples
    pub n_snps: usize,           // Number of SNPs
}

// Primary BED mode ==> SNP major, i.e. codify by individual for SNP
fn bed_snp_major(
    bed_file: &mut File,
    kept_samples: &[bool],
    samples_mat_idx: FxHashMap<usize, usize>,
    kept_sites: &[bool],
) -> Vec<Vec<u8>> {
    // Initialize genotype matrix S X I
    let n_kept_snps = kept_sites
        .iter()
        .filter(|&&v| v)
        .collect::<Vec<&bool>>()
        .len();
    let n_kept_samples = samples_mat_idx.len();
    let mut gts = vec![vec![0u8; n_kept_samples]; n_kept_snps];

    // Read genotypes
    let bytes_per_snp = kept_samples.len().div_ceil(4); // Ensure we have a full byte every time
    let mut snp_bytes = vec![0u8; bytes_per_snp];

    // Initiate progress bar
    log::info!("Load genotypes...");
    let mut kept_snp_idx = 0;
    for keep_site in kept_sites.iter() {
        bed_file
            .read_exact(&mut snp_bytes)
            .expect("Cannot read bytes");
        if *keep_site {
            for (sample_idx, kept_sample) in kept_samples.iter().enumerate() {
                if !*kept_sample {
                    continue;
                }
                let sample_idx_in_mat = samples_mat_idx.get(&sample_idx).expect("Missing index");
                let byte_index = sample_idx / 4; // Which byte (0, 1, 2, ...)
                let shift = 2 * (sample_idx % 4); // Shift to get the right pair
                let bits = (snp_bytes[byte_index] >> shift) & 0b11;

                gts[kept_snp_idx][*sample_idx_in_mat] = bits;
            }
            kept_snp_idx += 1;
        }
    }
    log::info!("Genotypes loaded.");
    log::info!("Matrix shape: {} x {}.", gts.len(), gts[0].len());
    gts
}

// Other BED mode ==> individual major, i.e. codify by SNPs for individual
fn bed_ind_major(
    bed_file: &mut File,
    kept_samples: &[bool],
    samples_mat_idx: FxHashMap<usize, usize>,
    kept_sites: &[bool],
) -> Vec<Vec<u8>> {
    // Initialize genotype matrix I X S
    let n_kept_snps = kept_sites
        .iter()
        .filter(|&&v| v)
        .collect::<Vec<&bool>>()
        .len();
    let n_kept_samples = samples_mat_idx.len();
    let mut gts = vec![vec![0u8; n_kept_samples]; n_kept_snps];

    // Read genotypes
    let bytes_per_ind = n_kept_snps.div_ceil(4);
    let mut ind_bytes = vec![0u8; bytes_per_ind];

    // Load n_snps at the same time, representing one individual
    for (sample_idx, kept_sample) in kept_samples.iter().enumerate() {
        bed_file
            .read_exact(&mut ind_bytes)
            .expect("Cannot read bytes");
        if !*kept_sample {
            continue;
        }
        let sample_idx_in_mat = samples_mat_idx.get(&sample_idx).expect("Missing index");
        let mut kept_snp_idx = 0;
        for (snp_idx, keep_site) in kept_sites.iter().enumerate() {
            if *keep_site {
                let byte_index = snp_idx / 4; // Which byte (0, 1, 2, ...)
                let shift = 2 * (snp_idx % 4); // Shift to get the right pair
                let bits = (ind_bytes[byte_index] >> shift) & 0b11;

                gts[kept_snp_idx][*sample_idx_in_mat] = bits;
                kept_snp_idx += 1;
            }
        }
    }
    log::info!("Genotypes loaded.");
    log::info!("Matrix shape: {} x {}.", gts.len(), gts[0].len());
    gts
}

/// Read PLINK BED/BIM/FAM files and return filtered `GenoData`.
///
/// # Examples
///
/// ```ignore
/// let data = xpclrs::plink::read_plink_files("data/plink", &s1, &s2, "1", 0, None, (None, None)).unwrap();
/// ```
pub fn read_plink_files(
    prefix: &str,
    s1: &[String],
    s2: &[String],
    chrom: &str,
    start: u64,
    end: Option<u64>,
    (phased, rrate): (Option<bool>, Option<f64>),
) -> io::Result<GenoData> {
    // Read .fam file to get number of samples and IDs
    let phased = phased.unwrap_or(false);
    if phased {
        log::warn!("Phased data requested, but PLINK input is always unphased. Proceeding with unphased data.");
    }
    let fam_path = format!("{}.fam", prefix);
    let fam_file = File::open(&fam_path)?;
    let fam_reader = BufReader::new(fam_file);
    let sample_ids_owned: Vec<Vec<u8>> = fam_reader
        .lines()
        .map(|line| line.expect("Cannot read line"))
        .filter_map(|line| {
            let fields = line.split_whitespace().collect::<Vec<&str>>();
            if fields.len() >= 2 {
                Some(fields[0..2].join("_").into_bytes())
            } else {
                None
            }
        })
        .collect();
    let sample_ids: Vec<&[u8]> = sample_ids_owned.iter().map(|id| id.as_slice()).collect();
    // Load sample lists as an u8 array
    let s1 = consolidate_list(&sample_ids, s1).expect("Failed to subset sampleA");
    let s2 = consolidate_list(&sample_ids, s2).expect("Failed to subset sampleB");
    // Fetch the indices of each sample in each list
    let i1_o = get_gt_index(&sample_ids, &s1).expect("Failed to get indeces of sampleA");
    let i2_o = get_gt_index(&sample_ids, &s2).expect("Failed to get indeces of sampleB");
    // Create a vector of indexes of samples to retain, and a
    // Second vector of unique samples to keep
    let mut individuals_ixs: Vec<usize> = [&i1_o[..], &i2_o[..]].concat();
    individuals_ixs.sort();
    individuals_ixs.dedup();
    let new_indexes = individuals_ixs
        .iter()
        .enumerate()
        .map(|(e, v)| (*v, e))
        .collect::<FxHashMap<usize, usize>>();
    // Create new group indexes
    let i1 = i1_o
        .iter()
        .map(|&idx| new_indexes.get(&idx).expect("Missing index"))
        .copied()
        .collect::<Vec<usize>>();
    let i2 = i2_o
        .iter()
        .map(|&idx| new_indexes.get(&idx).expect("Missing index"))
        .copied()
        .collect::<Vec<usize>>();
    // Generate a vector of retained samples
    let kept_samples: Vec<bool> = sample_ids
        .iter()
        .enumerate()
        .map(|(idx, _)| individuals_ixs.contains(&idx))
        .collect();

    // Print number of samples
    log::info!("Samples A: {}", i1.len());
    log::info!("Samples B: {}", i2.len());
    log::debug!("Old samples indexes A: {:?}", i1_o);
    log::debug!("New samples indexes A: {:?}", i1);
    log::debug!("Old samples indexes B: {:?}", i2_o);
    log::debug!("New samples indexes B: {:?}", i2);

    // Dies if no samples are retained
    if s1.is_empty() || s2.is_empty() {
        eprintln!("No samples found in the lists.");
        std::process::exit(1);
    }

    // Read .bim file to get number of SNPs and IDs
    let bim_path = format!("{}.bim", prefix);
    let bim_file = File::open(&bim_path)?;
    let bim_reader = BufReader::new(bim_file);
    let mut n_snps: usize = 0;
    let mut positions: Vec<usize> = vec![];
    let mut gd_data: Vec<f64> = vec![];
    let mut keep_vec: Vec<bool> = vec![];
    for line in bim_reader.lines() {
        n_snps += 1;
        let line = line?;
        let fields: Vec<&str> = line.split_whitespace().collect();
        if chrom == fields[0] {
            let pos: u64 = fields[3].parse().unwrap();
            if pos >= start && (end.is_none() || pos <= end.unwrap()) {
                positions.push(fields[3].parse::<usize>().unwrap());
                // Add genetic distance in the bim file, if provided, otherwise compute it
                // based on the recombination rate and physical position
                let gd = match fields[2].parse::<f64>() {
                    Ok(v) => {
                        if v != 0.0 {
                            v
                        } else {
                            pos as f64 * rrate.unwrap_or(1e-8)
                        }
                    }
                    Err(_) => pos as f64 * rrate.unwrap_or(1e-8),
                };
                gd_data.push(gd);
                keep_vec.push(true);
            } else {
                keep_vec.push(false);
            }
        } else {
            keep_vec.push(false);
        }
    }
    log::info!("Number of SNPs: {}", n_snps);

    // Read .bed file
    let bed_path = format!("{}.bed", prefix);
    let mut bed_file = File::open(&bed_path)?;

    // Check magic numbers
    let mut magic = [0u8; 3];
    bed_file.read_exact(&mut magic)?;
    if magic != [0x6c, 0x1b, 0x01] && magic != [0x6c, 0x1b, 0x00] {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid BED file magic numbers",
        ));
    }

    // Define if it is SNP major or individual major
    let genotypes: Vec<Vec<u8>> = if magic[2] == 0x01 {
        log::info!("Reading in SNP-major mode.");
        bed_snp_major(&mut bed_file, &kept_samples, new_indexes, &keep_vec)
    } else {
        log::info!("Reading in Individual-major mode.");
        bed_ind_major(&mut bed_file, &kept_samples, new_indexes, &keep_vec)
    };

    // Filter the dataset
    let mut monom_gt2 = 0;
    let mut miss_gt1 = 0;
    let mut miss_gt2 = 0;
    let mut pass = 0;
    let mut skipped = 0;
    let mut tot = 0;
    let (gt1_data, gt2_data, positions, gd_data): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
        izip!(genotypes.iter(), positions.iter(), gd_data.iter())
            .filter_map(|(snp, position, gd)| {
                // Prepare GT1 and GT2
                // Usually, 00, 10 and 11 mean 2, 1 and 0 minor allele counts respectively
                let gt1 = i1
                    .iter()
                    .map(|&i| match snp[i] {
                        0b00 => 2_i8,
                        0b10 => 1_i8,
                        0b11 => 0_i8,
                        _ => -9_i8,
                    })
                    .collect::<Vec<i8>>();
                let gt2 = i2
                    .iter()
                    .map(|&i| match snp[i] {
                        0b00 => 2_i8,
                        0b10 => 1_i8,
                        0b11 => 0_i8,
                        _ => -9_i8,
                    })
                    .collect::<Vec<i8>>();
                // Define if all variants are missing in GT1 or GT2
                let all_missing_g1 = gt1.iter().all(|i| *i == -9_i8);
                let all_missing_g2 = gt2.iter().all(|i| *i == -9_i8);
                // Get non-missing genotypes in GT2
                let non_missing = gt2
                    .iter()
                    .filter_map(|&i| if i != -9_i8 { Some(i as i64) } else { None })
                    .collect::<Vec<i64>>();
                // Count how many minor alleles are present
                let n_minor_alleles: i64 = non_missing.iter().sum();
                let n_major_alleles: i64 = (non_missing.iter().len() as i64 * 2) - n_minor_alleles;
                // If the minor allele count is <=1 or >=(n_samples * 2 -1),
                // the variants are monomorphic or singleton (only 1 major or minor)
                // alleles; discard these.
                let monom_or_singleton_gt2 = n_minor_alleles <= 1 || n_major_alleles <= 1;
                if all_missing_g1 || all_missing_g2 {
                    skipped += 1;
                    tot += 1;
                    miss_gt1 += if all_missing_g1 { 1 } else { 0 };
                    miss_gt2 += if all_missing_g2 { 1 } else { 0 };
                    None
                } else if monom_or_singleton_gt2 {
                    skipped += 1;
                    tot += 1;
                    monom_gt2 += 1;
                    None
                } else {
                    pass += 1;
                    tot += 1;
                    Some((gt1, gt2, *position, *gd))
                }
            })
            .multiunzip();

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
    log::info!(" - 0 multiallelic");
    log::info!(" - {monom_gt2} monomorphic/singleton in pop B");
    log::info!(" - {miss_gt1} all-missing in pop A");
    log::info!(" - {miss_gt2} all-missing in pop B");

    // Return the data structure
    Ok(GenoData {
        positions,
        gt1: gt1_data,
        gt2: gt2_data,
        gdistances: gd_data,
    })
}
