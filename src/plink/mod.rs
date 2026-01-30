/*
This module provides the binary plink-specific I/O functions.
*/
use crate::io::{consolidate_list, GenoData, get_gt_index};
use itertools::izip;
use itertools::MultiUnzip;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};


pub struct PlinkData {
    pub genotypes: Vec<Vec<u8>>, // Matrix of genotypes (0,1,2)
    pub snp_ids: Vec<Vec<u8>>,   // SNP IDs from .bim file
    pub n_samples: usize,        // Number of samples
    pub n_snps: usize,           // Number of SNPs
}

// Primary BED mode ==> SNP major, i.e. codify by individual for SNP
fn bed_snp_major(bed_file: &mut File, n_snps: usize, n_samples: usize) -> Vec<Vec<u8>> {
    // Initialize genotype matrix S X I
    let mut gts = vec![vec![0u8; n_samples]; n_snps];

    // Read genotypes
    let bytes_per_snp = n_samples.div_ceil(4); // Ensure we have a full byte every time
    let mut snp_bytes = vec![0u8; bytes_per_snp];

    // Initiate progress bar
    log::info!("Load genotypes...");

    for snp_row in gts.iter_mut() {
        bed_file
            .read_exact(&mut snp_bytes)
            .expect("Cannot read bytes");
        for sample_idx in 0..n_samples {
            let byte_index = sample_idx / 4; // Which byte (0, 1, 2, ...)
            let shift = 2 * (sample_idx % 4); // Shift to get the right pair
            let bits = (snp_bytes[byte_index] >> shift) & 0b11;

            snp_row[sample_idx] = bits;
        }
    }
    log::info!("Genotypes loaded.");
    gts
}

// Other BED mode ==> individual major, i.e. codify by SNPs for individual
fn bed_ind_major(bed_file: &mut File, n_snps: usize, n_samples: usize) -> Vec<Vec<u8>> {
    // Initialize genotype matrix I X S
    let mut gts = vec![vec![0u8; n_samples]; n_snps];

    // Read genotypes
    let bytes_per_ind = n_snps.div_ceil(4);
    let mut ind_bytes = vec![0u8; bytes_per_ind];

    // Load individual genotypes
    for sample_idx in 0..n_samples {
        bed_file
            .read_exact(&mut ind_bytes)
            .expect("Cannot read bytes");

        for snp_idx in 0..n_snps {
            let byte_index = snp_idx / 4; // Which byte (0, 1, 2, ...)
            let shift = 2 * (snp_idx % 4); // Shift to get the right pair
            let bits = (ind_bytes[byte_index] >> shift) & 0b11;

            gts[snp_idx][sample_idx] = bits;
        }
    }
    log::info!("Genotypes loaded.");
    gts
}

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
    let sample_ids: Vec<&[u8]> = sample_ids_owned
        .iter()
        .map(|id| id.as_slice())
        .collect();
    // Number of samples for later use.
    let n_samples = sample_ids.len();
    // Load sample lists as an u8 array
    let s1 = consolidate_list(&sample_ids, s1).expect("Failed to subset sampleA");
    let s2 = consolidate_list(&sample_ids, s2).expect("Failed to subset sampleB");
    // Fetch the indices of each sample in each list
    let i1 = get_gt_index(&sample_ids, &s1).expect("Failed to get indeces of sampleA");
    let i2 = get_gt_index(&sample_ids, &s2).expect("Failed to get indeces of sampleB");

    log::info!("Number of samples: {}", n_samples);

    // Read .bim file to get number of SNPs and IDs
    let bim_path = format!("{}.bim", prefix);
    let bim_file = File::open(&bim_path)?;
    let bim_reader = BufReader::new(bim_file);
    let mut n_snps: usize = 0;
    let mut positions: Vec<usize> = vec![];
    let mut keep_vec: Vec<bool> = vec![];
    for line in bim_reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split_whitespace().collect();
        if chrom == fields[0] {
            let pos: u64 = fields[3].parse().unwrap();
            if pos >= start && (end.is_none() || pos <= end.unwrap()) {
                positions.push(fields[3].parse::<usize>().unwrap());
                keep_vec.push(true);
                n_snps += 1;
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
        bed_snp_major(&mut bed_file, n_snps, n_samples)
    } else {
        log::info!("Reading in Individual-major mode.");
        bed_ind_major(&mut bed_file, n_snps, n_samples)
    }
        .into_iter()
        .zip(keep_vec.into_iter())
        .filter_map(|(gt_row, keep)| if keep { Some(gt_row) } else { None })
        .collect();

    // Filter the dataset
    let mut monom_gt2 = 0;
    let mut miss_gt1 = 0;
    let mut miss_gt2 = 0;
    let mut pass = 0;
    let mut skipped = 0;
    let mut tot = 0;
    let (gt1_data, gt2_data, positions, gd_data): (Vec<Vec<i8>>, Vec<Vec<i8>>, Vec<usize>, Vec<f64>) = izip!(
            genotypes.iter(),
            positions.iter()
        ).filter_map(|(snp, position)| {
            // Compute genetic distance if requested
            let gd = *position as f64 * rrate.unwrap_or(1e-8);
            // Prepare GT1 and GT2
            let gt1 = i1
                .iter()
                .map(|&i| match snp[i] {
                        0b00 => 0_i8,
                        0b10 => 1_i8,
                        0b11 => 2_i8,
                        _ => -9_i8,
                    }
                )
                .collect::<Vec<i8>>();
            let gt2 = i2
                .iter()
                .map(|&i| match snp[i] {
                        0b00 => 0_i8,
                        0b10 => 1_i8,
                        0b11 => 2_i8,
                        _ => -9_i8,
                    }
                )
                .collect::<Vec<i8>>();
            let all_missing_g1 = gt1.iter().all(|i| *i == -9_i8);
            let all_missing_g2 = gt2.iter().all(|i| *i == -9_i8);
            if all_missing_g1 || all_missing_g2 {
                skipped += 1;
                tot += 1;
                miss_gt1 += if all_missing_g1 { 1 } else { 0 };
                miss_gt2 += if all_missing_g2 { 1 } else { 0 };
                None
            } else if gt2.iter().filter(|&&i| i != -9_i8).all(|i| *i == 0_i8) || gt2.iter().filter(|&&i| i != -9_i8).all(|i| *i == 2_i8) {
                skipped += 1;
                tot += 1;
                monom_gt2 += 1;
                None
            } else {
                pass += 1;
                tot += 1;
                Some((gt1, gt2, *position, gd))
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
    log::info!(" - {monom_gt2} monomorphic in pop B");
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
