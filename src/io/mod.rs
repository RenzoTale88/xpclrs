/*
This module provides the shared I/O functions (e.g. VCF/BCF readers, text readers etc).
*/
use crate::methods::XPCLRResult;
use crate::plink::read_plink_files;
use crate::xcf::{indexed_xcf, readthrough_xcf};
use anyhow::Result;
use flate2::write;
use flate2::Compression;
use rust_htslib::bcf::record::Genotype;
use statistical::{mean, population_standard_deviation};
use std::{
    collections::HashSet,
    ffi::OsStr,
    fmt::Display,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

// Genotype data structure
pub struct GenoData {
    pub positions: Vec<usize>,
    // Per-site, per-sample alternate allele counts encoded as i8:
    // -9 = missing, 0/1/2 = alt allele count for diploids
    pub gt1: Vec<Vec<i8>>,
    pub gt2: Vec<Vec<i8>>,
    pub gdistances: Vec<f64>,
}

// Genotypes to counts (optimized single-pass implementation)
pub fn gt2gcount(gt: Genotype, ref_ix: u32) -> i8 {
    // Single-pass count without allocating a temporary allele vector.
    let mut seen = false;
    let mut alt_count: i8 = 0;
    for allele in gt.iter().filter_map(|a| a.index()) {
        seen = true;
        if allele != ref_ix {
            alt_count += 1;
        }
    }

    if !seen {
        -9_i8
    } else if alt_count == 0 {
        0
    } else {
        alt_count
    }
}

/// Read file, either compressed or uncompressed
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
pub fn consolidate_list(full_list: &Vec<&[u8]>, subset: &[String]) -> Result<Vec<String>> {
    let filtered = subset
        .iter()
        .filter(|&s| full_list.contains(&s.as_bytes()))
        .cloned()
        .collect::<Vec<String>>();
    Ok(filtered)
}

// Function to get the index of each sample in the lists
pub fn get_gt_index(full_list: &Vec<&[u8]>, subset: &[String]) -> Result<Vec<usize>> {
    let subset_set: HashSet<_> = subset.iter().collect();
    let indices: Vec<usize> = full_list
        .iter()
        .enumerate()
        .filter_map(|(idx, &val)| {
            if subset_set.contains(&String::from_utf8(val.to_owned()).unwrap()) {
                Some(idx)
            } else {
                None
            }
        })
        .collect();
    Ok(indices)
}

// Load the genotypes for the given samples
pub fn process_xcf(
    xcf_fn: String,
    s1: &[String],
    s2: &[String],
    chrom: &str,
    start: Option<u64>,
    end: Option<u64>,
    (phased, rrate, gdistkey, n_threads): (Option<bool>, Option<f64>, Option<String>, usize),
) -> Result<GenoData> {
    // Process the data depending on the presence of the index
    let start = start.unwrap_or(0);
    // Prepare the input VCF
    let tbi_path = format!("{xcf_fn}.tbi");
    let csi_path = format!("{xcf_fn}.csi");
    let has_index = Path::exists(Path::new(&tbi_path)) || Path::exists(Path::new(&csi_path));
    let g_data = match has_index {
        true => indexed_xcf(
            xcf_fn,
            s1,
            s2,
            chrom,
            start,
            end,
            (phased, rrate, gdistkey, n_threads),
        ),
        false => readthrough_xcf(
            xcf_fn,
            s1,
            s2,
            chrom,
            start,
            end,
            (phased, rrate, gdistkey, n_threads),
        ),
    }
    .expect("Failed to parse the VCF/BCF file");
    Ok(g_data)
}

// Load the genotypes for the given samples
pub fn process_plink(
    plink_root: String,
    s1: &[String],
    s2: &[String],
    chrom: &str,
    start: Option<u64>,
    end: Option<u64>,
    (phased, rrate): (Option<bool>, Option<f64>),
) -> Result<GenoData> {
    // Process the data depending on the presence of the index
    let start = start.unwrap_or(0);
    // Prepare the input VCF
    let g_data = read_plink_files(&plink_root, s1, s2, chrom, start, end, (phased, rrate))
        .expect("Failed to parse the BED/BIM/FAM file");
    Ok(g_data)
}

pub fn write_table(filename: &str) -> Box<dyn Write> {
    let path = Path::new(filename);
    let file = File::create(path)
        .unwrap_or_else(|why| panic!("couldn't open {}: {}", path.display(), why));

    let writer: Box<dyn Write> = if path.extension() == Some(OsStr::new("gz")) {
        // Create a GzEncoder which compresses data and writes it to the file
        let gz_encoder = write::GzEncoder::new(file, Compression::default());
        // Wrap the GzEncoder in a BufWriter for efficient buffering
        Box::new(BufWriter::with_capacity(128 * 1024, gz_encoder))
    } else {
        // Wrap the file in a BufWriter directly for uncompressed writing
        Box::new(BufWriter::with_capacity(128 * 1024, file))
    };

    // Create the fastq::Writer using the boxed writer
    writer
}

// The following results
// n, (start, stop, bpi, bpe, nsnps, avail), (model_li, null_li, selectionc)
// Map to:
// win index, (start and stop of window), (bpi and bpe are edges),
pub fn to_table(
    chrom: &str,
    xpclr_res: &[(usize, XPCLRResult)],
    xpclr_tsv: &mut Box<dyn std::io::Write>,
    outfmt: &str,
) -> Result<()> {
    let delim = match outfmt {
        "tsv" => "\t",
        "txt" => " ",
        "csv" => ",",
        _ => "\t",
    };
    // Write header
    writeln!(xpclr_tsv, "chrom{delim}start{delim}stop{delim}pos_start{delim}pos_stop{delim}modelL{delim}nullL{delim}sel_coef{delim}nSNPs{delim}nSNPs_avail{delim}xpclr{delim}xpclr_norm")?;

    // Compute normalizing factors
    let xpclr_values = xpclr_res
        .iter()
        .filter_map(|(_, r)| {
            if r.xpclr.is_nan() {
                None
            } else {
                Some(r.xpclr)
            }
        })
        .collect::<Vec<f64>>();
    let mean_xpclr = mean(&xpclr_values);
    let std_xpclr = population_standard_deviation(&xpclr_values, None);
    log::info!(
        "XP-CLR mean +/- st.d: {mean_xpclr} +/- {std_xpclr} (N={})",
        xpclr_values.len()
    );

    for (_n, res) in xpclr_res {
        let start = res.window.0;
        let stop = res.window.1;
        let bpi = res.window.2;
        let bpe = res.window.3;
        let nsnps = res.window.4;
        let avail = res.window.5;
        let model_li = res.ll_sel;
        let null_li = res.ll_neut;
        let selectionc = res.sel_coeff;
        let xpclr = res.xpclr;
        let xpclr_normalized = (res.xpclr - mean_xpclr) / std_xpclr;
        writeln!(xpclr_tsv, "{chrom}{delim}{start}{delim}{stop}{delim}{bpi}{delim}{bpe}{delim}{model_li}{delim}{null_li}{delim}{selectionc}{delim}{nsnps}{delim}{avail}{delim}{xpclr}{delim}{xpclr_normalized}")?;
    }

    Ok(())
}
