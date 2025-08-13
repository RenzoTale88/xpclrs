/*
This module provides the I/O functions (e.g. VCF/BCF readers, text readers etc).
*/
use anyhow::Result;
use counter::Counter;
use flate2::write;
use flate2::Compression;
use itertools::MultiUnzip;
use rust_htslib::{
    bcf::{
        self,
        record::{Genotype, GenotypeAllele},
        IndexedReader, Read, Reader,
    },
};
use statistical::{mean, standard_deviation};
use std::{
    collections::HashSet,
    ffi::OsStr,
    fmt::Display,
    fs::File,
    io::{
        BufRead,
        BufReader,
        BufWriter,
        Write
    },
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

// Process indexed XCF file
fn indexed_xcf(
    xcf_fn: String,
    s1: &[String],
    s2: &[String],
    chrom: &str,
    start: u64,
    end: Option<u64>,
    (rrate, _gdistkey): (Option<f32>, Option<String>),
) -> Result<(Vec<usize>, Vec<Vec<Genotype>>, Vec<Vec<Genotype>>, Vec<f32>)> {
    log::info!("Indexed reader.");
    // Prepare the indexed reader
    let mut reader = IndexedReader::from_path(xcf_fn).expect("Cannot load indexed BCF/VCF file");
    let rrate = rrate.unwrap_or(1e-8);

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
    // let mut monom_gt1 = 0;
    let mut monom_gt2 = 0;
    let mut miss_gt1 = 0;
    let mut miss_gt2 = 0;
    let mut pass = 0;
    let mut skipped = 0;
    let mut tot = 0;
    let (positions, gt1_data, gt2_data, gd_data): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) = reader
        .records()
        .filter_map(|r| {
            let record = r.ok()?;
            tot += 1;
            let genotypes = record.genotypes().expect("Cannot fetch the genotypes");
            let gt1 = i1
                .iter()
                .map(|i| genotypes.get(*i))
                .collect::<Vec<Genotype>>();
            let gt2 = i2
                .iter()
                .map(|i| genotypes.get(*i))
                .collect::<Vec<Genotype>>();
            // Check that the site is biallelic
            let alleles1: Counter<u32> = gt1
                .iter()
                .flat_map(|g| g.iter().filter_map(|&a: &GenotypeAllele| a.index()))
                .collect::<Counter<u32>>();

            let alleles2: Counter<u32> = gt2
                .iter()
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
                Some((record.pos() as usize, gt1, gt2, record.pos() as f32 * rrate))
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
    log::info!(" - {multiallelic} multiallelic");
    log::info!(" - {monom_gt2} monomorphic in pop B");
    log::info!(" - {miss_gt1} all-missing in pop A");
    log::info!(" - {miss_gt2} all-missing in pop B");
    Ok((positions, gt1_data, gt2_data, gd_data))
}

// Process unindexed XCF
fn readthrough_xcf(
    xcf_fn: String,
    s1: &[String],
    s2: &[String],
    chrom: &str,
    start: u64,
    end: Option<u64>,
    (rrate, _gdistkey): (Option<f32>, Option<String>),
) -> Result<(Vec<usize>, Vec<Vec<Genotype>>, Vec<Vec<Genotype>>, Vec<f32>)> {
    log::info!("Streamed reader.");
    log::info!("This is substantially slower than the indexed one.");
    log::info!("Consider generating an index for your BCF/VCF file.");
    // Prepare the indexed reader
    let mut reader = Reader::from_path(xcf_fn).expect("Cannot load indexed BCF/VCF file");
    let end = end.unwrap_or(999999999);
    let rrate = rrate.unwrap_or(1e-8);

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
    // let mut monom_gt1 = 0;
    let mut monom_gt2 = 0;
    let mut miss_gt1 = 0;
    let mut miss_gt2 = 0;
    let mut pass = 0;
    let mut skipped = 0;
    let mut tot = 0;
    let (positions, gt1_data, gt2_data, gd_data): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) = reader
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
            let gt1 = i1
                .iter()
                .map(|i| genotypes.get(*i))
                .collect::<Vec<Genotype>>();
            let gt2 = i2
                .iter()
                .map(|i| genotypes.get(*i))
                .collect::<Vec<Genotype>>();

            // Check that the site is biallelic
            let alleles1: Counter<u32> = gt1
                .iter()
                .flat_map(|g| g.iter().filter_map(|&a: &GenotypeAllele| a.index()))
                .collect::<Counter<u32>>();

            let alleles2: Counter<u32> = gt2
                .iter()
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
                Some((record.pos() as usize, gt1, gt2, record.pos() as f32 * rrate))
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
    log::info!(" - {multiallelic} multiallelic");
    log::info!(" - {monom_gt2} monomorphic in pop B");
    log::info!(" - {miss_gt1} all-missing in pop A");
    log::info!(" - {miss_gt2} all-missing in pop B");
    Ok((positions, gt1_data, gt2_data, gd_data))
}

// Load the genotypes for the given samples
pub fn process_xcf(
    xcf_fn: String,
    s1: &[String],
    s2: &[String],
    chrom: &str,
    start: Option<u64>,
    end: Option<u64>,
    (rrate, _gdistkey): (Option<f32>, Option<String>),
) -> Result<(Vec<usize>, Vec<Vec<Genotype>>, Vec<Vec<Genotype>>, Vec<f32>)> {
    // Process the data depending on the presence of the index
    let start = start.unwrap_or(0);
    // Prepare the input VCF
    let tbi_path = format!("{xcf_fn}.tbi");
    let csi_path = format!("{xcf_fn}.csi");
    let has_index = Path::exists(Path::new(&tbi_path)) || Path::exists(Path::new(&csi_path));
    let (positions, gt1_data, gt2_data, gd_data) = match has_index {
        true => indexed_xcf(xcf_fn, s1, s2, chrom, start, end, (rrate, None)),
        false => readthrough_xcf(xcf_fn, s1, s2, chrom, start, end, (rrate, None)),
    }
    .expect("Failed to parse the VCF/BCF file");
    Ok((positions, gt1_data, gt2_data, gd_data))
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
    xpclr_res: Vec<(usize, (usize, usize, usize, usize, usize, usize), (f32, f32, f32, f32))>,
    xpclr_tsv: &mut Box<dyn std::io::Write>
) -> Result<()> {
    // Write header
    writeln!(xpclr_tsv, "chrom,start,stop,pos_start,pos_stop,modelL,nullL,sel_coef,nSNPs,nSNPs_avail,xpclr,xpclr_norm")?;

    // Compute normalizing factors
    let xpclr_values = xpclr_res
        .iter()
        .filter_map(|&v| if v.2.3.is_nan() {
            None
        } else {
            Some(v.2.3)
        })
        .collect::<Vec<f32>>();
    let mean_xpclr = mean( &xpclr_values );
    let std_xpclr = standard_deviation( &xpclr_values, None );
    log::info!("XP-CLR mean +/- st.d: {mean_xpclr} (+/-{std_xpclr})");

    for (_n, (start, stop, bpi, bpe, nsnps, avail), (model_li, null_li, selectionc, xpclr)) in xpclr_res {
        let xpclr_normalized = (xpclr - mean_xpclr) / std_xpclr;
        writeln!(xpclr_tsv, "{chrom},{start},{stop},{bpi},{bpe},{model_li},{null_li},{selectionc},{nsnps},{avail},{xpclr},{xpclr_normalized}")?;
    }

    Ok(())
}
