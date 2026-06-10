use anyhow::Result;
use clap::{builder::PossibleValue, Arg, Command};
use env_logger::{self, Env};
use rustc_hash::FxHashMap;
use statistical::{mean, population_standard_deviation};
use std::fs::read_dir;
use xpclrs::{io::{read_file, table_writer}};

/*
 --format FORMAT, -F FORMAT
                       input expected. One of "vcf" (default), "hdf5", "zarr" or "txt"
 --map MAP             If using XPCLR-style text format. Input map file as per XPCLR specs (tab separated)
 --popA POPA           If using XPCLR-style text format. Filepath to population A genotypes (space separated)
 --popB POPB           If using XPCLR-style text format. Filepath to population B genotypes (space separated)
 --verbose VERBOSE, -V VERBOSE
                       How verbose to be in logging. 10=DEBUG, 20=INFO, 30=WARN, 40=ERROR, 50=CRITICAL
*/


// Define a local, chromosome-aware XPCLR result structure
pub struct XPCLRResult {
    pub chrom_idx: usize,
    pub window: (usize, usize, usize, usize, usize, usize), // (start, stop, bpi, bpe, nsnps, navail)
    pub ll_sel: f64,
    pub ll_neut: f64,
    pub sel_coeff: f64,
    pub xpclr: f64,
}

/// Write XP-CLR results to a delimited table.
/// The following results
/// n, (start, stop, bpi, bpe, nsnps, avail), (model_li, null_li, selectionc)
/// Map to:
/// win index, (start and stop of window), (bpi and bpe are edges),
///
/// # Examples
///
/// ```ignore
/// xpclrs::io::to_table("1", &results, &mut writer, "tsv").unwrap();
/// ```
pub fn to_table(
    chrom_hash: &FxHashMap<usize, String>,
    xpclr_res: &[XPCLRResult],
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
        .filter_map(|r| {
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

    for res in xpclr_res {
        let chrom = match chrom_hash.get(&res.chrom_idx){
            Some(c) => c,
            None => &String::from("0"),
        };
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


fn main() {
    let version = env!("CARGO_PKG_VERSION");
    let matches = Command::new("combine_xpclrs")
        .version(version)
        .author("Andrea Talenti <andrea.talenti@ed.ac.uk>")
        .about("Combine and normalize the XP-CLR outputs\n")
        .arg(
            Arg::new("INPUT")
                .short('I')
                .long("input")
                .required(true)
                .help("input directory"),
        )
        .arg(
            Arg::new("INFMT")
                .short('f')
                .long("informat")
                .required(false)
                .default_value("tsv")
                .value_parser([
                    PossibleValue::new("tsv"),
                    PossibleValue::new("txt"),
                    PossibleValue::new("csv"),
                ])
                .help("Format to save the output (csv, tsv, txt)"),
        )
        .arg(
            Arg::new("OUT")
                .short('O')
                .long("out")
                .required(true)
                .help("Output file name."),
        )
        .arg(
            Arg::new("OUTFMT")
                .short('F')
                .long("outformat")
                .required(false)
                .default_value("tsv")
                .value_parser([
                    PossibleValue::new("tsv"),
                    PossibleValue::new("txt"),
                    PossibleValue::new("csv"),
                ])
                .help("Format to save the output (csv, tsv, txt)"),
        )
        .arg(
            Arg::new("LOG")
                .short('l')
                .long("log")
                .required(false)
                .default_value("info")
                .value_parser([PossibleValue::new("info"), PossibleValue::new("debug")])
                .help("Logging level."),
        )
        .get_matches();

    // Fixed parameters
    // Get the input file
    let input_path = matches
        .get_one::<String>("INPUT")
        .expect("Input file is required")
        .to_owned();
    let in_fmt = matches
        .get_one::<String>("INFMT")
        .expect("Invalid input format")
        .to_owned();
    let delim = match in_fmt.as_str() {
        "tsv" => "\t",
        "txt" => " ",
        "csv" => ",",
        _ => "\t",
    };

    // Get the output path
    let out_path = matches
        .get_one::<String>("OUT")
        .expect("Output file is required")
        .to_owned();
    let out_fmt = matches
        .get_one::<String>("OUTFMT")
        .expect("Invalid output format")
        .to_owned();

    // set up logging
    let log_level = matches
        .get_one::<String>("LOG")
        .expect("Log level not valid")
        .to_owned();
    env_logger::Builder::from_env(Env::default().default_filter_or(log_level)).init();

    // Initial logging
    log::info!("combiner v{version}");

    // Import each file in the input directory and combine them into a single table
    let mut xpclr_res = Vec::new();
    // Define a chromosome index map
    let mut chromosome_indexes: FxHashMap<String, usize> = FxHashMap::default();
    // Process each file in the directory
    for entry in read_dir(input_path).expect("Failed to read input directory") {
        let entry = entry.expect("Failed to read directory entry");
        let path = entry.path();
        if path.is_file() {
            log::info!("Reading file: {}", path.display());
            // Process each record in the dataset
            let infile = path.to_str().expect("Cannot convert to string");
            let records = read_file(infile).expect("Cannot read file");
            for (line_n, record) in records.enumerate() {
                if line_n == 0 {
                    // Skip header
                    continue;
                }
                let record = record.expect("Failed to read record");
                let record = record.split(delim).collect::<Vec<&str>>();
                let chrom = record[0].to_string();
                let start = record[1].parse::<usize>().unwrap();
                let end = record[2].parse::<usize>().unwrap();
                let pos_start = record[3].parse::<usize>().unwrap();
                let pos_stop = record[4].parse::<usize>().unwrap();
                let nsnps = record[8].parse::<usize>().unwrap();
                let nsnps_avail = record[9].parse::<usize>().unwrap();
                let window = (start, end, pos_start, pos_stop, nsnps, nsnps_avail);
                let ll_sel = record[5].parse::<f64>().unwrap_or(f64::NAN);
                let ll_neut = record[6].parse::<f64>().unwrap_or(f64::NAN);
                let sel_coeff = record[7].parse::<f64>().unwrap_or(f64::NAN);
                let xpclr = record[10].parse::<f64>().unwrap_or(f64::NAN);
                let chrom_idx = if chromosome_indexes.is_empty(){
                    chromosome_indexes.insert(chrom.clone(), 0);
                    0
                } else {
                    if chromosome_indexes.contains_key(&chrom) {
                        chromosome_indexes.get(&chrom).copied().unwrap()
                    } else {
                        let idx = chromosome_indexes.len();
                        chromosome_indexes.insert(chrom.clone(), idx);
                        idx
                    }
                };
                let clrs_out = XPCLRResult{
                    chrom_idx,
                    window,
                    ll_sel,
                    ll_neut,
                    sel_coeff,
                    xpclr,
                };
                xpclr_res.push(clrs_out);
            }
        }
    }
    // Prepare the chromosome index conversion map
    let inverted_indexes = chromosome_indexes.iter().map(|(k, v)| (*v, k.clone())).collect::<FxHashMap<usize, String>>();
    // Save the output
    let mut xpclr_tsv = table_writer(&out_path);
    let _ = to_table(&inverted_indexes, &xpclr_res, &mut xpclr_tsv, &out_fmt);

    log::info!("Collected and normalized results.")
}
