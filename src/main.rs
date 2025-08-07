use clap::{value_parser, Arg, Command};
use rust_htslib::bcf::*;
use std::path::Path;
use xpclr::{
    io::{read_file, read_xcf, XcfReader},
    methods::{compute_harding_likelihood},
};

/*
 --format FORMAT, -F FORMAT
                       input expected. One of "vcf" (default), "hdf5", "zarr" or "txt"
 --map MAP             If using XPCLR-style text format. Input map file as per XPCLR specs (tab separated)
 --popA POPA           If using XPCLR-style text format. Filepath to population A genotypes (space separated)
 --popB POPB           If using XPCLR-style text format. Filepath to population B genotypes (space separated)
 --verbose VERBOSE, -V VERBOSE
                       How verbose to be in logging. 10=DEBUG, 20=INFO, 30=WARN, 40=ERROR, 50=CRITICAL
*/

fn main() {
    let version = env!("CARGO_PKG_VERSION");
    let matches = Command::new("xpclr")
        .version(version)
        .author("Andrea Talenti <andrea.talenti@ed.ac.uk>")
        .about("Compute the XP-CLR for a pair of populations from a VCF file.")
        .arg(
            Arg::new("VCF")
                .short('I')
                .long("input")
                .required(true)
                .help("input VCF file"),
        )
        .arg(
            Arg::new("OUT")
                .short('O')
                .long("out")
                .required(true)
                .help("Output file name."),
        )
        .arg(
            Arg::new("SAMPLES_A")
                .short('A')
                .long("samplesA")
                .required(true)
                .help("Samples in population A. Path to file with each ID on a line."),
        )
        .arg(
            Arg::new("SAMPLES_B")
                .short('B')
                .long("samplesB")
                .required(true)
                .help("Samples in population B. Path to file with each ID on a line."),
        )
        .arg(
            Arg::new("RECRATE")
                .long("rrate")
                .short('R')
                .required(false)
                .default_value("1e-8")
                .value_parser(value_parser!(f32))
                .help("Recombination rate per base."),
        )
        .arg(
            Arg::new("LDCUTOFF")
                .long("ld")
                .short('L')
                .required(false)
                .default_value("0.95")
                .value_parser(value_parser!(f32))
                .help("LD cutoff."),
        )
        .arg(
            Arg::new("MAXSNPS")
                .long("maxsnps")
                .short('m')
                .required(false)
                .default_value("200")
                .value_parser(value_parser!(i32))
                .help("Max SNPs in a window."),
        )
        .arg(
            Arg::new("MINSNPS")
                .long("minsnps")
                .short('N')
                .required(false)
                .default_value("10")
                .value_parser(value_parser!(i32))
                .help("Min SNPs in a window."),
        )
        .arg(
            Arg::new("SIZE")
                .long("size")
                .required(false)
                .default_value("20000")
                .value_parser(value_parser!(i32))
                .help("Sliding window size."),
        )
        .arg(
            Arg::new("START")
                .long("start")
                .required(false)
                .default_value("1")
                .value_parser(value_parser!(i32))
                .help("Start position for the sliding windows."),
        )
        .arg(
            Arg::new("STOP")
                .long("stop")
                .required(false)
                .help("Stop position for the sliding windows."),
        )
        .arg(
            Arg::new("STEP")
                .long("step")
                .required(false)
                .default_value("20000")
                .value_parser(value_parser!(i32))
                .help("Step size for the sliding windows."),
        )
        .arg(
            Arg::new("CHROM")
                .short('C')
                .long("chr")
                .required(false)
                .help("Chromosome to analyse."),
        )
        .arg(
            Arg::new("DISTKEYS")
                .long("gdistkey")
                .required(false)
                .help("Key for genetic position in variants table of hdf5/VCF"),
        )
        .arg(
            Arg::new("NTHREADS")
                .short('t')
                .long("threads")
                .required(false)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Number of threads to use"),
        )
        .get_matches();

    // Get the VCF
    let xcf_path = matches
        .get_one::<String>("VCF")
        .expect("VCF file is required");
    let tbi_path = format!("{xcf_path}.tbi");
    let csi_path = format!("{xcf_path}.csi");
    let has_index = Path::exists(Path::new(&tbi_path)) || Path::exists(Path::new(&csi_path));
    let input_xcf = read_xcf(xcf_path, has_index).expect("Failed to read VCF file");
    let xcf_header = match &input_xcf {
        XcfReader::Indexed(reader) => reader.header(),
        XcfReader::Readthrough(reader) => reader.header(),
    };
    let sample_list = xcf_header.samples();
    println!("Samples in VCF: {}", xcf_header.sample_count());

    // Get the output path
    let _out_path = matches
        .get_one::<String>("OUT")
        .expect("Output file is required");
    let n_threads = matches
        .get_one::<usize>("NTHREADS")
        .expect("Number of threads invalid");

    // Get the sample lists
    let sample_a = matches
        .get_one::<String>("SAMPLES_A")
        .expect("Samples A file is required");
    let sample_b = matches
        .get_one::<String>("SAMPLES_B")
        .expect("Samples B file is required");

    // Load sample lists as an u8 array
    let samples_a = read_file(sample_a)
        .expect("Invalid file for samples A")
        .map(|s| s.expect("Failed to read sample A line"))
        .filter(|s| sample_list.contains(&s.as_bytes()))
        .collect::<Vec<String>>();
    let samples_b = read_file(sample_b)
        .expect("Invalid file for samples B")
        .map(|s| s.expect("Failed to read sample B line"))
        .filter(|s| sample_list.contains(&s.as_bytes()))
        .collect::<Vec<String>>();

    // Print number of samples
    println!("Samples A: {}", samples_a.len());
    println!("Samples B: {}", samples_b.len());

    // Dies if no samples are retained
    if samples_a.is_empty() || samples_b.is_empty() {
        eprintln!("No samples found in the lists.");
        std::process::exit(1);
    }

    // Define thread pool
    let pool = rayon::ThreadPoolBuilder::new().num_threads(*n_threads).build().unwrap();

    // Demo
    let p1: Vec<f32> = vec![0.001, 0.0002, 0.01, 0.4, 0.9];
    pool.install(|| {
        // Code here runs using at most num_threads threads
        // For example, a parallel iterator:
        let likelihood = compute_harding_likelihood(1, 100, 0.00528, 0.2, 0.7).expect("Cannot compute likelihood.");
        println!("Harding likelihood: {likelihood}");
    });

}
