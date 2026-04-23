use clap::{builder::PossibleValue, value_parser, Arg, ArgAction, Command};
use env_logger::{self, Env};
use rayon::current_num_threads;
use xpclrs::{
    io::{process_plink, process_xcf, read_file, to_table, write_table},
    methods::clr,
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
        .about("Compute the XP-CLR for a pair of populations from a VCF file.\nMethods presented by Chen H, Patterson N, Reich D. Population differentiation as a test for selective sweeps. Genome Res. 2010 Mar;20(3):393-402. doi: 10.1101/gr.100545.109. Epub 2010 Jan 19. PMID: 20086244; PMCID: PMC2840981.\nOriginal implementation is available at https://github.com/hardingnj/xpclr/\n")
        .arg(
            Arg::new("INPUT")
                .short('I')
                .long("input")
                .required(true)
                .help("input file(s)"),
        )
        .arg(
            Arg::new("OUT")
                .short('O')
                .long("out")
                .required(true)
                .help("Output file name."),
        )
        .arg(
            Arg::new("SAMPLES")
                .short('A')
                .long("samples")
                .required(true)
                .help("Samples in population. Path to file with each ID on a line."),
        )
        .arg(
            Arg::new("RECRATE")
                .long("rrate")
                .short('R')
                .required(false)
                .default_value("1e-8")
                .value_parser(value_parser!(f64))
                .help("Recombination rate per base."),
        )
        .arg(
            Arg::new("MAXSNPS")
                .long("maxsnps")
                .short('M')
                .required(false)
                .default_value("200")
                .value_parser(value_parser!(u64))
                .help("Max SNPs in a window."),
        )
        .arg(
            Arg::new("MINSNPS")
                .long("minsnps")
                .short('N')
                .required(false)
                .default_value("10")
                .value_parser(value_parser!(u64))
                .help("Min SNPs in a window."),
        )
        .arg(
            Arg::new("SIZE")
                .long("size")
                .required(false)
                .default_value("20000")
                .value_parser(value_parser!(u64))
                .help("Sliding window size."),
        )
        .arg(
            Arg::new("START")
                .long("start")
                .required(false)
                .default_value("1")
                .value_parser(value_parser!(u64))
                .help("Start position for the sliding windows."),
        )
        .arg(
            Arg::new("STOP")
                .long("stop")
                .required(false)
                .value_parser(value_parser!(u64))
                .help("Stop position for the sliding windows."),
        )
        .arg(
            Arg::new("STEP")
                .long("step")
                .required(false)
                .default_value("20000")
                .value_parser(value_parser!(u64))
                .help("Step size for the sliding windows."),
        )
        .arg(
            Arg::new("PHASED")
                .short('P')
                .long("phased")
                .required(false)
                .action(ArgAction::SetTrue)
                .help("Whether data is phased for more precise r2 calculation (does not work with --plink)."),
        )
        .arg(
            Arg::new("CHROM")
                .short('C')
                .long("chr")
                .required(true)
                .help("Chromosome to analyse."),
        )
        .arg(Arg::new("DISTKEYS").long("gdistkey").required(false).help(
            "Key in INFO field providing the genetic position of each variant in the VCF file",
        ))
        .arg(
            Arg::new("NTHREADS")
                .short('t')
                .long("threads")
                .required(false)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Number of threads to use"),
        )
        .arg(
            Arg::new("OUTFMT")
                .short('f')
                .long("format")
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
            Arg::new("PLINK")
                .long("plink")
                .required(false)
                .action(ArgAction::SetTrue)
                .help("Input is in PLINK binary format (.bed/.bim/.fam) rather than VCF/BCF; EXPERIMENTAL."),
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
    let chrom = matches
        .get_one::<String>("CHROM")
        .expect("Invalid chromosome code")
        .to_owned();
    let start = matches.get_one::<u64>("START").copied();
    let step = matches.get_one::<u64>("STEP").copied().unwrap();
    let size = matches.get_one::<u64>("SIZE").copied().unwrap();
    let end = matches.get_one::<u64>("STOP").copied();
    let minsnps = matches.get_one::<u64>("MINSNPS").copied().unwrap();
    let maxsnps = matches.get_one::<u64>("MAXSNPS").copied().unwrap();
    let rrate = matches.get_one::<f64>("RECRATE").copied();
    let distkey = matches.get_one::<String>("DISTKEYS").cloned();
    let phased = matches.get_one::<bool>("PHASED").copied();
    let plink = matches.get_one::<bool>("PLINK").copied().unwrap_or(false);

    // Get the input file
    let input_path = matches
        .get_one::<String>("INPUT")
        .expect("Input file is required")
        .to_owned();

    // Get the output path
    let out_path = matches
        .get_one::<String>("OUT")
        .expect("Output file is required")
        .to_owned();
    let out_fmt = matches
        .get_one::<String>("OUTFMT")
        .expect("Invalid output format")
        .to_owned();
    let n_threads = match matches
        .get_one::<usize>("NTHREADS")
        .expect("Number of threads invalid")
    {
        0 => current_num_threads(),
        v => *v,
    };

    // set up logging
    let log_level = matches
        .get_one::<String>("LOG")
        .expect("Log level not valid")
        .to_owned();
    env_logger::Builder::from_env(Env::default().default_filter_or(log_level)).init();

    // Initial logging
    log::info!("xpclrs v{version}");

    // Get the sample lists
    let sample = matches
        .get_one::<String>("SAMPLES")
        .expect("Sample file is required");

    // Load sample lists as an u8 array
    let samples = read_file(sample)
        .expect("Invalid file for samples A")
        .map(|s| s.expect("Invalid string"))
        .collect::<Vec<String>>();

    // Load VCF file
    let g_data = if plink {
        process_plink(
            input_path,
            &samples,
            None,
            &chrom,
            start,
            end,
            (phased, rrate),
        )
    } else {
        process_xcf(
            input_path,
            &samples,
            None,
            &chrom,
            start,
            end,
            (phased, rrate, distkey, n_threads),
        )
    }
    .expect("Failed loading genotype data");

    // Establish the windows
    let end = match end {
        Some(v) => v,
        None => *g_data.positions.iter().max().unwrap() as u64,
    };

    // Validate that the window size is appropriate for the genomic region.
    let start_pos = start.unwrap();
    if end < start_pos {
        eprintln!(
            "Invalid genomic region: end position ({}) is less than start position ({}).",
            end, start_pos
        );
        std::process::exit(1);
    }
    // Use the smaller value between the user-defined size and the region length
    // (avoid failing when the region is smaller than the window size).
    let size = size.min(end - start_pos);
    if size == 0 {
        eprintln!("Invalid window size: size must be greater than zero.");
        std::process::exit(1);
    }

    // Prepare the windows.
    let windows = (start_pos..end)
        .step_by(step as usize)
        .map(|v| (v as usize, (v + size - 1) as usize))
        .collect::<Vec<(usize, usize)>>();

    // Define thread pool
    log::info!("Using {n_threads} threads.");
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .unwrap();

    // Prepare output file
    // Run CLR
    pool.install(|| {
        let clr_res = clr(g_data, windows, maxsnps as usize, minsnps as usize, phased)
            .expect("Failed running the CLR function");
        // Write output
        log::info!("Writing output file to {out_path}...");
        let mut clr_tsv = write_table(&format!("{out_path}.{chrom}.clr"));
        let _ = to_table(&chrom, &clr_res, &mut clr_tsv, &out_fmt);
    });

    log::info!("CLR computation completed.")
}
