# xpclrs
A rust implementation of the XP-CLR method.
This implementation achieves near identical results in a fraction of the run time.
The software analyses chromosome 24 of the VarGoat dataset (777,865 total variants, 236,145 used for the analysis, with two groups of 32 and 22 individuals, respectively)  in 00m:26s and using 963Mb of memory, versus 55m:20s and of memory for the original implementation.


## Installation
Clone and compile the code:
```
git clone https://www.github.com/RenzoTale88/xpclrs
cd xpclrs
cargo build --release
```

## Input
The software requires the following mandatory options:
1. Input genotypes in VCF(.GZ)/BCF format with `-I-`/`--input`.
2. The lists of individuals in each group (one individual per line) with `-A`/`--samplesA` and `-B`/`--samplesB`.
3. The sequence to analyse with `-C`/`--chr`.

The VCF can optionally include a genetic distance key, that can be specified with the `--gdistkey [NAME]`. Alternatively, users can provide the recombination rate with the `-R`/`--rrate` option.

The analysis can be further sped up using multithreading with the `--threads`/`-t` option, followed by the number of threads to use. If set to 0, the software will try to use all the threads available.

An example command and the list of available options can be seen using `--help`:
```
$ xpclrs --help
Compute the XP-CLR for a pair of populations from a VCF file.

Usage: xpclrs [OPTIONS] --input <VCF> --out <OUT> --samplesA <SAMPLES_A> --samplesB <SAMPLES_B> --chr <CHROM>

Options:
  -I, --input <VCF>           input VCF file
  -O, --out <OUT>             Output file name.
  -A, --samplesA <SAMPLES_A>  Samples in population A. Path to file with each ID on a line.
  -B, --samplesB <SAMPLES_B>  Samples in population B. Path to file with each ID on a line.
  -R, --rrate <RECRATE>       Recombination rate per base. [default: 1e-8]
  -L, --ld <LDCUTOFF>         LD cutoff. [default: 0.95]
  -m, --maxsnps <MAXSNPS>     Max SNPs in a window. [default: 200]
  -N, --minsnps <MINSNPS>     Min SNPs in a window. [default: 10]
      --size <SIZE>           Sliding window size. [default: 20000]
      --start <START>         Start position for the sliding windows. [default: 1]
      --stop <STOP>           Stop position for the sliding windows.
      --step <STEP>           Step size for the sliding windows. [default: 20000]
  -C, --chr <CHROM>           Chromosome to analyse.
      --gdistkey <DISTKEYS>   Key in INFO field providing the genetic position of each variant in the VCF file
  -t, --threads <NTHREADS>    Number of threads to use [default: 1]
  -f, --format <OUTFMT>       Format to save the output (csv, tsv, txt) [default: tsv] [possible values: tsv, txt, csv]
  -l, --log <LOG>             Number of threads to use [default: info] [possible values: info, debug]
  -h, --help                  Print help
  -V, --version               Print version
```

### Demo data
Can test with the demo data in the original xpclr repository [here](https://github.com/hardingnj/xpclr/tree/master/fixture).