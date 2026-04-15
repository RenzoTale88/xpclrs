# xpclrs
A rust implementation of the XP-CLR method.
This implementation achieves near identical results in a fraction of the run time.
The software analyses chromosome 24 of the VarGoat dataset (777,865 total variants, 236,145 used for the analysis, with two groups of 32 and 22 individuals, respectively) in 00m:26s and using 77Mb of memory, versus 55m:20s and 321Mb of the original implementation.

## Installation
The compilation of the software requires the following packages to be installed:
1. [openblas](https://github.com/OpenMathLib/OpenBLAS)
2. [libclang](https://clang.llvm.org/)
3. [curl](https://curl.se/)
4. [rust](https://rust-lang.org/)

Then, install with cargo:
```
cargo install xpclrs
```

Check that the package is successfully installed with:
```
xpclrs --help
```

### Install with Docker/Singularity
The software is also available as a docker container in dockerhub. You can install it by pulling the image with docker:
```
docker pull tale88/xpclrs:latest
```
Or with singularity:
```
singularity build xpclrs.sif tale88/xpclrs:latest
```

## Input
The software requires the following mandatory options:
1. Input genotypes in VCF(.GZ)/BCF format with `-I-`/`--input`.
    * PLINK binary files (BED/BIM/FAM) are also supported by providing the root of the file name with the same `-I`/`--input` option and adding the `--plink` flag..
    * Loading in plink file is substantially faster than using the VCF format, but worth noticing that it can lead to different results due to the variants being coded as major/minor rather than REF/ALT (XP-CLR relies on allele frequencies).
2. The lists of individuals in each group (one individual per line) with `-A`/`--samplesA` and `-B`/`--samplesB`.
    * PLINK samples are loaded as `FID_IID`. So if your sample in the FAM file is `POP1 SAMP1 0 0 0 -9`, the sample will be listed as `POP1_SAMP1` in the group of individuals.
3. The sequence to analyse with `-C`/`--chr`.

The VCF can optionally include a genetic distance key (provided with `--gdistkey [NAME]`). Alternatively, users can provide the recombination rate with the `-R`/`--rrate` option. 
For PLINK inputs, the software will automatically detect the presence of a genetic position in the dataset and use that; if the value is equal to 0, the script will compute the genetic position based on the physical position and the recombination rate. Ensure that there are not gaps in the genetic position (i.e. a 0 following a known genetic position).

## Running the tool
The list of available options for `xpclrs` can be seen using `--help`:
```
$ xpclrs --help
Compute the XP-CLR for a pair of populations from a VCF file.
Methods presented by Chen H, Patterson N, Reich D. Population differentiation as a test for selective sweeps. Genome Res. 2010 Mar;20(3):393-402. doi: 10.1101/gr.100545.109. Epub 2010 Jan 19. PMID: 20086244; PMCID: PMC2840981.
Original implementation is available at https://github.com/hardingnj/xpclr/


Usage: xpclrs [OPTIONS] --input <INPUT> --out <OUT> --samplesA <SAMPLES_A> --samplesB <SAMPLES_B> --chr <CHROM>

Options:
  -I, --input <INPUT>         input file(s)
  -O, --out <OUT>             Output file name.
  -A, --samplesA <SAMPLES_A>  Samples in population A. Path to file with each ID on a line.
  -B, --samplesB <SAMPLES_B>  Samples in population B. Path to file with each ID on a line.
  -R, --rrate <RECRATE>       Recombination rate per base. [default: 1e-8]
  -L, --ld <LDCUTOFF>         LD cutoff. [default: 0.95]
  -M, --maxsnps <MAXSNPS>     Max SNPs in a window. [default: 200]
  -N, --minsnps <MINSNPS>     Min SNPs in a window. [default: 10]
      --size <SIZE>           Sliding window size. [default: 20000]
      --start <START>         Start position for the sliding windows. [default: 1]
      --stop <STOP>           Stop position for the sliding windows.
      --step <STEP>           Step size for the sliding windows. [default: 20000]
  -P, --phased                Whether data is phased for more precise r2 calculation (does not work with --plink).
  -C, --chr <CHROM>           Chromosome to analyse.
      --gdistkey <DISTKEYS>   Key in INFO field providing the genetic position of each variant in the VCF file
  -t, --threads <NTHREADS>    Number of threads to use [default: 1]
  -f, --format <OUTFMT>       Format to save the output (csv, tsv, txt) [default: tsv] [possible values: tsv, txt, csv]
  -F, --fast                  Run analysis in fast mode (faster integration, but gives results that are less accurate compared with the original tool)
      --plink                 Input is in PLINK binary format (.bed/.bim/.fam) rather than VCF/BCF; EXPERIMENTAL.
  -l, --log <LOG>             Logging level. [default: info] [possible values: info, debug]
  -h, --help                  Print help
  -V, --version               Print version```

Users can perform a trial run on the demo data provided by this repository with the command:
```
xpclrs --input test/test.vcf.gz --out test --samplesA test/samplesA.txt --samplesB test/samplesB.txt --chr chr1
```

It is possible to run the same analysis with multiple cores by setting `--threads`/`-t` to a higher integer value (if set to 0, the software will try to use all the threads available):
```
xpclrs --input test/test.vcf.gz --out test --samplesA test/samplesA.txt --samplesB test/samplesB.txt --chr chr1 --threads 4
```

The software can consider the phase when computing the linkage disequilibrium by providing the `--phased` option.

When providing inputs in plink binary format, users need to provide the appropriate `--plink` option:
```
xpclrs --input test/plink --plink --out test --samplesA test/samplesA_plink.txt --samplesB test/samplesB_plink.txt --chr 1
```

Finally, the software provides a `--fast` option, that disable the adaptive integration and provides approximate results. This speeds up the software significantly, but results may vary when compared with the original implementation.


### Demo data
The tool comes with a demo data generated from the [1000GP dataset](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). 
It is also possible to test the tool using the demo data in the [original xpclr repository](https://github.com/hardingnj/xpclr/tree/master/fixture).

## Citation
If you use the tool, please cite:
> Chen H, Patterson N, Reich D. Population differentiation as a test for selective sweeps. Genome Res. 2010 Mar;20(3):393-402. doi: 10.1101/gr.100545.109. Epub 2010 Jan 19. PMID: 20086244; PMCID: PMC2840981.
> The original xpclr tool [here](https://github.com/hardingnj/xpclr)
> Talenti A. XPCLRS: Fast Selection Signature Detection Using Cross-Population Composite Likelihood Ratio. bioRxiv 2026.02.27.708459. doi: 10.64898/2026.02.27.708459.
