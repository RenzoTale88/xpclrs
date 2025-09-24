/*
This module provides the functions required to compute the XP-CLR.
*/
use anyhow::Result;
use counter::Counter;
use itertools::Itertools;
use rand::prelude::IndexedRandom;
use rayon::prelude::*;
use rgsl::IntegrationWorkspace;
use rust_htslib::bcf::record::{Genotype, GenotypeAllele};
use scirs2_integrate::quad::{quad, QuadOptions};
use statistical::mean;
use statrs::distribution::{Binomial, Discrete};
use std::f64::consts::PI;

// Drafted bisect left and right
pub struct Bisector<'a, T> {
    data: &'a [T],
}

impl<'a, T: Ord> Bisector<'a, T> {
    pub fn new(data: &'a [T]) -> Self {
        Bisector { data }
    }

    pub fn bisect_left(&self, x: &T) -> usize {
        let mut low = 0;
        let mut high = self.data.len();
        while low < high {
            let mid = (low + high) / 2;
            if &self.data[mid] < x {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        low
    }

    pub fn bisect_right(&self, x: &T) -> usize {
        let mut low = 0;
        let mut high = self.data.len();
        while low < high {
            let mid = (low + high) / 2;
            if &self.data[mid] <= x {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        low
    }
}

// PartialOrd bisector
pub struct PartialBisector<'a, T> {
    data: &'a [T],
}

impl<'a, T: PartialOrd> PartialBisector<'a, T> {
    pub fn new(data: &'a [T]) -> Self {
        PartialBisector { data }
    }

    pub fn bisect_left(&self, x: &T) -> usize {
        let mut low = 0;
        let mut high = self.data.len();
        while low < high {
            let mid = (low + high) / 2;
            if self.data[mid].partial_cmp(x) == Some(std::cmp::Ordering::Less) {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        low
    }

    pub fn bisect_right(&self, x: &T) -> usize {
        let mut low = 0;
        let mut high = self.data.len();
        while low < high {
            let mid = (low + high) / 2;
            if self.data[mid].partial_cmp(x) != Some(std::cmp::Ordering::Greater) {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        low
    }
}

// Compute the omega
pub fn est_omega(q1: &[f64], q2: &[f64]) -> Result<f64> {
    if q2.contains(&0.0) || q2.contains(&1.0) {
        eprintln!("No SNPs in p2 can be fixed.");
        std::process::exit(1);
    };
    // Compute the omega
    let w = mean(
        &q1.iter()
            .zip(q2)
            .map(|(p, q)| ((p - q).powi(2)) / (q * (1f64 - q)))
            .collect::<Vec<f64>>(),
    );
    Ok(w)
}

// Measure variability for each SNP
pub fn var_estimate(w: f64, q2: f64) -> Result<f64> {
    log::debug!(
        "var_estimate {w} {q2} {} {} {}",
        1.0_f64 - q2,
        q2 * (1.0_f64 - q2),
        w * (q2 * (1.0_f64 - q2))
    );
    Ok(w * (q2 * (1.0_f64 - q2)))
}

// Rounding function
fn round_to(x: f64, digits: u32) -> f64 {
    let factor = 10_f64.powi(digits as i32);
    (x * factor).round() / factor
}

// c is a proxy for selection strength.
// c is [0,1] the closer to 1, the weaker the influence of selection.
pub fn compute_c(
    r: f64,
    s: f64,
    ne: Option<u64>,
    minrd: Option<f64>,
    sf: Option<u32>,
) -> Result<f64> {
    let ne = ne.unwrap_or(20000) as f64;
    let minrd = minrd.unwrap_or(1e-7);
    let sf = sf.unwrap_or(5);
    if s <= 0.0 {
        Ok(1.0)
    } else {
        let x = -((2.0 * ne).ln()) * (r.max(minrd)) / s;
        Ok(round_to(1.0 - (x.exp()), sf))
    }
}

/// Compute pdf(p1 | c, p2, var) matching the Python logic for scalar p1.
fn pdf_scalar(p1: f64, c: f64, p2: f64, var: f64) -> f64 {
    let p1 = vec![p1];

    // A-term
    let a_term: f64 = (2f64 * PI * var).sqrt().powf(-1.0);

    // Create the target vector
    let mut r: Vec<f64> = vec![0f64; p1.len()];

    // Extract values where p1 is greater then 1-c
    let bisector = PartialBisector::new(&p1);
    let left = bisector.bisect_left(&c);
    let right = bisector.bisect_right(&(1f64 - c));

    // left hand side
    let b_term_l = &p1[0..left]
        .iter()
        .map(|i| (c - i) / (c.powf(2f64)))
        .collect::<Vec<f64>>();

    let c_term_l = &p1[0..left]
        .iter()
        .map(|i| (i - (c * p2)).powf(2f64) / (2f64 * c.powf(2f64) * var))
        .collect::<Vec<f64>>();

    let l_slice = &mut r[..left];
    for ((l_i, &b), &c) in l_slice.iter_mut().zip(b_term_l).zip(c_term_l) {
        *l_i += a_term * b * (-c).exp();
    }

    // Repeat for right term
    let b_term_r = &p1[right..]
        .iter()
        .map(|i| (i + c - 1.0_f64) / (c.powf(2f64)))
        .collect::<Vec<f64>>();

    let c_term_r = &p1[right..]
        .iter()
        .map(|i| (i + c - 1.0_f64 - (c * p2)).powf(2f64) / (2f64 * c.powf(2f64) * var))
        .collect::<Vec<f64>>();

    let r_slice = &mut r[right..];
    for ((r_i, &b), &c) in r_slice.iter_mut().zip(b_term_r).zip(c_term_r) {
        *r_i += a_term * b * (-c).exp();
    }
    r[0]
}

/// pdf_integral(p1) = pdf(p1) * BinomPMF(xj | nj, p1)
fn pdf_integral_scalar(p1: f64, xj: u64, nj: u64, c: f64, p2: f64, var: f64) -> f64 {
    let dens = pdf_scalar(p1, c, p2, var);
    let binom = Binomial::new(p1, nj).expect("Cannot generate binomial distr.");
    // let pmf = binom.pmf(xj);
    let logpmf = binom.ln_pmf((xj as i64).try_into().unwrap());
    let pmf = logpmf.exp();
    dens * pmf
}

fn _integrate_qags_scirs2<F>(
    f: F,
    a: f64,
    b: f64,
    epsabs: f64,
    epsrel: f64,
    limit: usize,
) -> (f64, f64)
where
    F: Fn(f64) -> f64,
{
    // QAGS: adaptive Gauss–Kronrod with singularity handling
    let qopts = QuadOptions {
        abs_tol: epsabs,
        rel_tol: epsrel,
        max_evals: limit,
        use_abs_error: false,
        use_simpson: true
    };
    let result = quad(|x: f64| f(x), a, b, Some(qopts)).expect("Cannot converge");
    (result.value, result.abs_error)
}

fn _integrate_qags_gsl<F>(
    f: F,
    a: f64,
    b: f64,
    epsabs: f64,
    epsrel: f64,
    limit: usize,
) -> (f64, f64)
where
    F: Fn(f64) -> f64,
{
    // QAGS: adaptive Gauss–Kronrod with singularity handling
    let mut w = IntegrationWorkspace::new(10_000).expect("GSL workspace alloc failed");
    let func = |x: f64| f(x);
    let (res, err) = w
        .qags(func, a, b, epsabs, epsrel, limit)
        .expect("GSL qags failed");
    (res, err)
}

fn _integrate_qagp_gsl<F>(
    f: F,
    a: f64,
    b: f64,
    interior_points: &[f64], // known breakpoints (excluding endpoints)
    epsabs: f64,
    epsrel: f64,
    limit: usize,
) -> (f64, f64)
where
    F: Fn(f64) -> f64,
{
    // QAGP: adaptive Gauss–Kronrod with user-supplied breakpoints
    // pts must include endpoints a and b, sorted
    let mut pts = Vec::with_capacity(interior_points.len() + 2);
    pts.push(a);
    pts.extend_from_slice(interior_points);
    pts.push(b);
    pts.sort_by(|x, y| x.partial_cmp(y).unwrap());

    let mut w = IntegrationWorkspace::new(limit).expect("GSL workspace alloc failed");
    let func = |x: f64| f(x);
    let (res, err) = w
        .qagp(func, &mut pts, epsabs, epsrel, limit)
        .expect("GSL qagp failed");
    (res, err)
}

fn _compute_chen_likelihood(
    xj: u64,
    nj: u64,
    c: f64,
    p2: f64,
    var: f64,
    mode: Option<&str>,
) -> Result<f64> {
    // Integral bounds and tolerances matching SciPy quad
    let a = 0.001;
    let b = 0.999;
    // Use practical SciPy-like settings used in the Python implementation
    let epsabs = 0.0;      // SciPy default
    let epsrel = 0.001;     // Matches xpclr tolerance used elsewhere
    let limit = 50; // number of subintervals (QUADPACK-style, like SciPy)

    // Recommended: use QAGP with breakpoints at c and 1-c
    let breaks = [c, 1.0 - c];

    // Integral with binomial factor
    let (like_i, _err_i) = match mode {
        Some("qagp") =>  _integrate_qagp_gsl(
            |p1| pdf_integral_scalar(p1, xj, nj, c, p2, var),
            a,
            b,
            &breaks,
            epsabs,
            epsrel,
            limit,
        ),
        Some("qags") => _integrate_qags_gsl(
            |p1| pdf_integral_scalar(p1, xj, nj, c, p2, var),
            a,
            b,
            epsabs,
            epsrel,
            limit,
        ),
        _ =>  _integrate_qags_scirs2(
            |p1| pdf_integral_scalar(p1, xj, nj, c, p2, var),
            a,
            b,
            epsabs,
            epsrel,
            limit,
        ),
    };

    // Base integral of the pdf (denominator)
    let (like_b, _err_b) = match mode {
        Some("qagp") =>  _integrate_qagp_gsl(
            |p1| pdf_scalar(p1, c, p2, var),
            a,
            b,
            &breaks,
            epsabs,
            epsrel,
            limit,
        ),
        Some("qags") => _integrate_qags_gsl(
            |p1| pdf_scalar(p1, c, p2, var),
            a,
            b,
            epsabs,
            epsrel,
            limit,
        ),
        _ =>  _integrate_qags_scirs2(
            |p1| pdf_scalar(p1, c, p2, var),
            a,
            b,
            epsabs,
            epsrel,
            limit,
        ),
    };

    // Return the right value
    let ratio = if like_i > 0.0 && like_b > 0.0 {
        like_i.ln() - like_b.ln()
    } else {
        -1800.0
    };
    Ok(ratio)
}

// Compute composite likelihood
fn compute_complikelihood(
    sc_m: (f64, Option<usize>),
    xs: &[u64],
    ns: &[u64],
    rds: &[f64],
    p2freqs: &[f64],
    weights: &[f64],
    omegas: &[f64],
) -> Result<f64> {
    let sc = sc_m.0;
    let _method = sc_m.1.unwrap_or(0);
    if !(0.0..1.0).contains(&sc) {
        Ok(f64::INFINITY)
    } else {
        let marginall = itertools::izip!(xs, ns, rds, p2freqs, weights, omegas)
            .map(|(xj, nj, r, p2, weight, omega)| {
                // compute the variance
                let var = var_estimate(*omega, *p2).expect("Cannot compute variance");
                // Compute C
                let c = compute_c(*r, sc, None, None, Some(5_u32)).expect("Cannot compute C");
                // Compute likelihood
                // Use GSL QAGP with breakpoints to mirror SciPy's quad behavior
                let cl = _compute_chen_likelihood(*xj, *nj, c, *p2, var, Some("qags"))
                    .expect("Cannot compute the likelihood");
                // Return the weighted margin
                cl * *weight
            })
            .collect::<Vec<f64>>();
        // final value
        let ml: f64 = marginall.iter().sum();
        // println!("compute_complikelihood sc={sc} marginall={marginall:?} ml={ml}");
        Ok(-ml)
    }
}

fn compute_xpclr(
    counts: (&[u64], &[u64]),
    rds: &[f64],
    p2freqs: &[f64],
    weights: &[f64],
    omegas: &[f64],
    sel_coeffs: &[f64],
    method: Option<usize>,
) -> Result<(f64, f64, f64)> {
    // Extract counts
    let xs = counts.0;
    let ns = counts.1;

    // Define objectives
    let mut maximum_li = f64::INFINITY;
    let mut maxli_sc = 0.0f64;
    let mut null_model_li = f64::INFINITY;

    // Define selection coefficient
    for (counter, sc) in sel_coeffs.iter().enumerate() {
        // Compute ll
        let ll = compute_complikelihood((*sc, method), xs, ns, rds, p2freqs, weights, omegas)
            .expect("Cannot infer composite likelihood");
        if counter == 0 {
            null_model_li = ll;
        }
        log::debug!("compute_xpclr {counter} {sc} {ll}");
        // Replace values
        if ll < maximum_li {
            maximum_li = ll;
            maxli_sc = *sc;
        }
    }
    log::debug!("compute_xpclr {maximum_li} {null_model_li} {maxli_sc}\n\n");
    Ok((-maximum_li, -null_model_li, maxli_sc))
}

// Define indexes of variants in the window
fn get_window(
    pos: &[usize],
    start: usize,
    stop: usize,
    max_pos_size: usize,
) -> Result<(Vec<usize>, usize)> {
    // Define start index
    let start_ix = pos.binary_search(&start).unwrap_or_else(|i| i);
    let stop_ix = pos.binary_search(&stop).unwrap_or_else(|i| i);
    if (stop_ix - start_ix) > max_pos_size {
        // The window has too many sites; randomly select some
        let mut ix: Vec<usize> = (start_ix..stop_ix)
            .collect::<Vec<usize>>()
            .choose_multiple(&mut rand::rng(), max_pos_size)
            .cloned()
            .collect::<Vec<usize>>();
        ix.sort();
        Ok((ix, stop_ix - start_ix))
    } else {
        let ix = (start_ix..stop_ix).step_by(1).collect::<Vec<usize>>();
        Ok((ix, stop_ix - start_ix))
    }
}

// Compute A1/A2 counts and A2 frequency
fn pair_gt_to_af(
    gt1_m: &[Vec<Genotype>],
    gt2_m: &[Vec<Genotype>],
) -> Result<(
    Vec<u32>,
    Vec<u64>,
    Vec<u64>,
    Vec<f64>,
    Vec<u64>,
    Vec<u64>,
    Vec<f64>,
)> {
    let vals: Vec<(u32, u64, u64, f64, u64, u64, f64)> = gt1_m
        .iter()
        .zip(gt2_m)
        .map(|(gts1, gts2)| {
            let counts1 = gts1
                .iter()
                .flat_map(|gt| gt.iter().filter_map(|a| a.index()))
                .collect::<Counter<u32>>();
            let counts2 = gts2
                .iter()
                .flat_map(|gt| gt.iter().filter_map(|a| a.index()))
                .collect::<Counter<u32>>();
            let mut all_alleles: Counter<_> = counts1.clone();
            all_alleles.extend(&counts2);

            let ref_allele = all_alleles
                .keys()
                .min()
                .expect("Can't compute the alt allele count");
            let tot_counts1: usize = counts1.values().sum();
            let tot_counts2: usize = counts2.values().sum();
            let ref_counts1 = counts1.get(ref_allele).unwrap_or(&0);
            let ref_counts2 = counts2.get(ref_allele).unwrap_or(&0);
            let alt_counts1 = (tot_counts1 - ref_counts1) as f64;
            let alt_counts2 = (tot_counts2 - ref_counts2) as f64;
            (
                *ref_allele,
                tot_counts1 as u64,
                alt_counts1 as u64,
                alt_counts1 / (tot_counts1 as f64),
                tot_counts2 as u64,
                alt_counts2 as u64,
                alt_counts2 / (tot_counts2 as f64),
            )
        })
        .collect();
    Ok(Itertools::multiunzip(vals.into_iter()))
}

// Attempt to compute the LD using the same method as scikit-allele
// [here](https://github.com/cggh/scikit-allel/blob/master/allel/opt/stats.pyx#L90)
// and [here]()
fn gn_pairwise_corrcoef_int8(gn: &[Vec<i8>]) -> Result<Vec<Vec<f64>>> {
    let n = gn.len();
    // Precompute gn_sq[i][k] = gn[i][k]^2
    let gn_sq: Vec<Vec<i8>> = gn
        .iter()
        .map(|row| row.iter().map(|&v| v * v).collect())
        .collect();

    // Create square matrix
    let mut out = vec![vec![0.0_f64; n]; n];

    // Iterate through the variants
    for i in 0..(n - 1) {
        for j in (i + 1)..n {
            let gn0 = &gn[i];
            let gn1 = &gn[j];
            let gn0_sq = &gn_sq[i];
            let gn1_sq = &gn_sq[j];

            let r = gn_corrcoef_int8(gn0, gn1, gn0_sq, gn1_sq);
            out[i][j] = r.powi(2);
            out[j][i] = r.powi(2);
        }
    }

    Ok(out)
}

// Example reimplementation of gn_corrcoef_int8 in Rust
fn gn_corrcoef_int8(a: &[i8], b: &[i8], a_sq: &[i8], b_sq: &[i8]) -> f64 {
    // Convert to f64 and compute Pearson correlation
    let mut m0: f64 = 0.0;
    let mut m1: f64 = 0.0;
    let mut v0: f64 = 0.0;
    let mut v1: f64 = 0.0;
    let mut cov: f64 = 0.0;
    let mut n: f64 = 0.0;

    // perform sums
    for i in 0..a.len() {
        let x = a[i];
        let y = b[i];
        if x >= 0 && y >= 0 {
            n += 1.0f64;
            m0 += x as f64;
            m1 += y as f64;
            v0 += a_sq[i] as f64;
            v1 += b_sq[i] as f64;
            cov += (x * y) as f64;
        }
    }

    if n == 0.0 || v0 == 0.0 || v1 == 0.0 {
        return f64::NAN;
    }

    // Reproduce logic
    m0 /= n;
    m1 /= n;
    v0 /= n;
    v1 /= n;
    cov /= n;
    cov -= m0 * m1;
    v0 -= m0 * m0;
    v1 -= m1 * m1;

    // Return r
    cov / (v0 * v1).sqrt()
}

// Convert genotypes in the appropriate form
fn gt2haplotypes(gt_m: Vec<&Vec<Genotype>>) -> Vec<Vec<i8>> {
    gt_m.iter()
        .map(|gts| {
            gts.iter()
                .flat_map(|gt| {
                    gt.iter()
                        .map(|a| match a {
                            GenotypeAllele::PhasedMissing => -9_i8,
                            GenotypeAllele::UnphasedMissing => -9_i8,
                            _ => a.index().unwrap() as i8,
                        })
                        .collect::<Vec<i8>>()
                })
                .collect::<Vec<i8>>()
        })
        .collect::<Vec<Vec<i8>>>()
}

// Genotypes to counts
fn gt2gcounts(gt_m: Vec<&Vec<Genotype>>, ref_all: Vec<u32>) -> Vec<Vec<i8>> {
    gt_m.iter()
        .zip(ref_all.iter()) // iter over (Vec<Genotype>, &u32 ref allele index)
        .map(|(gts, &ref_ix)| {
            gts.iter()
                .map(|gt| {
                    // Extract allele indices, ignoring missing
                    let alleles: Vec<u32> = gt
                        .iter()
                        .filter_map(|a| a.index()) // skip missing
                        .collect();

                    if alleles.is_empty() {
                        // All missing
                        -9_i8
                    } else {
                        // Count how many are NOT the ref allele
                        let alt_count = alleles.iter().filter(|&&ix| ix != ref_ix).count() as i8;

                        if alt_count == 0 {
                            0
                        } else {
                            alt_count
                        }
                    }
                })
                .collect()
        })
        .collect()
}

// LD cutoff
fn apply_cutoff(matrix: &[Vec<f64>], cutoff: f64) -> Vec<Vec<bool>> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|&val| {
                    let cond1 = val > cutoff; // ld**2 > cutoff
                    let cond2 = val.is_nan(); // np.isnan(ld)
                    cond1 || cond2 // elementwise OR
                })
                .collect()
        })
        .collect()
}

// Compute the variant weight
fn compute_weights(
    gt_m: Vec<&Vec<Genotype>>,
    ref_all: Vec<u32>,
    ldcutoff: f64,
    isphased: bool,
) -> Result<Vec<f64>> {
    // First create a matrix that can be used to compute the LD
    let d = if isphased {
        gt2haplotypes(gt_m)
    } else {
        gt2gcounts(gt_m, ref_all)
    };

    // Compute the R2
    let ld = gn_pairwise_corrcoef_int8(&d).expect("Cannot compute LD");

    // Apply cutoff
    let above_cut = apply_cutoff(&ld, ldcutoff);

    // Return weight for each site
    let weights = above_cut
        .iter()
        .map(|v| {
            let summa: i32 = v
                .iter()
                .map(|b| match b {
                    true => 1,
                    false => 0,
                })
                .sum();
            1_f64 / ((summa + 1) as f64)
        })
        .collect::<Vec<f64>>();
    Ok(weights)
}

// Main XP-CLR caller
pub fn xpclr(
    (gt1, gt2): (Vec<Vec<Genotype>>, Vec<Vec<Genotype>>), // Genotypes
    (bpositions, geneticd, windows): (Vec<usize>, Vec<f64>, Vec<(usize, usize)>), // Positions
    (ldcutoff, phased): (Option<f64>, Option<bool>),      // LD-related
    (maxsnps, minsnps): (usize, usize),                   // Size/count filters
) -> Result<
    Vec<(
        usize,
        (usize, usize, usize, usize, usize, usize),
        (f64, f64, f64, f64),
    )>,
> {
    let sel_coeffs = vec![
        0.0, 0.00001, 0.00005, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.003, 0.005, 0.01,
        0.05, 0.08, 0.1, 0.15,
    ];

    let ldcutoff = ldcutoff.unwrap_or(0.95f64);
    let isphased = phased.unwrap_or(false);

    // Get the allele frequencies first
    let (ar, t1, a1, q1, _t2, _a2, q2) =
        pair_gt_to_af(&gt1, &gt2).expect("Failed to copmute the AF for pop 1");

    // Then, let's compute the omega
    let w = est_omega(&q1, &q2).expect("Cannot compute omega");
    log::info!("Omega: {w}");

    // Process each window
    let mut results: Vec<(usize, (usize, usize, usize, usize, usize, usize), (f64, f64, f64, f64))> = windows
        // Parallellize by window
        .par_iter()
        .enumerate()
        .map(|(n, (start, stop))| {
            let (ix, n_avail) = get_window(&bpositions, *start, *stop, maxsnps).expect("Cannot find the window");
            let n_snps = ix.len();
            let max_ix = ix.iter().last().unwrap_or(&0_usize).to_owned();
            log::debug!("xpclr Window idx: {n}; Window BP interval: {start}-{stop}; N SNPs selected: {n_snps}; N SNP available: {n_avail}");
            if n_snps < minsnps {
                let xpclr_vals = (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
                (n, (*start, *stop, *start, *stop, n_snps, n_avail), xpclr_vals)
            } else {
                let bpi = bpositions[ix[0]];
                let bpe = bpositions[max_ix];
                // Do not clone, just refer to them
                let (gt_range, ar_range, gd_range, a1_range, t1_range, p2freqs): (Vec<&Vec<Genotype>>, Vec<u32>, Vec<f64>, Vec<u64>, Vec<u64>, Vec<f64>) = Itertools::multiunzip(ix.iter().map(|&i| (&gt2[i], &ar[i], &geneticd[i], &a1[i], &t1[i], &q2[i])));
                // Compute distances from the average gen. dist.
                let mdist = mean(&gd_range);
                let rds = gd_range.iter().map(|d| (d - mdist).abs()).collect::<Vec<f64>>();

                // Compute the weights
                let weights = compute_weights(gt_range, ar_range, ldcutoff, isphased).expect("Failed to compute the weights");
                let omegas = vec![w; rds.len()];
                // Compute XP-CLR
                log::debug!("xpclr {start} {stop} {} {p2freqs:?}", rds.len());
                let xpclr_res = compute_xpclr(
                    (&a1_range, &t1_range),
                    &rds,
                    &p2freqs,
                    &weights,
                    &omegas,
                    &sel_coeffs,
                    None
                ).expect("Failed computing XP-CLR for window");
                let xpclr_v = 2.0_f64 * (xpclr_res.0 - xpclr_res.1);
                (n, (*start, *stop, bpi, bpe, n_snps, n_avail), (xpclr_res.0, xpclr_res.1, xpclr_res.2, xpclr_v))
            }
        })
        .collect();
    results.sort_by_key(|item| item.0);
    Ok(results)
}
