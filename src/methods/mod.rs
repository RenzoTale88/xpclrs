/*
This module provides the functions required to compute the XP-CLR.
*/
use anyhow::Result;
use counter::Counter;
use itertools::MultiUnzip;
use quad_rs::{EvaluationError, Integrable, Integrator};
use rand::prelude::IndexedRandom;
use rayon::prelude::*;
use rust_htslib::bcf::record::{Genotype, GenotypeAllele};
use scirs2_integrate::romberg::{romberg, RombergOptions};
use statistical::mean;
use statrs::distribution::{Binomial, Discrete};
use std::f32::{consts::PI};

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
pub fn est_omega(q1: &[f32], q2: &[f32]) -> Result<f32> {
    if q2.contains(&0.0) || q2.contains(&1.0) {
        eprintln!("No SNPs in p2 can be fixed.");
        std::process::exit(1);
    };
    // Compute the omega
    let w = mean(
        &q1.iter()
            .zip(q2)
            .map(|(p, q)| ((p - q) * (p - q)) / (q * (1f32 - q)))
            .collect::<Vec<f32>>(),
    );
    Ok(w)
}

// Measure variability for each SNP
pub fn var_estimate(w: f32, q2: f32) -> Result<f32> {
    Ok(w * (q2 * (1f32 - q2)))
}

// c is a proxy for selection strength.
// c is [0,1] the closer to 1, the weaker the influence of selection.
pub fn compute_c(
    r: f32,
    s: f32,
    ne: Option<u64>,
    minrd: Option<f32>,
    sf: Option<u64>,
) -> Result<f32> {
    let ne = ne.unwrap_or(20000) as f32;
    let minrd = minrd.unwrap_or(1e-7);
    let _sf = sf.unwrap_or(5) as f32;
    if s <= 0.0 {
        Ok(1.0)
    } else {
        let x = -((2.0 * ne).ln()) * (r.max(minrd)) / s;
        Ok(1.0 - (x.exp()))
    }
}

// Standard function
fn compute_pdens(p1: &[f32], c: f32, p2: f32, var: f32) -> Result<Vec<f32>> {
    // First term
    let a_term: f32 = (2f32 * PI * var).sqrt().powf(-1.0);

    // Create the target vector
    let mut r: Vec<f32> = vec![0f32; p1.len()];

    // Extract values where p1 is greater then 1-c
    let bisector = PartialBisector::new(p1);
    let left = bisector.bisect_left(&c);
    let right = bisector.bisect_right(&(1f32 - c));

    // left hand side
    let b_term_l = &p1[0..left]
        .iter()
        .map(|i| (c - i) / (c.powf(2f32)))
        .collect::<Vec<f32>>();
    let c_term_l = &p1[0..left]
        .iter()
        .map(|i| (i - (c * p2)).powf(2f32) / (2f32 * c.powf(2f32) * var))
        .collect::<Vec<f32>>();
    let l_slice = &mut r[..left];
    for ((l_i, &b), &c) in l_slice.iter_mut().zip(b_term_l).zip(c_term_l) {
        *l_i += a_term * b * (-c).exp();
    }

    // Repeat for right term
    let b_term_r = &p1[right..]
        .iter()
        .map(|i| (c - i) / (c.powf(2f32)))
        .collect::<Vec<f32>>();
    let c_term_r = &p1[right..]
        .iter()
        .map(|i| (i - (c * p2)).powf(2f32) / (2f32 * c.powf(2f32) * var))
        .collect::<Vec<f32>>();
    let r_slice = &mut r[right..];
    for ((r_i, &b), &c) in r_slice.iter_mut().zip(b_term_r).zip(c_term_r) {
        *r_i += a_term * b * (-c).exp();
    }

    Ok(r)
}

pub fn pdens_binomial(p1: &[f32], xj: u64, nj: u64, c: f32, p2: f32, var: f32) -> Result<Vec<f32>> {
    // Compute dens first
    let dens = compute_pdens(p1, c, p2, var).expect("Can't compute dens");

    // Apply binomial function
    let binomials = p1
        .iter()
        .zip(dens)
        .map(|(p, d)| {
            d * Binomial::new(*p as f64, nj)
                .expect("Can't compute binomial distribution")
                .pmf(xj) as f32
        })
        .collect::<Vec<f32>>();
    Ok(binomials)
}

pub fn pden_binomial(p1_val: f32, xj: u64, nj: u64, c: f32, p2: f32, var: f32) -> Result<f32> {
    // a_term = 1 / sqrt(2 * PI * var)
    let a_term = 1.0 / (2.0 * PI * var).sqrt();

    // decide which branch depending on p1_val in relation to c and 1-c
    // You may adapt the bisect logic or simply partition manually here for scalar input

    // For example, approximate the b_term and c_term as per formulas:
    let (b_term, c_term) = if p1_val < c {
        let b = (c - p1_val) / (c * c);
        let c_t = (p1_val - c * p2).powi(2) / (2.0 * c.powi(2) * var);
        (b, c_t)
    } else if p1_val > 1.0 - c {
        let b = (p1_val + c - 1.0) / (c * c);
        let c_t = (p1_val + c - 1.0 - c * p2).powi(2) / (2.0 * c.powi(2) * var);
        (b, c_t)
    } else {
        // Between c and 1 - c, density is zero (or negligible)
        return Ok(0.0);
    };

    // Compute Gaussian-like term
    let density = a_term * b_term * (-c_term).exp();

    // Apply binomial pmf weight
    // Use statrs crate Binomial for pmf
    use statrs::distribution::Binomial;
    let binom = Binomial::new(p1_val as f64, nj).expect("Cannot compute binomial");
    let pmf_val = binom.pmf(xj) as f32;

    let result = density * pmf_val;

    Ok(result)
}

struct Pdensity {
    c: f32,
    p2: f32,
    var: f32,
}

impl Integrable for Pdensity {
    type Input = f32;
    type Output = f32;

    fn integrand(&self, input: &Self::Input) -> Result<Self::Output, EvaluationError<Self::Input>> {
        let c = self.c;
        let p2 = self.p2;
        let var = self.var;

        // Compute p_density component at this input point "input"
        // We rewrite compute_pdens logic for a single scalar input
        let p1_val = *input;

        // a_term = 1 / sqrt(2 * PI * var)
        let a_term = 1.0 / (2.0 * PI * var).sqrt();

        // decide which branch depending on p1_val in relation to c and 1-c
        // You may adapt the bisect logic or simply partition manually here for scalar input

        // For example, approximate the b_term and c_term as per formulas:
        let (b_term, c_term) = if p1_val < c {
            let b = (c - p1_val) / (c * c);
            let c_t = (p1_val - c * p2).powi(2) / (2.0 * c.powi(2) * var);
            (b, c_t)
        } else if p1_val > 1.0 - c {
            let b = (p1_val + c - 1.0) / (c * c);
            let c_t = (p1_val + c - 1.0 - c * p2).powi(2) / (2.0 * c.powi(2) * var);
            (b, c_t)
        } else {
            // Between c and 1 - c, density is zero (or negligible)
            return Ok(0.0);
        };

        // Compute Gaussian-like term
        let density = a_term * b_term * (-c_term).exp();

        Ok(density)
    }
}

struct PdensityBinom {
    c: f32,
    p2: f32,
    var: f32,
    xj: u64,
    nj: u64,
}

impl Integrable for PdensityBinom {
    type Input = f32;
    type Output = f32;

    fn integrand(&self, input: &Self::Input) -> Result<Self::Output, EvaluationError<Self::Input>> {
        let c = self.c;
        let p2 = self.p2;
        let var = self.var;
        let xj = self.xj;
        let nj = self.nj;

        // Compute p_density component at this input point "input"
        // We rewrite compute_pdens logic for a single scalar input
        let p1_val = *input;

        // a_term = 1 / sqrt(2 * PI * var)
        let a_term = 1.0 / (2.0 * PI * var).sqrt();

        // decide which branch depending on p1_val in relation to c and 1-c
        // You may adapt the bisect logic or simply partition manually here for scalar input

        // For example, approximate the b_term and c_term as per formulas:
        let (b_term, c_term) = if p1_val < c {
            let b = (c - p1_val) / (c * c);
            let c_t = (p1_val - c * p2).powi(2) / (2.0 * c.powi(2) * var);
            (b, c_t)
        } else if p1_val > 1.0 - c {
            let b = (p1_val + c - 1.0) / (c * c);
            let c_t = (p1_val + c - 1.0 - c * p2).powi(2) / (2.0 * c.powi(2) * var);
            (b, c_t)
        } else {
            // Between c and 1 - c, density is zero (or negligible)
            return Ok(0.0);
        };

        // Compute Gaussian-like term
        let density = a_term * b_term * (-c_term).exp();

        // Apply binomial pmf weight
        // Use statrs crate Binomial for pmf
        use statrs::distribution::Binomial;
        let binom = Binomial::new(p1_val as f64, nj).expect("Cannot compute binomial");
        let pmf_val = binom.pmf(xj) as f32;

        let result = density * pmf_val;

        Ok(result)
    }
}

fn compute_chen_likelihood(xj: u64, nj: u64, c: f32, p2: f32, var: f32) -> Result<f32> {
    // Default integrator
    let integrator = Integrator::default().relative_tolerance(0.001);

    // Prepare integrands
    let integrand_pdf = Pdensity { c, p2, var };
    let integrand_bin = PdensityBinom { c, p2, var, xj, nj };

    // Integration range on p1 values, for example [0.0, 1.0]
    let like_b = integrator
        .integrate(integrand_pdf, 0.001f32..0.999)
        .expect("Invalid integral")
        .result
        .result
        .expect("Not a valid number");
    let like_i = integrator
        .integrate(integrand_bin, 0.001f32..0.999)
        .expect("Invalid integral")
        .result
        .result
        .expect("Not a valid number");

    // Return the right value
    let ratio = match (like_i, like_b) {
        (0.0, _) => -1800f32,
        (_, 0.0) => -1800f32,
        (_, _) => like_i.ln() - like_b.ln(),
    };
    Ok(ratio)
}

// Second likelihood method
fn compute_romberg_likelihood(xj: u64, nj: u64, c: f32, p2: f32, var: f32) -> Result<f32> {
    let options = RombergOptions {
        max_iters: 50,                 // corresponds to divmax in SciPy
        abs_tol: 1.48e-6,              // absolute tolerance
        rel_tol: 1.48e-6,              // relative tolerance
        max_true_dimension: 3,         // default
        min_monte_carlo_samples: 1000, // fallback param
    };

    // Wrapper to run the function correctly
    let pdf_integral = |p1_val: f32| {
        // Use the captured `values` here to compute the integrand at x
        // For example, you can call your Rust implementation of pdf_integral:
        pden_binomial(p1_val, xj, nj, c, p2, var).expect("Cannot compute pdensity")
    };

    // Compute romberg integration
    let cl = romberg(pdf_integral, 0.0, 1.0, Some(options))?.value.ln();
    // Return acceptable value
    if cl.is_finite() {
        Ok(cl)
    } else {
        Ok(-1800f32)
    }
}

// Compute composite likelihood
fn compute_complikelihood(
    sc_m: (f32, Option<usize>),
    xs: &[u64],
    ns: &[u64],
    rds: &[f32],
    p2freqs: &[f32],
    weights: &[f32],
    omegas: &[f32],
) -> Result<f32> {
    let sc = sc_m.0;
    let method = sc_m.1.unwrap_or(0);
    if !(0.0..1.0).contains(&sc) {
        Ok(f32::INFINITY)
    } else {
        let marginall = &(xs, ns, rds, p2freqs, weights, omegas)
            .into_par_iter()
            .map(|(xj, nj, r, p2, weight, omega)| {
                // compute the variance
                let var = var_estimate(*omega, *p2).expect("Cannot compute variance");
                // Compute C
                let c = compute_c(*r, sc, None, None, None).expect("Cannot compute C");
                // Compute likelihood
                let cl = match method {
                    0 => compute_chen_likelihood(*xj, *nj, c, *p2, var),
                    1 => compute_romberg_likelihood(*xj, *nj, c, *p2, var),
                    _ => compute_chen_likelihood(*xj, *nj, c, *p2, var),
                }
                .expect("Cannot compute the likelihood");
                // Return the weighted margin
                cl * *weight
            })
            .collect::<Vec<f32>>();
        // final value
        let ml: f32 = marginall.iter().sum();
        Ok(-ml)
    }
}

fn compute_xpclr(
    counts: (&[u64], &[u64]),
    rds: &[f32],
    p2freqs: &[f32],
    weights: &[f32],
    omegas: &[f32],
    sel_coeffs: &[f32],
    method: Option<usize>,
) -> Result<(f32, f32, f32)> {
    // Extract counts
    let xs = counts.0;
    let ns = counts.1;

    // Define objectives
    let mut maximum_li = f32::INFINITY;
    let mut maxli_sc = 0.0f32;
    let mut null_model_li = f32::INFINITY;

    // Define selection coefficient
    for (counter, sc) in sel_coeffs.iter().enumerate() {
        // Compute ll
        let ll = compute_complikelihood((*sc, method), xs, ns, rds, p2freqs, weights, omegas)
            .expect("Cannot infer composite likelihood");
        if counter == 0 {
            null_model_li = ll;
        }
        // Replace values
        if ll < maximum_li {
            maximum_li = ll;
            maxli_sc = *sc;
        } else {
            break;
        }
    }
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
fn gt_to_af(gt_m: &[Vec<Genotype>]) -> Result<(Vec<u32>, Vec<u64>, Vec<u64>, Vec<f32>)> {
    let vals: Vec<(u32, u64, u64, f32)> = gt_m
        .iter()
        .map(|gts| {
            let counts = gts
                .iter()
                .flat_map(|gt| gt.iter().filter_map(|a| a.index()))
                .collect::<Counter<u32>>();
            let tot_counts: usize = counts.values().sum();
            let ref_allele = counts
                .keys()
                .min()
                .expect("Can't compute the minor allele count");
            let ref_counts = counts
                .get(ref_allele)
                .expect("Cannot compute ref allele frequency");
            let alt_counts = (tot_counts - ref_counts) as f32;
            (
                *ref_allele,
                tot_counts as u64,
                alt_counts as u64,
                alt_counts / (tot_counts as f32),
            )
        })
        .collect();
    Ok(vals.into_iter().multiunzip())
}

// Attempt to compute the LD using the same method as scikit-allele 
// [here](https://github.com/cggh/scikit-allel/blob/master/allel/opt/stats.pyx#L90)
// and [here]()
fn gn_pairwise_corrcoef_int8(gn: &[Vec<i8>]) -> Result<Vec<Vec<f32>>> {
    let n = gn.len();
    // Precompute gn_sq[i][k] = gn[i][k]^2
    let gn_sq: Vec<Vec<i8>> = gn
        .iter()
        .map(|row| row.iter().map(|&v| v * v).collect())
        .collect();

    // Create square matrix
    let mut out = vec![vec![1.0_f32; n]; n];

    // Iterate through the variants
    for i in 0..(n-1) {
        for j in (i + 1)..n {
            let gn0 = &gn[i];
            let gn1 = &gn[j];
            let gn0_sq = &gn_sq[i];
            let gn1_sq = &gn_sq[j];

            let r = gn_corrcoef_int8(gn0, gn1, gn0_sq, gn1_sq);
            out[i][j] = r;
            out[j][i] = r;
        }
    };

    Ok(out)
}

// Example reimplementation of gn_corrcoef_int8 in Rust
fn gn_corrcoef_int8(a: &[i8], b: &[i8], a_sq: &[i8], b_sq: &[i8]) -> f32 {
    // Convert to f32 and compute Pearson correlation
    let mut m0: f32 = 0.0;
    let mut m1: f32 = 0.0;
    let mut v0: f32 = 0.0;
    let mut v1: f32 = 0.0;
    let mut cov: f32 = 0.0;
    let mut n: f32 = 0.0;

    // perform sums
    for i in 0..a.len() {
        let x = a[i];
        let y = b[i];
        if x >= 0 && y >= 0 {
            n += 1.0f32;
            m0 += x as f32;
            m1 += y as f32;
            v0 += a_sq[i] as f32;
            v1 += b_sq[i] as f32;
            cov += (x * y) as f32;
        }
    };

    if n == 0.0 || v0 == 0.0 || v1 == 0.0 {
        return f32::NAN;
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
fn apply_cutoff(matrix: &[Vec<f32>], cutoff: f32) -> Vec<Vec<bool>> {
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
    ldcutoff: f32,
    isphased: bool,
) -> Result<Vec<f32>> {
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
            let summa: f32 = v
                .iter()
                .map(|b| match b {
                    true => 1_f32,
                    false => 0_f32,
                })
                .sum();
            1_f32 / (summa + 1.0)
        })
        .collect::<Vec<f32>>();
    Ok(weights)
}

// Main XP-CLR caller
pub fn xpclr(
    (gt1, gt2): (Vec<Vec<Genotype>>, Vec<Vec<Genotype>>), // Genotypes
    (bpositions, geneticd, windows): (Vec<usize>, Vec<f32>, Vec<(usize, usize)>), // Positions
    (ldcutoff, phased): (Option<f32>, Option<bool>), // LD-related
    (maxsnps, minsnps): (usize, usize),                   // Size/count filters
) -> Result<Vec<(usize, (usize, usize, usize, usize, usize, usize), (f32, f32, f32, f32))>> {
    let sel_coeffs = vec![
        0.0, 0.00001, 0.00005, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.003, 0.005, 0.01,
        0.05, 0.08, 0.1, 0.15,
    ];
    let ldcutoff = ldcutoff.unwrap_or(0.95f32);
    let isphased = phased.unwrap_or(false);

    // Get the allele frequencies first
    let (ar, t1, a1, q1) = gt_to_af(&gt1).expect("Failed to copmute the AF for pop 1");
    let (_, _t2, _a2, q2) = gt_to_af(&gt2).expect("Failed to copmute the AF for pop 1");

    // Compute genetic distances

    // Then, let's compute the omega
    let w = est_omega(&q1, &q2).expect("Cannot compute omega");
    log::info!("Omega: {w}");

    // Process each window
    let mut results: Vec<(usize, (usize, usize, usize, usize, usize, usize), (f32, f32, f32, f32))> = windows
        // Parallellize by window
        .par_iter()
        .enumerate()
        .map(|(n, (start, stop))| {
            let (ix, n_avail) = get_window(&bpositions, *start, *stop, maxsnps).expect("Cannot find the window");
            let n_snps = ix.len();
            let max_ix = ix.iter().last().unwrap_or(&0_usize).to_owned();
            log::debug!("Window idx: {n}; Window BP interval: {start}-{stop}; N SNPs selected: {n_snps}; N SNP available: {n_avail}");
            if n_snps < minsnps {
                let xpclr_vals = (f32::NAN, f32::NAN, f32::NAN, f32::NAN);
                (n, (*start, *stop, *start, *stop, n_snps, n_avail), xpclr_vals)
            } else {
                let bpi = bpositions[ix[0]];
                let bpe = bpositions[max_ix];
                // Do not clone, just refer to them
                let (gt_range, ar_range, gd_range, a1_range, t1_range, p2freqs): (Vec<&Vec<Genotype>>, Vec<u32>, Vec<f32>, Vec<u64>, Vec<u64>, Vec<f32>) = ix.iter().map(|&i| (&gt2[i], &ar[i], &geneticd[i], &t1[i], &a1[i], &q2[i])).multiunzip();
                // Compute distances from the average gen. dist.
                let mdist = mean(&gd_range);
                let rds = gd_range.iter().map(|d| d - mdist ).collect::<Vec<f32>>();
                // Compute the weights
                let weights = compute_weights(gt_range, ar_range, ldcutoff, isphased).expect("Failed to compute the weights");
                let omegas = vec![w; rds.len()];
                // Compute XP-CLR
                let xpclr_res = compute_xpclr(
                    (&a1_range, &t1_range),
                    &rds,
                    &p2freqs,
                    &weights,
                    &omegas,
                    &sel_coeffs,
                    Some(0)
                ).expect("Failed computing XP-CLR for window");
                let xpclr_v = 2.0_f32 * (xpclr_res.0 - xpclr_res.1);
                (n, (*start, *stop, bpi, bpe, n_snps, n_avail), (xpclr_res.0, xpclr_res.1, xpclr_res.2, xpclr_v))
            }
        })
        .collect();

    results.sort_by_key(|item| item.0);
    Ok(results)
}

