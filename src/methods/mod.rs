/*
This module provides the functions required to compute the XP-CLR.
*/
use anyhow::Result;
use quad_rs::{EvaluationError, Integrable, Integrator};
use rayon::prelude::*;
use statistical::mean;
use statrs::distribution::{Binomial, Discrete};
use std::f32::consts::PI;

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
pub fn omega_est(q1: &[f32], q2: &[f32]) -> Result<f32> {
    if q2.contains(&0.0) || q2.contains(&1.0) {
        eprintln!("No SNPs in p2 can be fixed.");
        std::process::exit(1);
    };
    // Compute the omega
    let w = mean(
        &q1.par_iter()
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
        .par_iter()
        .map(|i| (c - i) / (c.powf(2f32)))
        .collect::<Vec<f32>>();
    let c_term_l = &p1[0..left]
        .par_iter()
        .map(|i| (i - (c * p2)).powf(2f32) / (2f32 * c.powf(2f32) * var))
        .collect::<Vec<f32>>();
    let l_slice = &mut r[..left];
    for ((l_i, &b), &c) in l_slice.iter_mut().zip(b_term_l).zip(c_term_l) {
        *l_i += a_term * b * (-c).exp();
    }

    // Repeat for right term
    let b_term_r = &p1[right..]
        .par_iter()
        .map(|i| (c - i) / (c.powf(2f32)))
        .collect::<Vec<f32>>();
    let c_term_r = &p1[right..]
        .par_iter()
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
        .par_iter()
        .zip(dens)
        .map(|(p, d)| {
            d * Binomial::new(*p as f64, nj)
                .expect("Can't compute binomial distribution")
                .pmf(xj) as f32
        })
        .collect::<Vec<f32>>();
    Ok(binomials)
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

pub fn compute_chen_likelihood(
    n_alt_obs: u64,
    tot_alt_obs: u64,
    c: f32,
    p2: f32,
    var: f32,
) -> Result<f32> {
    // Default integrator
    let integrator = Integrator::default()
        .relative_tolerance(0.001);

    // Prepare integrands
    let integrand_pdf = Pdensity { c, p2, var };
    let integrand_bin = PdensityBinom {
        c,
        p2,
        var,
        xj: n_alt_obs,
        nj: tot_alt_obs,
    };

    // Integration range on p1 values, for example [0.0, 1.0]
    let like_b  = integrator
        .integrate(integrand_pdf, 0.001f32..0.999)
        .expect("Invalid integral")
        .result.result.expect("Not a valid number");
    let like_i = integrator
        .integrate(integrand_bin, 0.001f32..0.999)
        .expect("Invalid integral")
        .result.result.expect("Not a valid number");

    // Return the right value
    let ratio = match (like_i, like_b) {
        (0.0, _) => -1800f32,
        (_, 0.0) => -1800f32,
        (_, _) => like_i.ln() - like_b.ln()
    };
    Ok(ratio)
}

/*
# This is recoded from the implementation of xpclr. I do not think it represents
# what is discussed in the paper. Here we take the integral of eq 5 in the paper
# then take the ratio of it to eq 5 without the binomial component.
# additionally they neglect the probability of p1 being 0 or 1, I presume to
# allow the romberg integration to converge ok.
# This calculates the likelihood of a given SNP
@lru_cache(maxsize=2**16)
def chen_likelihood(values):

    """
    :param values: is an array of nobserved alt alleles, total obs alleles, c,
    p2freq, and var.
    :return: The likelihood ratio of the two likelihoods.
    """

    with warnings.catch_warnings(record=True) as w:

        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")

        # These 2 function calls are major speed bottlenecks.
        # to do: http://docs.scipy.org/doc/scipy/reference/tutorial/integrate
        # .html#faster-integration-using-ctypes
        i_likl = quad(pdf_integral, a=0.001, b=0.999, args=(values,),
                      epsrel=0.001, epsabs=0, full_output=1)

        i_base = quad(pdf, a=0.001, b=0.999, args=(values[2:],),
                      epsrel=0.001, epsabs=0, full_output=1)

        if w:
            logger.warning(w[-1].message, i_likl, i_base)

    like_i, like_b = i_likl[0], i_base[0]

    if like_i == 0.0:
        ratio = -1800
    elif like_b == 0.0:
        ratio = -1800
    else:
        ratio = np.log(like_i) - np.log(like_b)

    return ratio
*/
