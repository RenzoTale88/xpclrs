/*
This module provides the functions required to compute the XP-CLR.
*/
use anyhow::Result;
use statistical::mean;
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

// given a set of potential p1s, we wish to know the probability density at each
// value of p1.
pub fn pdf(p1: &[f32], c: f32, p2: f32, var: f32) -> Result<Vec<f32>> {
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

/*
def pdf(p1, data):

    """
    param p1: x values for which we want pdens
    :param data is an array of c, p2, var
    :return:
    """

    # right hand side
    b_term_r = (p1[ixr:] + c - 1)/(c**2)
    c_term_r = ((p1[ixr:] + c - 1 - (c*p2))**2) / (2*(c**2)*var)
    r[ixr:] += a_term * b_term_r * np.exp(-c_term_r)

    return r
*/
