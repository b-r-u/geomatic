#[derive(Copy, Clone, Debug)]
pub struct EllipsoidParams {
    /// Semi-major axis
    pub semi_major_axis: f64,
    /// Inverse flattening (1 / f)
    pub f_inv: f64,
}

impl EllipsoidParams {
    #[inline(always)]
    pub fn flattening(&self) -> f64 {
        self.f_inv.recip()
    }

    #[inline(always)]
    pub fn eccentricity(&self) -> f64 {
        let f = self.flattening();
        (2.0 * f - f * f).sqrt()
    }
}

pub trait Ellipsoid {
    const PARAMS: EllipsoidParams;
}

#[derive(Copy, Clone, Debug)]
pub enum WGS84 {}

impl Ellipsoid for WGS84 {
    const PARAMS: EllipsoidParams = EllipsoidParams {
        semi_major_axis: 6378137.0,
        f_inv: 298.257223563,
    };
}

#[derive(Copy, Clone, Debug)]
pub enum GRS1980 {}

impl Ellipsoid for GRS1980 {
    const PARAMS: EllipsoidParams = EllipsoidParams {
        semi_major_axis: 6378137.0,
        f_inv: 298.2572221008827,
    };
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::{approx_eq, assert_approx_eq};


    #[test]
    fn eccentricity() {
        assert_approx_eq(WGS84::PARAMS.eccentricity(), 0.08181919084262, 1e-12);
        assert_approx_eq(GRS1980::PARAMS.eccentricity(), 0.0818191910435, 1e-12);
        assert!(!approx_eq(WGS84::PARAMS.eccentricity(), GRS1980::PARAMS.eccentricity(), 1e-12));
    }
}

