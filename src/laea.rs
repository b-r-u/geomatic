use crate:: coord::{Point, Point3035, Point4326, CrsLaea};
use crate::ellipsoid::Ellipsoid;


pub struct LaeaParams {
    /// Center of projection latitude
    pub center_lat: f64,
    /// Center of projection longitude
    pub center_lon: f64,
    /// False origin easting
    pub false_easting: f64,
    /// False origin northing
    pub false_northing: f64,
}

#[inline(always)]
pub fn forward_gen<C: CrsLaea>(input: Point4326) -> Point<C> {
    let phi = input.lat().to_radians();
    let lam = input.lon().to_radians();

    let params = C::PARAMS;
    let ellips = C::Ellipsoid::PARAMS;

    let phi0 = params.center_lat.to_radians();
    let lam0 = params.center_lon.to_radians();
    let feast = params.false_easting;
    let fnorth = params.false_northing;
    let a = ellips.semi_major_axis;
    let e = ellips.eccentricity();

    let qp = (1.0 - e * e) * ((1.0 / (1.0 - e * e)) - ((1.0 / (2.0 * e)) * ((1.0 - e) / (1.0 + e)).ln()));
    let q0 = (1.0 - e * e) * (phi0.sin() / (1.0 - e*e*phi0.sin().powi(2)) - (1.0 / (2.0 * e) * ((1.0 - e * phi0.sin()) / (1.0 + e * phi0.sin())).ln()));
    let q = (1.0 - e * e) * (phi.sin() / (1.0 - e*e*phi.sin().powi(2)) - (1.0 / (2.0 * e) * ((1.0 - e * phi.sin()) / (1.0 + e * phi.sin())).ln()));

    let beta = (q / qp).asin();
    let beta0 = (q0 / qp).asin();
    let rq = a * (0.5 * qp).sqrt();

    let d = a * (phi0.cos() / (1.0 - e*e * phi0.sin().powi(2)).sqrt()) / (rq * beta0.cos());
    let b = rq * (2.0 / (1.0 + beta0.sin() * beta.sin() + (beta0.cos() * beta.cos() * (lam - lam0).cos()))).sqrt();

    let east = feast + (b * d) * (beta.cos() * (lam - lam0).sin());
    let north = fnorth + (b / d) * (beta0.cos() * beta.sin() - (beta0.sin() * beta.cos() * (lam - lam0).cos()));
    Point::<C>::new(east, north)
}


pub fn forward_ref(input: Point4326) -> Point3035 {
    forward_gen(input)
}

pub fn forward(input: Point4326) -> Point3035 {
    // constants for EPSG:3035
    let lam0 = 0.174532925_f64; // == 10°
    let feast = 4321000.0_f64;
    let fnorth = 3210000.0_f64;
    // GRS 1980 ellipsoid
    let e = 0.081819191_f64;

    // const computed
    let qp = 1.995531087485621_f64;
    let rq = 6371007.180890992_f64;
    let beta0 = 0.9053975167810719_f64;
    let d = 1.0004253945276074_f64;


    // inputs: e, qp, rq, beta0, lam0, d, feast, fnorth


    // input dependent
    let phi = input.lat().to_radians();
    let phi_sin = phi.sin();
    let q = (1.0 - e * e) * (phi_sin / (1.0 - e*e*phi_sin.powi(2)) - (1.0 / (2.0 * e) * ((1.0 - e * phi_sin) / (1.0 + e * phi_sin)).ln()));
    let beta = (q / qp).asin();

    let lam = input.lon().to_radians();
    let b = rq * (2.0 / (1.0 + beta0.sin() * beta.sin() + (beta0.cos() * beta.cos() * (lam - lam0).cos()))).sqrt();

    let east = feast + (b * d) * (beta.cos() * (lam - lam0).sin());
    let north = fnorth + (b / d) * (beta0.cos() * beta.sin() - (beta0.sin() * beta.cos() * (lam - lam0).cos()));
    Point3035::new(east, north)
}

#[inline(always)]
pub fn backward_gen<C: CrsLaea>(input: Point<C>) -> Point4326 {
    let params = C::PARAMS;
    let ellips = C::Ellipsoid::PARAMS;

    let phi0 = params.center_lat.to_radians();
    let lam0 = params.center_lon.to_radians();
    let feast = params.false_easting;
    let fnorth = params.false_northing;
    let a = ellips.semi_major_axis;
    let e = ellips.eccentricity();

    let qp = (1.0 - e * e) * ((1.0 / (1.0 - e * e)) - ((1.0 / (2.0 * e)) * ((1.0 - e) / (1.0 + e)).ln()));
    let q0 = (1.0 - e * e) * (phi0.sin() / (1.0 - e*e*phi0.sin().powi(2)) - (1.0 / (2.0 * e) * ((1.0 - e * phi0.sin()) / (1.0 + e * phi0.sin())).ln()));

    let beta0 = (q0 / qp).asin();
    let rq = a * (0.5 * qp).sqrt();

    let d = a * (phi0.cos() / (1.0 - e*e * phi0.sin().powi(2)).sqrt()) / (rq * beta0.cos());

    let east = input.coords.0;
    let north = input.coords.1;

    let p = (((east - feast)/d).powi(2) + (d * (north - fnorth)).powi(2)).sqrt();
    let c = 2.0 * (p / (2.0 * rq)).asin();
    let beta2 = ((c.cos() * beta0.sin()) + ((d * (north - fnorth) * c.sin() * beta0.cos()) / p)).asin();

    let phi = beta2 + ((e*e/3.0 + 31.0*e.powi(4)/180.0 + 517.0*e.powi(6)/5040.0) * (2.0 * beta2).sin()) + ((23.0 * e.powi(4)/360.0 + 251.0*e.powi(6)/3780.0) * (4.0 * beta2).sin()) + ((761.0 * e.powi(6) / 45360.0) * (6.0 * beta2).sin());
    let lam = lam0 + ((east - feast) * c.sin()).atan2(d*p*beta0.cos() * c.cos() - d*d*(north - fnorth) * beta0.sin()*c.sin());

    Point4326::new(phi.to_degrees(), lam.to_degrees())
}

pub fn backward_ref(input: Point3035) -> Point4326 {
    backward_gen(input)
}

pub fn backward(input: Point3035) -> Point4326 {
    // constants for EPSG:3035
    let phi0 = 0.907571211_f64; // == 52°
    let lam0 = 0.174532925_f64; // == 10°
    let feast = 4321000.0_f64;
    let fnorth = 3210000.0_f64;
    // GRS 1980 ellipsoid
    let a = 6378137.0_f64;
    let e = 0.081819191_f64;

    let qp = (1.0 - e * e) * ((1.0 / (1.0 - e * e)) - ((1.0 / (2.0 * e)) * ((1.0 - e) / (1.0 + e)).ln()));
    let q0 = (1.0 - e * e) * (phi0.sin() / (1.0 - e*e*phi0.sin().powi(2)) - (1.0 / (2.0 * e) * ((1.0 - e * phi0.sin()) / (1.0 + e * phi0.sin())).ln()));

    let beta0 = (q0 / qp).asin();
    let rq = a * (0.5 * qp).sqrt();

    let d = a * (phi0.cos() / (1.0 - e*e * phi0.sin().powi(2)).sqrt()) / (rq * beta0.cos());

    let p = (((input.east() - feast)/d).powi(2) + (d * (input.north() - fnorth)).powi(2)).sqrt();
    let c = 2.0 * (p / (2.0 * rq)).asin();
    let beta2 = ((c.cos() * beta0.sin()) + ((d * (input.north() - fnorth) * c.sin() * beta0.cos()) / p)).asin();

    let phi = beta2 + ((e*e/3.0 + 31.0*e.powi(4)/180.0 + 517.0*e.powi(6)/5040.0) * (2.0 * beta2).sin()) + ((23.0 * e.powi(4)/360.0 + 251.0*e.powi(6)/3780.0) * (4.0 * beta2).sin()) + ((761.0 * e.powi(6) / 45360.0) * (6.0 * beta2).sin());
    let lam = lam0 + ((input.east() - feast) * c.sin()).atan2(d*p*beta0.cos() * c.cos() - d*d*(input.north() - fnorth) * beta0.sin()*c.sin());

    Point4326::new(phi.to_degrees(), lam.to_degrees())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::assert_approx_eq;

    #[test]
    fn test_forward() {
        let input = Point4326::new(50.0, 5.0);
        let output = forward(input);

        assert_approx_eq(output.east(), 3962799.45, 0.01);
        assert_approx_eq(output.north(), 2999718.85, 0.01);
    }

    #[test]
    fn test_backward() {
        let input = Point3035::new(3962799.45, 2999718.85);
        let output = backward(input);

        assert_approx_eq(output.lat(), 50.0, 1e-6);
        assert_approx_eq(output.lon(), 5.0, 1e-6);
    }

    #[test]
    fn test_roundtrip() {
        let roundtrip = |lat, lon| {
            let input = Point4326::new(lat, lon);
            let output = backward(forward(input));
            assert_approx_eq(input.lat(), output.lat(), 1e-6);
            assert_approx_eq(input.lon(), output.lon(), 1e-6);
        };
        roundtrip(52.0, 10.0);
        roundtrip(54.234, 9.123);
        roundtrip(47.234, 6.123);
        roundtrip(56.234, 14.123);
    }

    #[test]
    fn test_compare() {
        let compare = |lat, lon| {
            let a = forward_ref(Point4326::new(lat, lon));
            let b = forward(Point4326::new(lat, lon));
            assert_approx_eq(a.east(), b.east(), 1e-3);
            assert_approx_eq(a.north(), b.north(), 1e-3);
        };
        compare(52.0, 10.0);
        compare(54.234, 9.123);
        compare(47.234, 6.123);
        compare(56.234, 14.123);
    }

    #[test]
    fn test_compare_backward() {
        let compare = |east, north| {
            let a = backward_ref(Point3035::new(east, north));
            let b = backward(Point3035::new(east, north));
            assert_approx_eq(a.coords.0, b.coords.0, 1e-6);
            assert_approx_eq(a.coords.1, b.coords.1, 1e-6);
        };
        compare(3962799.45, 2999718.85);
        compare(3963799.45, 2998718.85);
    }
}

