pub fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
    (a-b).abs() < eps
}

pub fn assert_approx_eq(a: f64, b: f64, eps: f64) {
    dbg!((a, b));
    assert!(approx_eq(a, b, eps));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_approx_eq() {
        assert_eq!(approx_eq(1.0, 1.0, 1e-5), true);
        assert_eq!(approx_eq(-1.0, -1.0, 1e-5), true);
        assert_eq!(approx_eq(0.0, 1e-6, 1e-5), true);
        assert_eq!(approx_eq(0.0, -1e-6, 1e-5), true);
        assert_eq!(approx_eq(2.0, 2.1, 1e-5), false);
        assert_eq!(approx_eq(2.0, 2.1, 0.5), true);
        assert_approx_eq(2.0, 2.1, 0.5);
    }
}
