//! Outward-rounded finite arithmetic for planar coverage certificates.

use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug)]
pub(super) struct Interval {
    pub lo: f64,
    pub hi: f64,
}

impl Interval {
    pub fn point(value: f64) -> Self {
        if !value.is_finite() {
            return Self::unbounded();
        }
        Self {
            lo: value,
            hi: value,
        }
    }

    fn unbounded() -> Self {
        Self {
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        }
    }

    fn valid(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite() && self.lo <= self.hi
    }

    pub fn abs_upper(self) -> f64 {
        self.lo.abs().max(self.hi.abs())
    }

    pub fn radius_from(self, value: f64) -> f64 {
        if !self.valid() || !value.is_finite() {
            return f64::INFINITY;
        }
        if self.lo == value && self.hi == value {
            return 0.0;
        }
        (self.lo - value)
            .abs()
            .max((self.hi - value).abs())
            .next_up()
    }
}

impl Add for Interval {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        if !self.valid() || !other.valid() {
            return Self::unbounded();
        }
        if self.lo == 0.0 && self.hi == 0.0 {
            return other;
        }
        if other.lo == 0.0 && other.hi == 0.0 {
            return self;
        }
        Self {
            lo: (self.lo + other.lo).next_down(),
            hi: (self.hi + other.hi).next_up(),
        }
    }
}

impl Neg for Interval {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
}

impl Sub for Interval {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self + -other
    }
}

impl Mul for Interval {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        if !self.valid() || !other.valid() {
            return Self::unbounded();
        }
        // Exact identities avoid manufacturing an uncertainty interval for the
        // zero and signed-unit coefficients of a cardinal projection frame.
        for (a, b) in [(self, other), (other, self)] {
            if a.lo == a.hi {
                if a.lo == 1.0 {
                    return b;
                }
                if a.lo == -1.0 {
                    return -b;
                }
                if a.lo == 0.0 && b.lo.is_finite() && b.hi.is_finite() && b.lo <= b.hi {
                    return Self::point(0.0);
                }
            }
        }
        let values = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        Self {
            lo: values.into_iter().fold(f64::INFINITY, f64::min).next_down(),
            hi: values
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        }
    }
}

impl Div for Interval {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        if !self.valid() || !other.valid() {
            return Self::unbounded();
        }
        if other.lo <= 0.0 && other.hi >= 0.0 {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
        }
        if other.lo == 1.0 && other.hi == 1.0 {
            return self;
        }
        if other.lo == -1.0 && other.hi == -1.0 {
            return -self;
        }
        self * Self {
            lo: (1.0 / other.hi).next_down(),
            hi: (1.0 / other.lo).next_up(),
        }
    }
}

pub(super) fn dot(a: [Interval; 3], b: [Interval; 3]) -> Interval {
    (0..3).fold(Interval::point(0.0), |sum, k| sum + a[k] * b[k])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_zero_and_signed_unit_identities_keep_invalid_intervals_invalid() {
        let a = Interval { lo: -0.7, hi: 1.3 };
        for x in [
            a + Interval::point(0.0),
            a * Interval::point(1.0),
            a / Interval::point(1.0),
        ] {
            assert_eq!((x.lo, x.hi), (a.lo, a.hi));
        }
        let zero = a * Interval::point(0.0);
        assert_eq!((zero.lo, zero.hi), (0.0, 0.0));
        assert_eq!(zero.radius_from(0.0), 0.0);
        let invalid = Interval {
            lo: f64::INFINITY,
            hi: f64::INFINITY,
        };
        let product = invalid * Interval::point(0.0);
        assert!(!product.lo.is_finite() || !product.hi.is_finite());
        let quotient = Interval::point(0.0) / Interval { lo: -1.0, hi: 1.0 };
        assert!(!quotient.lo.is_finite() && !quotient.hi.is_finite());
    }
}
