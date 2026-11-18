//! Outward-rounded finite arithmetic for planar coverage certificates.

use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug)]
pub(super) struct Interval {
    pub lo: f64,
    pub hi: f64,
}

impl Interval {
    pub fn point(value: f64) -> Self {
        Self {
            lo: value,
            hi: value,
        }
    }

    pub fn abs_upper(self) -> f64 {
        self.lo.abs().max(self.hi.abs())
    }

    pub fn radius_from(self, value: f64) -> f64 {
        (self.lo - value)
            .abs()
            .max((self.hi - value).abs())
            .next_up()
    }
}

impl Add for Interval {
    type Output = Self;
    fn add(self, other: Self) -> Self {
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
        if other.lo <= 0.0 && other.hi >= 0.0 {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
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
