use std::{ops::{Div, Neg, Mul}, mem::{take, swap}};

use num::{BigInt, BigRational, One, Zero, Signed, Integer};

use crate::Polynomial;

/// Performs checked inverse.
pub trait CheckedInv {
    fn checked_inv(&self) -> Option<Self>
    where
        Self: Sized;
}

/// Asserts this value is a unit (invertible element) in a ring.
pub struct AssertUnit<T>(T);

impl<T> AssertUnit<T> {
    pub fn into_inner(self) -> T {
        self.0
    }
}

pub trait CommutativeRing: Zero + One + Neg<Output = Self> + Clone {
    fn sub(self, other: Self) -> Self {
        self.add(other.neg())
    }
    fn invert(x: &AssertUnit<Self>) -> AssertUnit<Self>;
    fn assert_is_unit(self) -> AssertUnit<Self>;
    /// Returns whether this element is nilpotent (i.e. there exists some n such that x^n == 0).
    /// Reduced rings has no non-zero nilpotent elements. Integral domains are an example of reduced rings.
    fn is_nilpotent(&self) -> bool;
}

pub trait CoefficientDomain: CommutativeRing{
    fn unit_and_normal(self) -> (AssertUnit<Self>, Self);
    fn gcd(&self, other: &Self) -> Self;
}

pub trait FromUsize {
    fn from_usize(n: usize) -> Self;
}

impl FromUsize for BigInt {
    fn from_usize(n: usize) -> Self {
        Self::from(n)
    }
}

impl FromUsize for BigRational {
    fn from_usize(n: usize) -> Self {
        Self::from(BigInt::from(n))
    }
}

/// any implementors of this trait have their set of field elements represented
/// by the possible values the implementor type can take.
pub trait Field: CommutativeRing + CheckedInv + Div<Output = Self> {
    fn div(self, other: Self) -> Option<Self> {
        other.checked_inv().map(|b| self.mul(b))
    }
}

impl CheckedInv for BigRational {
    fn checked_inv(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(self.recip())
        }
    }
}

/// The ring of integers (`Z`)
impl CommutativeRing for BigInt {
    fn assert_is_unit(self) -> AssertUnit<Self> {
        assert!(self.abs().is_one());
        AssertUnit(self)
    }
    fn invert(x: &AssertUnit<Self>) -> AssertUnit<Self> {
        (-&x.0).assert_is_unit()
    }
    fn is_nilpotent(&self) -> bool {
        self.is_zero()
    }
}

impl CoefficientDomain for BigInt {
    fn unit_and_normal(self) -> (AssertUnit<Self>, Self) {
        if self.is_negative() {
            ((-BigInt::one()).assert_is_unit(), -self)
        } else {
            (BigInt::one().assert_is_unit(), self)
        }
    }
    fn gcd(&self, other: &Self) -> Self {
        Integer::gcd(self, other)
    }
}

impl CommutativeRing for i64 {
    fn assert_is_unit(self) -> AssertUnit<Self> {
        assert_eq!(1, self.abs());
        AssertUnit(self)
    }
    fn invert(x: &AssertUnit<Self>) -> AssertUnit<Self> {
        (-x.0).assert_is_unit()
    }
    fn is_nilpotent(&self) -> bool {
        *self == 0
    }
}

impl CoefficientDomain for i64 {
    fn unit_and_normal(self) -> (AssertUnit<Self>, Self) {
        if self < 0 {
            ((-1).assert_is_unit(), -self)
        } else {
            (1.assert_is_unit(), self)
        }
    }
    fn gcd(&self, other: &Self) -> Self {
        Integer::gcd(self, other)
    }
}

/// The ring of rationals (`Q`)
impl CommutativeRing for BigRational {
    fn invert(x: &AssertUnit<Self>) -> AssertUnit<Self> {
        x.0.recip().assert_is_unit()
    }
    fn assert_is_unit(self) -> AssertUnit<Self> {
        assert!(!self.is_zero());
        AssertUnit(self)
    }
    fn is_nilpotent(&self) -> bool {
        self.is_zero()
    }
}

impl CoefficientDomain for BigRational {
    /// in rationals, all non-zero numbers are units, and the unit normal can
    /// either be zero or one.
    fn unit_and_normal(self) -> (AssertUnit<Self>, Self) {
        if self.is_zero() {
            (Self::one().assert_is_unit(), Self::zero())
        } else {
            (self.assert_is_unit(), Self::one())
        }
    }
    fn gcd(&self, other: &Self) -> Self {
        if self.is_zero() && other.is_zero() {
            panic!("0 gcd 0");
        }

        Self::one()
    }
}

/// The field of rationals (`Q`)
impl Field for BigRational {}

/// The ring of polynomials over a ring (`R[x]`)
impl<Ring: CommutativeRing> CommutativeRing for Polynomial<Ring> {
    fn assert_is_unit(self) -> AssertUnit<Self> {
        assert!(!self.is_zero());
        for (degree, c) in self.coeffs.iter().enumerate() {
            if degree == 0 {
                c.clone().assert_is_unit();
            } else {
                assert!(c.is_nilpotent())
            }
        }
        AssertUnit(self)
    }
    fn invert(x: &AssertUnit<Self>) -> AssertUnit<Self> {
        todo!()
    }
    fn is_nilpotent(&self) -> bool {
        todo!()
    }
}

impl<K: CoefficientDomain> CoefficientDomain for Polynomial<K> {
    fn unit_and_normal(mut self) -> (AssertUnit<Self>, Self) {
        if self.is_zero() {
            return (Self::one().assert_is_unit(), Self::zero());
        }
        let (unit, _) = self.leading_coefficient_cloned().unit_and_normal();
        for coeff in &mut self.coeffs {
            let mut c = K::zero();
            swap(coeff, &mut c);
            *coeff = c.mul(K::invert(&unit).into_inner());
        }
        (Polynomial::from_elem_with_degree(unit.into_inner(), 1).assert_is_unit(), self)
    }
    fn gcd(&self, other: &Self) -> Self {
        todo!()
    }
}