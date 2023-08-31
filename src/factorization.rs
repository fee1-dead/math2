use std::num::NonZeroUsize;

use num::{One, Zero};

use crate::Polynomial;
use crate::traits::{Field, FromUsize};

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct SquareFreeFactorization<F> {
    pub leading_coeff: F,
    pub factors: Vec<(Polynomial<F>, NonZeroUsize)>,
}


impl<F: Field> Polynomial<F> {
    pub fn square_free_factorization(self) -> SquareFreeFactorization<F>
    where
        F: FromUsize + PartialEq,
    {
        if self.is_zero() {
            return SquareFreeFactorization { leading_coeff: F::zero(), factors: Vec::new() };
        }

        let leading_coeff = self.leading_coefficient_cloned();
        let u = self.scalar_mul(leading_coeff.clone().checked_inv().unwrap());
        let mut factors = Vec::new();
        let mut r = u.clone().gcd(u.clone().derivative());
        let mut f = u.div_rem(r.clone()).0;
        let mut j = NonZeroUsize::new(1).unwrap();
        while !r.is_one() {
            let g = r.clone().gcd(f.clone());
            let s = f.div_rem(g.clone()).0;
            if !s.is_one() {
                factors.push((s, j));
            }
            r = r.div_rem(g.clone()).0;
            f = g;
            j = j.saturating_add(1);
        }
        if !f.is_one() {
            factors.push((f, j));
        }
        SquareFreeFactorization { leading_coeff, factors }
    }
}