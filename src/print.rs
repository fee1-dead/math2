use std::fmt::{self, Display};

use num::Signed;

use crate::factorization::SquareFreeFactorization;
use crate::traits::CommutativeRing;
use crate::Polynomial;

pub trait PrintableCoeff: Display + CommutativeRing + PartialEq + Signed {}

impl<X: Display + CommutativeRing + PartialEq + Signed> PrintableCoeff for X {}

fn print_if_not_one(x: &impl PrintableCoeff, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    if !x.is_one() {
        write!(f, "{x}")?;
    }
    Ok(())
}

fn print_as_factor(x: &impl PrintableCoeff, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    if x.is_negative() {
        f.write_str("-")?;
    }
    print_if_not_one(&x.abs(), f)
}

impl<T: PrintableCoeff> Polynomial<T> {
    pub fn print_with_var<'a>(&'a self, var: &'a str) -> PrintWithVar<'a, Polynomial<T>> {
        PrintWithVar {
            var: var.into(),
            thing: self,
        }
    }
}

impl<T: PrintableCoeff> SquareFreeFactorization<T> {
    pub fn print_with_var<'a>(
        &'a self,
        var: &'a str,
    ) -> PrintWithVar<'a, SquareFreeFactorization<T>> {
        PrintWithVar {
            var: var.into(),
            thing: self,
        }
    }
}

pub struct PrintWithVar<'a, F> {
    var: &'a str,
    thing: &'a F,
}

impl<T: PrintableCoeff> Display for PrintWithVar<'_, Polynomial<T>> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut first = true;
        let Self { var, .. } = self;
        for (degree, coeff) in self.thing.coeffs.iter().enumerate().rev() {
            if coeff.is_zero() {
                continue;
            }

            if !first {
                f.write_str(if coeff.is_negative() { " - " } else { " + " })?;
            }

            let abs;
            let coeff = if first {
                coeff
            } else {
                abs = coeff.abs();
                &abs
            };

            first = false;

            if degree == 0 {
                write!(f, "{coeff}")?;
            } else {
                print_if_not_one(coeff, f)?;
                if degree == 1 {
                    write!(f, "{var}")?;
                } else {
                    write!(f, "{var}^{degree}")?;
                }
            }
        }

        Ok(())
    }
}

impl<T: PrintableCoeff> Display for PrintWithVar<'_, SquareFreeFactorization<T>> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        print_as_factor(&self.thing.leading_coeff, f)?;
        for (poly, exp) in &self.thing.factors {
            write!(f, "({})", poly.print_with_var(self.var))?;
            if exp.get() > 1 {
                write!(f, "^{}", exp)?;
            }
        }
        Ok(())
    }
}
