use std::iter::Product;
use std::ops::Mul;
use std::ops::BitXor;

pub trait Group: Sized + Mul + Mul<Output=Self> + Product + BitXor<i64> {
    fn identity() -> Self;
    fn inverse_of_gen(self) -> Self;

    // When customarily implementing the 'inverse(),' one needs to be sure that the inverse function is consistent with the 'inverse_of_gen()'
    fn inverse(self) -> Self {
        let gens = self.in_terms_of_gens();
        gens.into_iter().rev().map(|x| x.inverse_of_gen()).product()
    }

    // in_terms_of_gen can be implemented if there is a nice set of generators
    fn in_terms_of_gens(self) -> Vec<Self> {
        vec![self]
    }
}

pub mod permutation;

// Group action
pub trait Gset<G: Group>: Sized {
    // In this function, we assume that 'elem' is a generator of the group
    fn gen_action_by(&mut self, elem: G);

    fn action_by(&mut self, elem: G) {
        let gens = elem.in_terms_of_gens();
        for gen in gens.into_iter() {
            self.gen_action_by( gen )
        };
    }
}

