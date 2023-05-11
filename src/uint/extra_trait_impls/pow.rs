use crate::modular::runtime_mod::DynResidue;
use crate::{Limb, Pow, Uint, U128, U192, U256, U320};

macro_rules! impl_pow_cross_sizes_multiple {
    (($first_type:ident, $first_bits:expr), ($(($second_type:ident, $second_bits:expr)),+ $(,)?)) => {
        $(
            impl Pow<$second_type> for DynResidue<{nlimbs!($first_bits)}> {
                fn pow(&self, exponent: &$second_type) -> DynResidue<{nlimbs!($first_bits)}> {
                    let mut i = 0;
                    let mut shifted_self = self.clone();
                    let mut res = DynResidue::<{nlimbs!($first_bits)}>::one(self.params().clone());

                    while i < $second_bits / $first_bits {
                        let mut limbs = [Limb::ZERO; nlimbs!($first_bits)];
                        let mut j = 0;

                        while j < nlimbs!($first_bits) {
                            limbs[j] = exponent.limbs[i*nlimbs!($first_bits) + j];
                            j += 1;
                        }

                        let part: $first_type = Uint { limbs };

                        res *= shifted_self.pow(&part);
                        i += 1;

                        if i == $second_bits / $first_bits {
                            break;
                        } else {
                            shifted_self = shifted_self.pow(&$first_type::MAX) * shifted_self; // Shift by Self::bits
                        }
                    }

                    res
               }
            }
        )+
    };
}

macro_rules! impl_pow_cross_sizes {
    (($first_type:ident, $first_bits:expr), ($(($second_type:ident, $second_bits:expr)),+ $(,)?)) => {
        $(
            impl Pow<$second_type> for DynResidue<{nlimbs!($first_bits)}> {
                fn pow(&self, exponent: &$second_type) -> DynResidue<{nlimbs!($first_bits)}> {
                    let mut i = 0;
                    let mut shifted_self = self.clone();
                    let mut res = DynResidue::<{nlimbs!($first_bits)}>::one(self.params().clone());

                    while i < $second_bits / $first_bits {
                        let mut limbs = [Limb::ZERO; nlimbs!($first_bits)];
                        let mut j = 0;

                        while j < nlimbs!($first_bits) {
                            limbs[j] = exponent.limbs[i*nlimbs!($first_bits) + j];
                            j += 1;
                        }

                        let part: $first_type = Uint { limbs };

                        res *= shifted_self.pow(&part);
                        shifted_self = shifted_self.pow(&$first_type::MAX) * shifted_self; // Shift by Self::bits

                        i += 1;
                    }

                    let mut limbs = [Limb::ZERO; nlimbs!($second_bits % $first_bits)];
                    let mut j = 0;

                    while j < nlimbs!($second_bits % $first_bits) {
                        limbs[j] = exponent.limbs[i*nlimbs!($first_bits) + j];
                        j += 1;
                    }

                    let part: $first_type = Uint { limbs }.into();

                    res *= shifted_self.pow_bounded_exp(
                        &part,
                        $second_bits % $first_bits,
                    );

                    res
               }
            }
        )+
    };
}

impl_pow_cross_sizes_multiple! {
    (U128, 128),
    (
        (U256, 256),
    )
}

impl_pow_cross_sizes! {
    (U128, 128),
    (
        (U192, 192),
        (U320, 320)
    )
}
