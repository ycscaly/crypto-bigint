use crate::modular::runtime_mod::DynResidue;
use crate::{Limb, Pow, Uint, U128, U192, U256, U320};

macro_rules! impl_pow_cross_sizes_multiple {
    (($first_type:ident, $first_bits:expr), ($(($second_type:ident, $second_bits:expr)),+ $(,)?)) => {
        $(
            impl Pow<$second_type> for DynResidue<{nlimbs!($first_bits)}> {
                fn pow(&self, exponent: &$second_type) -> DynResidue<{nlimbs!($first_bits)}> {
                    let mut i = 0;
                    let mut shifted_base = self.clone();
                    let mut res = DynResidue::<{nlimbs!($first_bits)}>::one(self.params().clone());

                    while i < $second_bits / $first_bits {
                        let mut limbs = [Limb::ZERO; nlimbs!($first_bits)];
                        let mut j = 0;

                        while j < nlimbs!($first_bits) {
                            limbs[j] = exponent.limbs[i*nlimbs!($first_bits) + j];
                            j += 1;
                        }

                        let part: $first_type = Uint { limbs };

                        res *= shifted_base.pow(&part);
                        i += 1;

                        if i == $second_bits / $first_bits {
                            break;
                        } else {
                            shifted_base = shifted_base.pow(&$first_type::MAX) * shifted_base; // Shift by Self::bits
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
                    let mut shifted_base = self.clone();
                    let mut res = DynResidue::<{nlimbs!($first_bits)}>::one(self.params().clone());

                    while i < $second_bits / $first_bits {
                        let mut limbs = [Limb::ZERO; nlimbs!($first_bits)];
                        let mut j = 0;

                        while j < nlimbs!($first_bits) {
                            limbs[j] = exponent.limbs[i*nlimbs!($first_bits) + j];
                            j += 1;
                        }

                        let exponent_i: $first_type = Uint { limbs }; // The ith part of the exponent
                        res *= shifted_base.pow(&exponent_i);
                        shifted_base = shifted_base.pow(&$first_type::MAX) * shifted_base; // Shift by Self::bits

                        i += 1;
                    }

                    let mut limbs = [Limb::ZERO; nlimbs!($second_bits % $first_bits)];
                    let mut j = 0;

                    while j < nlimbs!($second_bits % $first_bits) {
                        limbs[j] = exponent.limbs[i*nlimbs!($first_bits) + j];
                        j += 1;
                    }

                    let part: $first_type = Uint { limbs }.into();

                    res *= shifted_base.pow_bounded_exp(
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

#[cfg(test)]
mod tests {
    use crate::modular::runtime_mod::{DynResidue, DynResidueParams};
    use crate::{Pow, U128, U192, U256, U320};

    #[test]
    fn pow_zero_and_one_cross_sizes() {
        let params = DynResidueParams::new(&U128::MAX);
        let zero = DynResidue::zero(params);
        let one = DynResidue::one(params);

        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U192>>::pow(&zero, &U192::ZERO).retrieve(),
            U128::ONE
        );
        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U256>>::pow(&zero, &U256::ZERO).retrieve(),
            U128::ONE
        );
        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U320>>::pow(&zero, &U320::ZERO).retrieve(),
            U128::ONE
        );
        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U192>>::pow(&zero, &U192::ONE).retrieve(),
            U128::ZERO
        );
        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U256>>::pow(&zero, &U256::ONE).retrieve(),
            U128::ZERO
        );
        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U320>>::pow(&zero, &U320::ONE).retrieve(),
            U128::ZERO
        );
        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U192>>::pow(&one, &U192::ZERO).retrieve(),
            U128::ONE
        );
        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U256>>::pow(&one, &U256::ZERO).retrieve(),
            U128::ONE
        );
        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U320>>::pow(&one, &U320::ZERO).retrieve(),
            U128::ONE
        );
        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U192>>::pow(&one, &U192::ONE).retrieve(),
            U128::ONE
        );
        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U256>>::pow(&one, &U256::ONE).retrieve(),
            U128::ONE
        );
        assert_eq!(
            <DynResidue::<{ U128::LIMBS }> as Pow<U320>>::pow(&one, &U320::ONE).retrieve(),
            U128::ONE
        );
    }
}
