use crate::modular::runtime_mod::DynResidue;

use crate::{
    Limb, Pow, Uint, U1024, U128, U1280, U1536, U1792, U192, U2048, U256, U3072, U320, U3584, U384,
    U4096, U4224, U4352, U448, U512, U576, U6144, U64, U640, U768, U8192, U896,
};

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
    (U64, 64),
    (
        (U128, 128),
        (U192, 192),
        (U256, 256),
        (U320, 320),
        (U384, 384),
        (U448, 448),
        (U512, 512),
        (U576, 576),
        (U640, 640),
        (U768, 768),
        (U896, 896),
        (U1024, 1024),
        (U1280, 1280),
        (U1536, 1536),
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes_multiple! {
    (U128, 128),
    (
        (U256, 256),
        (U384, 384),
        (U512, 512),
        (U640, 640),
        (U768, 768),
        (U896, 896),
        (U1024, 1024),
        (U1280, 1280),
        (U1536, 1536),
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes_multiple! {
    (U192, 192),
    (
        (U384, 384),
        (U576, 576),
        (U768, 768),
        (U1536, 1536),
        (U3072, 3072),
        (U4224, 4224),
        (U6144, 6144),
    )
}

impl_pow_cross_sizes_multiple! {
    (U256, 256),
    (
        (U512, 512),
        (U768, 768),
        (U1024, 1024),
        (U1280, 1280),
        (U1536, 1536),
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes_multiple! {
    (U320, 320),
    (
        (U640, 640),
        (U1280, 1280),
    )
}

impl_pow_cross_sizes_multiple! {
    (U384, 384),
    (
        (U768, 768),
        (U1536, 1536),
        (U3072, 3072),
        (U4224, 4224),
        (U6144, 6144),
    )
}

impl_pow_cross_sizes_multiple! {
    (U448, 448),
    (
        (U896, 896),
        (U1792, 1792),
        (U3584, 3584),
    )
}

impl_pow_cross_sizes_multiple! {
    (U512, 512),
    (
        (U1024, 1024),
        (U1536, 1536),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes_multiple! {
    (U640, 640),
    (
        (U1280, 1280),
    )
}

impl_pow_cross_sizes_multiple! {
    (U768, 768),
    (
        (U1536, 1536),
        (U3072, 3072),
        (U6144, 6144),
    )
}

impl_pow_cross_sizes_multiple! {
    (U896, 896),
    (
        (U1792, 1792),
        (U3584, 3584),
    )
}

impl_pow_cross_sizes_multiple! {
    (U1024, 1024),
    (
        (U2048, 2048),
        (U3072, 3072),
        (U4096, 4096),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes_multiple! {
    (U1536, 1536),
    (
        (U3072, 3072),
        (U6144, 6144),
    )
}

impl_pow_cross_sizes_multiple! {
    (U1792, 1792),
    (
        (U3584, 3584),
    )
}

impl_pow_cross_sizes_multiple! {
    (U2048, 2048),
    (
        (U4096, 4096),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes_multiple! {
    (U3072, 3072),
    (
        (U6144, 6144),
    )
}

impl_pow_cross_sizes_multiple! {
    (U4096, 4096),
    (
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U128, 128),
    (
        (U192, 192),
        (U320, 320),
        (U448, 448),
        (U576, 576),
    )
}

impl_pow_cross_sizes! {
    (U192, 192),
    (
        (U256, 256),
        (U320, 320),
        (U448, 448),
        (U512, 512),
        (U640, 640),
        (U896, 896),
        (U1024, 1024),
        (U1280, 1280),
        (U1792, 1792),
        (U2048, 2048),
        (U3584, 3584),
        (U4096, 4096),
        (U4352, 4352),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U256, 256),
    (
        (U320, 320),
        (U384, 384),
        (U448, 448),
        (U576, 576),
        (U640, 640),
        (U896, 896),
        (U4224, 4224),
    )
}

impl_pow_cross_sizes! {
    (U320, 320),
    (
        (U384, 384),
        (U448, 448),
        (U512, 512),
        (U576, 576),
        (U768, 768),
        (U896, 896),
        (U1024, 1024),
        (U1536, 1536),
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U384, 384),
    (
        (U448, 448),
        (U512, 512),
        (U576, 576),
        (U640, 640),
        (U896, 896),
        (U1024, 1024),
        (U1280, 1280),
        (U1792, 1792),
        (U2048, 2048),
        (U3584, 3584),
        (U4096, 4096),
        (U4352, 4352),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U448, 448),
    (
        (U512, 512),
        (U576, 576),
        (U640, 640),
        (U768, 768),
        (U1024, 1024),
        (U1280, 1280),
        (U1536, 1536),
        (U2048, 2048),
        (U3072, 3072),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U512, 512),
    (
        (U576, 576),
        (U640, 640),
        (U768, 768),
        (U896, 896),
        (U1280, 1280),
        (U1792, 1792),
        (U4224, 4224),
        (U4352, 4352),
    )
}

impl_pow_cross_sizes! {
    (U576, 576),
    (
        (U640, 640),
        (U768, 768),
        (U896, 896),
        (U1024, 1024),
        (U1280, 1280),
        (U1536, 1536),
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U640, 640),
    (
        (U768, 768),
        (U896, 896),
        (U1024, 1024),
        (U1536, 1536),
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U768, 768),
    (
        (U896, 896),
        (U1024, 1024),
        (U1280, 1280),
        (U1792, 1792),
        (U2048, 2048),
        (U3584, 3584),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U896, 896),
    (
        (U1024, 1024),
        (U1280, 1280),
        (U1536, 1536),
        (U2048, 2048),
        (U3072, 3072),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U1024, 1024),
    (
        (U1280, 1280),
        (U1536, 1536),
        (U1792, 1792),
        (U3584, 3584),
        (U4224, 4224),
        (U4352, 4352),
    )
}

impl_pow_cross_sizes! {
    (U1280, 1280),
    (
        (U1536, 1536),
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U1536, 1536),
    (
        (U1792, 1792),
        (U2048, 2048),
        (U3584, 3584),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U1792, 1792),
    (
        (U2048, 2048),
        (U3072, 3072),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U2048, 2048),
    (
        (U3072, 3072),
        (U3584, 3584),
        (U4224, 4224),
        (U4352, 4352),
    )
}

impl_pow_cross_sizes! {
    (U3072, 3072),
    (
        (U3584, 3584),
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U3584, 3584),
    (
        (U4096, 4096),
        (U4224, 4224),
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U4096, 4096),
    (
        (U4224, 4224),
        (U4352, 4352),
        (U6144, 6144),
    )
}

impl_pow_cross_sizes! {
    (U4224, 4224),
    (
        (U4352, 4352),
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U4352, 4352),
    (
        (U6144, 6144),
        (U8192, 8192),
    )
}

impl_pow_cross_sizes! {
    (U6144, 6144),
    (
        (U8192, 8192),
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
