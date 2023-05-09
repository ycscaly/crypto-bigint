use crate::{
    Uint, U1024, U128, U1280, U1536, U1792, U192, U2048, U256, U3072, U320, U3584, U384, U4096,
    U448, U512, U576, U6144, U64, U640, U768, U8192, U896,
};

use crate::Limb;

macro_rules! convert_internal {
    ($x:ident, $source_bits:expr, $target_bits:expr) => {{
        let mut limbs = [Limb::ZERO; nlimbs!($target_bits)];
        let mut i = 0;

        while i < nlimbs!($source_bits) {
            limbs[i] = $x.limbs[i];
            i += 1;
        }

        Uint { limbs }
    }};
}

macro_rules! impl_convert {
    (($source_type:ident, $source_bits:expr), ($(($target_type:ident, $target_bits:expr)),+ $(,)?)) => {
        $(
            impl From<$source_type> for $target_type {
                fn from(x: $source_type) -> Self {
                    convert_internal!(x, $source_bits, $target_bits)
                }
            }

            impl From<&$source_type> for $target_type {
                fn from(x: &$source_type) -> Self {
                    convert_internal!(x, $source_bits, $target_bits)
                }
            }
        )+
    };
}

impl_convert! {
    (U64, 64),
    (
    (U128, 128),
    (U192, 192),
    (U256, 256),
    (U320, 320),
    (U384, 384),
    (U448, 448),
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
    (U6144, 6144),
    (U8192, 8192)
    )
}

impl_convert! {
    (U128, 128),
    (
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
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U192, 192),
    (
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
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U256, 256),
    (
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
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U320, 320),
    (
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
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U384, 384),
    (
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
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U448, 448),
    (
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
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U512, 512),
    (
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
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
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
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U640, 640),
    (
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
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U768, 768),
    (
        (U896, 896),
        (U1024, 1024),
        (U1280, 1280),
        (U1536, 1536),
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U896, 896),
    (
        (U1024, 1024),
        (U1280, 1280),
        (U1536, 1536),
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U1024, 1024),
    (
        (U1280, 1280),
        (U1536, 1536),
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U1280, 1280),
    (
        (U1536, 1536),
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U1536, 1536),
    (
        (U1792, 1792),
        (U2048, 2048),
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U2048, 2048),
    (
        (U3072, 3072),
        (U3584, 3584),
        (U4096, 4096),
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U3072, 3072),
    (
        (U3584, 3584),
        (U4096, 4096),
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U3584, 3584),
    (
        (U4096, 4096),
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U4096, 4096),
    (
        (U6144, 6144),
        (U8192, 8192)
    )
}

impl_convert! {
    (U6144, 6144),
    (
        (U8192, 8192)
    )
}

#[cfg(test)]
mod tests {
    use crate::{U128, U384, U64};

    #[test]
    fn concat_zero_equals_convert() {
        let x = U64::from_u64(0x0011223344556677);

        assert_eq!(U128::from(&x), U128::from_u64(0x0011223344556677));
        assert_eq!(U128::from(x), U128::from_u64(0x0011223344556677));
    }

    #[test]
    fn converts_non_concatable_types() {
        assert_eq!(U384::ONE, U384::from(&U128::ONE));
        assert_eq!(U384::ONE, U384::from(U128::ONE));
    }
}
