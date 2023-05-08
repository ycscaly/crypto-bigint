use crate::{
    Uint, U1024, U128, U1280, U1536, U1792, U192, U2048, U256, U3072, U320, U3584, U384, U4096,
    U448, U512, U576, U6144, U64, U640, U768, U8192, U896,
};

use crate::{Pow, PowBoundedExp};

macro_rules! impl_pow_cross_sizes {
    (($first_type:ident, $first_bits:expr), ($(($second_type:ident, $second_bits:expr)),+ $(,)?)) => {
        $(
            impl Pow<$second_type> for DynResidue<{nlimbs!($first_bits)}> {
                fn pow(&self, exponent: &$second_type) -> DynResidue<{nlimbs!($first_bits)}> {
                    let d = $second_bits / $first_bits;
                    let r = $second_bits % $first_bits;

                    let parts = if r > 0 {
                        d + 1
                    } else {
                        d
                    };

                    let mut i = 0;
                    let mut res = DynResidue<{nlimbs!($first_bits)}>::one(self.params());
                    while i < d {
                        let mut limbs = [Limb::ZERO; nlimbs!($first_bits)];
                        let mut j = 0;

                        while j < nlimbs!($first_bits) && j < (nlimbs!($second_bits) - i*nlimbs!($first_bits)) {
                            limbs[j] = exponent.limbs[i*nlimbs!($first_bits) + j];
                            j += 1;
                        }

                        let part: $first_type = Uint { limbs };

                        let x = self.pow(&part);
                        let res = (res * (h_hat.pow(&U4096::MAX)) * h_hat).retrieve(); // TODO: *=

                        i += 1;
                    }
               }
            }

            // impl Mul<&$second_type> for $first_type  {
            //     type Output = Uint<{nlimbs!($first_bits) + nlimbs!($second_bits)}>;
            //
            //     fn mul(self, rhs: &$second_type) -> Uint<{nlimbs!($first_bits) + nlimbs!($second_bits)}> {
            //         self.mul_wide(rhs).into()
            //     }
            // }
            //
            // impl Mul<$second_type> for &$first_type  {
            //     type Output = Uint<{nlimbs!($first_bits) + nlimbs!($second_bits)}>;
            //
            //     fn mul(self, rhs: $second_type) -> Uint<{nlimbs!($first_bits) + nlimbs!($second_bits)}> {
            //         self.mul_wide(&rhs).into()
            //     }
            // }
            //
            // impl Mul<&$second_type> for &$first_type {
            //     type Output = Uint<{nlimbs!($first_bits) + nlimbs!($second_bits)}>;
            //
            //     fn mul(self, rhs: &$second_type) -> Uint<{nlimbs!($first_bits) + nlimbs!($second_bits)}> {
            //         self.mul_wide(rhs).into()
            //     }
            // }
            //
            // impl Mul<$first_type> for $second_type {
            //     type Output = Uint<{nlimbs!($second_bits) + nlimbs!($first_bits)}>;
            //
            //     fn mul(self, rhs: $first_type) -> Uint<{nlimbs!($second_bits) + nlimbs!($first_bits)}> {
            //         self.mul_wide(&rhs).into()
            //     }
            // }
            //
            // impl Mul<&$first_type> for $second_type  {
            //     type Output = Uint<{nlimbs!($second_bits) + nlimbs!($first_bits)}>;
            //
            //     fn mul(self, rhs: &$first_type) -> Uint<{nlimbs!($second_bits) + nlimbs!($first_bits)}> {
            //         self.mul_wide(rhs).into()
            //     }
            // }
            //
            // impl Mul<$first_type> for &$second_type  {
            //     type Output = Uint<{nlimbs!($second_bits) + nlimbs!($first_bits)}>;
            //
            //     fn mul(self, rhs: $first_type) -> Uint<{nlimbs!($second_bits) + nlimbs!($first_bits)}> {
            //         self.mul_wide(&rhs).into()
            //     }
            // }
            //
            // impl Mul<&$first_type> for &$second_type {
            //     type Output = Uint<{nlimbs!($second_bits) + nlimbs!($first_bits)}>;
            //
            //     fn mul(self, rhs: &$first_type) -> Uint<{nlimbs!($second_bits) + nlimbs!($first_bits)}> {
            //         self.mul_wide(rhs).into()
            //     }
            // }
        )+
    };
}
