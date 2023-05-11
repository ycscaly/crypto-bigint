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
                    let mut i = 0;
                    let mut shifted_self = self;
                    let mut res = DynResidue<{nlimbs!($first_bits)}>::one(self.params());

                    while i < $second_bits / $first_bits {
                        let mut limbs = [Limb::ZERO; nlimbs!($first_bits)];
                        let mut j = 0;

                        while j < nlimbs!($first_bits) {
                            limbs[j] = exponent.limbs[i*nlimbs!($first_bits) + j];
                            j += 1;
                        }

                        let part: $first_type = Uint { limbs };

                        let x = self.pow(&part);
                        let res = (res * (h_hat.pow(&U4096::MAX)) * h_hat).retrieve(); // TODO: *=

                        res *= shifted_self.pow(&part);
                        shifted_self = shifted_self.pow(&Self::MAX) * shifted_self; // Shift by Self::bits

                        i += 1;
                    }

                    if $second_bits % $first_bits > 0 {
                        let mut limbs = [Limb::ZERO; nlimbs!($second_bits % $first_bits)];
                        let mut j = 0;

                        while j < nlimbs!($second_bits % $first_bits) {
                            limbs[j] = exponent.limbs[i*nlimbs!($first_bits) + j];
                            j += 1;
                        }

                        let part: $first_type = Uint { limbs };

                        let hi = shifted_self.pow_bounded_exp(
                            &part,
                            $second_bits % $first_bits,
                        );
                    }

                    res
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
