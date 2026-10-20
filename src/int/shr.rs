//! [`Int`] bitwise right shift operations.

use crate::{Choice, CtOption, Int, Limb, ShrVartime, Uint, WrappingShr};
use core::ops::{Shr, ShrAssign};

impl<const LIMBS: usize> Int<LIMBS> {
    /// Computes `self >> shift`.
    ///
    /// Note, this is _signed_ shift right, i.e., the value shifted in on the left is equal to
    /// the most significant bit.
    ///
    /// # Panics
    /// - if `shift >= Self::BITS`.
    #[inline(always)]
    #[must_use]
    pub const fn shr(&self, shift: u32) -> Self {
        let sign_bits = Self::select(&Self::ZERO, &Self::MINUS_ONE, self.is_negative());
        let res = Uint::shr(&self.0, shift);
        Self::from_bits(res.bitor(&sign_bits.0.unbounded_shl(Self::BITS - shift)))
    }

    /// Computes `self >> shift` in variable time.
    ///
    /// Note, this is _signed_ shift right, i.e., the value shifted in on the left is equal to
    /// the most significant bit.
    ///
    /// NOTE: this operation is variable time with respect to `shift` *ONLY*.
    ///
    /// When used with a fixed `shift`, this function is constant-time with respect to `self`.
    ///
    /// # Panics
    /// - if `shift >= Self::BITS`.
    #[inline(always)]
    #[must_use]
    #[track_caller]
    pub const fn shr_vartime(&self, shift: u32) -> Self {
        self.overflowing_shr_vartime(shift)
            .expect("`shift` within the bit size of the integer")
    }

    /// Computes `self >> shift` in variable time: the `v0.7.0-rc.9` `Int::shr_vartime`, kept in
    /// this fork line for line so that downstream code that called it keeps calling the same
    /// routine. It carries the `_rc9` suffix because 0.7.5 has a `shr_vartime` of its own, with a
    /// different body (it builds the sign fill with a constant-time shift of a full-width value).
    ///
    /// Note, this is _signed_ shift right, i.e., the value shifted in on the left is equal to
    /// the most significant bit.
    ///
    /// # Panics
    /// - if `shift >= Self::BITS`.
    #[inline(always)]
    #[must_use]
    pub const fn shr_vartime_rc9(&self, shift: u32) -> Self {
        self.overflowing_shr_vartime_rc9(shift)
            .expect("`shift` within the bit size of the integer")
    }

    /// The `v0.7.0-rc.9` `Int::overflowing_shr_vartime`, line for line (its `ConstCtOption`
    /// became `CtOption`, and its `Limb::select` is the crate-private mask select of 0.7.5).
    ///
    /// NOTE: this operation is variable time with respect to `shift` *ONLY*.
    #[inline(always)]
    #[must_use]
    #[allow(clippy::integer_division_remainder_used, reason = "rc.9 body")]
    pub const fn overflowing_shr_vartime_rc9(&self, shift: u32) -> CtOption<Self> {
        let is_negative = self.is_negative();

        if shift >= Self::BITS {
            return CtOption::none(Self::select(&Self::ZERO, &Self::MINUS_ONE, is_negative));
        }

        // Select the base limb, based on the sign of this value.
        let base = Limb::select(Limb::ZERO, Limb::MAX, is_negative);
        let mut limbs = [base; LIMBS];

        let shift_num = (shift / Limb::BITS) as usize;
        let rem = shift % Limb::BITS;

        let mut i = 0;
        while i < LIMBS - shift_num {
            limbs[i] = self.0.limbs[i + shift_num];
            i += 1;
        }

        if rem == 0 {
            return CtOption::some(Self(Uint::new(limbs)));
        }

        // construct the carry s.t. the `rem`-most significant bits of `carry` are 1 when self
        // is negative, i.e., shift in 1s when the msb is 1.
        let mut carry = Limb::select(Limb::ZERO, Limb::MAX, is_negative);
        carry = carry.bitxor(carry.shr(rem)); // logical shift right; shifts in zeroes.

        while i > 0 {
            i -= 1;
            let shifted = limbs[i].shr(rem);
            let new_carry = limbs[i].shl(Limb::BITS - rem);
            limbs[i] = shifted.bitor(carry);
            carry = new_carry;
        }

        CtOption::some(Self(Uint::new(limbs)))
    }

    /// Computes `self >> shift`.
    ///
    /// Note, this is _signed_ shift right, i.e., the value shifted in on the left is equal to
    /// the most significant bit.
    ///
    /// Returns `None` if `shift >= Self::BITS`.
    #[inline(always)]
    #[must_use]
    #[allow(clippy::integer_division_remainder_used, reason = "needs triage")]
    pub const fn overflowing_shr(&self, shift: u32) -> CtOption<Self> {
        let in_range = Choice::from_u32_lt(shift, Self::BITS);
        let adj_shift = in_range.select_u32(0, shift);
        CtOption::new(self.shr(adj_shift), in_range)
    }

    /// Computes `self >> shift`.
    ///
    /// NOTE: this is _signed_ shift right, i.e., the value shifted in on the left is equal to
    /// the most significant bit.
    ///
    /// Returns `None` if `shift >= Self::BITS`.
    ///
    /// NOTE: this operation is variable time with respect to `shift` *ONLY*.
    ///
    /// When used with a fixed `shift`, this function is constant-time with respect to `self`.
    #[inline(always)]
    #[must_use]
    pub const fn overflowing_shr_vartime(&self, shift: u32) -> Option<Self> {
        if shift < Self::BITS {
            Some(self.unbounded_shr_vartime(shift))
        } else {
            None
        }
    }

    /// Computes `self >> shift` in a panic-free manner.
    ///
    /// If the shift exceeds the precision, returns
    /// - `0` when `self` is non-negative, and
    /// - `-1` when `self` is negative.
    #[inline(always)]
    #[must_use]
    pub const fn unbounded_shr(&self, shift: u32) -> Self {
        let default = Self::select(&Self::ZERO, &Self::MINUS_ONE, self.is_negative());
        ctutils::unwrap_or!(self.overflowing_shr(shift), default, Self::select)
    }

    /// Computes `self >> shift` in variable-time in a panic-free manner.
    ///
    /// If the shift exceeds the precision, returns
    /// - `0` when `self` is non-negative, and
    /// - `-1` when `self` is negative.
    ///
    /// NOTE: this operation is variable time with respect to `shift` *ONLY*.
    ///
    /// When used with a fixed `shift`, this function is constant-time with respect to `self`.
    #[inline(always)]
    #[must_use]
    pub const fn unbounded_shr_vartime(&self, shift: u32) -> Self {
        let sign_bits = Self::select(&Self::ZERO, &Self::MINUS_ONE, self.is_negative());
        if let Some(res) = self.0.overflowing_shr_vartime(shift) {
            Self::from_bits(res.bitor(&sign_bits.0.unbounded_shl(Self::BITS - shift)))
        } else {
            sign_bits
        }
    }

    /// Computes `self >> shift` in a panic-free manner.
    ///
    /// If the shift exceeds the precision, returns
    /// - `0` when `self` is non-negative, and
    /// - `-1` when `self` is negative.
    ///
    /// This is the `v0.7.0-rc.9` behaviour, kept in this fork; upstream 0.7.5 reduces the shift
    /// modulo the type's width instead.
    #[inline]
    #[must_use]
    pub const fn wrapping_shr(&self, shift: u32) -> Self {
        self.unbounded_shr(shift)
    }

    /// Computes `self >> shift` in variable-time in a panic-free manner.
    ///
    /// If the shift exceeds the precision, returns
    /// - `0` when `self` is non-negative, and
    /// - `-1` when `self` is negative.
    ///
    /// This is the `v0.7.0-rc.9` behaviour, kept in this fork; upstream 0.7.5 reduces the shift
    /// modulo the type's width instead.
    ///
    /// NOTE: this operation is variable time with respect to `shift` *ONLY*.
    ///
    /// When used with a fixed `shift`, this function is constant-time with respect to `self`.
    #[inline]
    #[must_use]
    pub const fn wrapping_shr_vartime(&self, shift: u32) -> Self {
        self.unbounded_shr_vartime(shift)
    }
}

macro_rules! impl_shr {
    ($($shift:ty),+) => {
        $(
            impl<const LIMBS: usize> Shr<$shift> for Int<LIMBS> {
                type Output = Int<LIMBS>;

                #[inline]
                fn shr(self, shift: $shift) -> Int<LIMBS> {
                    <&Self>::shr(&self, shift)
                }
            }

            impl<const LIMBS: usize> Shr<$shift> for &Int<LIMBS> {
                type Output = Int<LIMBS>;

                #[inline]
                fn shr(self, shift: $shift) -> Int<LIMBS> {
                    Int::<LIMBS>::shr(self, u32::try_from(shift).expect("invalid shift"))
                }
            }

            impl<const LIMBS: usize> ShrAssign<$shift> for Int<LIMBS> {
                fn shr_assign(&mut self, shift: $shift) {
                    *self = self.shr(shift)
                }
            }
        )+
    };
}

impl_shr!(i32, u32, usize);

impl<const LIMBS: usize> WrappingShr for Int<LIMBS> {
    fn wrapping_shr(&self, shift: u32) -> Int<LIMBS> {
        self.wrapping_shr(shift)
    }
}

impl<const LIMBS: usize> ShrVartime for Int<LIMBS> {
    fn overflowing_shr_vartime(&self, shift: u32) -> Option<Self> {
        self.overflowing_shr_vartime(shift)
    }

    fn unbounded_shr_vartime(&self, shift: u32) -> Self {
        self.unbounded_shr_vartime(shift)
    }

    fn wrapping_shr_vartime(&self, shift: u32) -> Self {
        self.wrapping_shr_vartime(shift)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn shr_vartime_rc9_agrees_with_shr_vartime() {
        use crate::I256;
        let values = [I256::ZERO, I256::ONE, I256::MINUS_ONE, I256::MAX, I256::MIN, I256::from_i64(-123_456_789)];
        for v in values {
            for shift in [0, 1, 2, 63, 64, 65, 127, 128, 200, 255] {
                assert_eq!(v.shr_vartime_rc9(shift), v.shr_vartime(shift), "{v:?} >> {shift}");
            }
        }
    }

    use core::ops::Div;

    use crate::{I256, ShrVartime};

    const N: I256 =
        I256::from_be_hex("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141");

    const N_2: I256 =
        I256::from_be_hex("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5D576E7357A4501DDFE92F46681B20A0");

    #[test]
    fn shr0() {
        assert_eq!(I256::MAX >> 0, I256::MAX);
        assert_eq!(I256::MIN >> 0, I256::MIN);
    }

    #[test]
    fn shr1() {
        assert_eq!(N >> 1, N_2);
        assert_eq!(ShrVartime::overflowing_shr_vartime(&N, 1), Some(N_2));
        assert_eq!(ShrVartime::wrapping_shr_vartime(&N, 1), N_2);
    }

    #[test]
    fn shr5() {
        assert_eq!(
            I256::MAX >> 5,
            I256::MAX.div(I256::from(32).to_nz().unwrap()).unwrap()
        );
        assert_eq!(
            I256::MIN >> 5,
            I256::MIN.div(I256::from(32).to_nz().unwrap()).unwrap()
        );
    }

    #[test]
    fn shr7_vartime() {
        assert_eq!(
            I256::MAX.shr_vartime(7),
            I256::MAX.div(I256::from(128).to_nz().unwrap()).unwrap()
        );
        assert_eq!(
            I256::MIN.shr_vartime(7),
            I256::MIN.div(I256::from(128).to_nz().unwrap()).unwrap()
        );
    }

    #[test]
    fn shr256_const() {
        assert!(N.overflowing_shr(256).is_none().to_bool());
        assert!(ShrVartime::overflowing_shr_vartime(&N, 256).is_none());
    }

    #[test]
    #[should_panic(expected = "`shift` exceeds upper bound")]
    fn shr_bounds_panic() {
        let _ = N >> 256;
    }

    #[test]
    fn unbounded_shr_vartime_zero_shift() {
        assert_eq!(I256::MAX.unbounded_shr_vartime(0), I256::MAX);
        assert_eq!(I256::MIN.unbounded_shr_vartime(0), I256::MIN);
        assert_eq!(I256::ONE.unbounded_shr_vartime(0), I256::ONE);
        assert_eq!(I256::MINUS_ONE.unbounded_shr_vartime(0), I256::MINUS_ONE);
        assert_eq!(I256::ZERO.unbounded_shr_vartime(0), I256::ZERO);
    }

    #[test]
    fn overflowing_shr_vartime_zero_shift() {
        let values = [I256::MAX, I256::MIN, I256::ONE, I256::MINUS_ONE, I256::ZERO];
        for &val in &values {
            assert_eq!(val.overflowing_shr_vartime(0), Some(val));
        }
    }

    #[test]
    fn shr_vartime_zero_shift() {
        let values = [I256::MAX, I256::MIN, I256::ONE, I256::MINUS_ONE, I256::ZERO];
        for &val in &values {
            assert_eq!(val.shr_vartime(0), val);
        }
    }

    #[test]
    fn wrapping_shr_vartime_at_or_past_the_width_fills_with_the_sign() {
        let values = [I256::MAX, I256::MIN, I256::ONE, I256::MINUS_ONE, I256::ZERO];
        for &val in &values {
            // Shift by 0 is the identity.
            assert_eq!(val.wrapping_shr_vartime(0), val);
            // Any shift of the full width or more leaves only the sign (rc.9 behaviour).
            let sign_fill = if val.is_negative().to_bool() {
                I256::MINUS_ONE
            } else {
                I256::ZERO
            };
            for i in 1..4 {
                assert_eq!(val.wrapping_shr_vartime(i * I256::BITS), sign_fill);
            }
        }
    }

    #[test]
    fn unbounded_shr() {
        assert_eq!(I256::MAX.unbounded_shr(257), I256::ZERO);
        assert_eq!(I256::MIN.unbounded_shr(257), I256::MINUS_ONE);
        assert_eq!(
            ShrVartime::unbounded_shr_vartime(&I256::MAX, 257),
            I256::ZERO
        );
        assert_eq!(
            ShrVartime::unbounded_shr_vartime(&I256::MIN, 257),
            I256::MINUS_ONE
        );
    }

    #[test]
    fn wrapping_shr() {
        assert_eq!(I256::MAX.wrapping_shr(257), I256::ZERO);
        assert_eq!(I256::MIN.wrapping_shr(257), I256::MINUS_ONE);
        assert_eq!(
            ShrVartime::wrapping_shr_vartime(&I256::MAX, 257),
            I256::ZERO
        );
        assert_eq!(
            ShrVartime::wrapping_shr_vartime(&I256::MIN, 257),
            I256::MINUS_ONE
        );
    }
}
