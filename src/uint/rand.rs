//! Random number generator support

use super::{Uint, Word};
use crate::{CtLt, Limb, NonZero, Random, RandomBits, RandomBitsError, RandomMod};
use rand_core::{Rng, TryRng};

impl<const LIMBS: usize> Random for Uint<LIMBS> {
    fn try_random_from_rng<R: TryRng + ?Sized>(rng: &mut R) -> Result<Self, R::Error> {
        let mut limbs = [Limb::ZERO; LIMBS];

        for limb in &mut limbs {
            *limb = Limb::try_random_from_rng(rng)?;
        }

        Ok(limbs.into())
    }
}

/// Fill the given limbs slice with random bits.
///
/// NOTE: Assumes that the limbs in the given slice are zeroed!
///
/// When combined with a platform-independent "4-byte sequential" `rng`, this function is
/// platform-independent. We consider an RNG "`X`-byte sequential" whenever
/// `rng.fill_bytes(&mut bytes[..i]); rng.fill_bytes(&mut bytes[i..])` constructs the same `bytes`,
/// as long as `i` is a multiple of `X`.
/// Note that the `TryRng` trait does _not_ require this behaviour from `rng`.
///
/// This is the `v0.7.0-rc.9` algorithm, kept in this fork so that values derived from a seeded
/// RNG stay the same as on the deployed revision: it draws whole words (4 or 8 bytes) per limb,
/// whereas upstream 0.7.5 draws exactly `bit_length.div_ceil(8)` bytes.
#[allow(clippy::integer_division_remainder_used, reason = "public parameter")]
pub(crate) fn random_bits_core<R: TryRng + ?Sized>(
    rng: &mut R,
    zeroed_limbs: &mut [Limb],
    bit_length: u32,
) -> Result<(), R::Error> {
    if bit_length == 0 {
        return Ok(());
    }

    let buffer: Word = 0;
    let mut buffer = buffer.to_be_bytes();

    let nonzero_limbs = bit_length.div_ceil(Limb::BITS) as usize;
    let partial_limb = bit_length % Limb::BITS;
    let mask = Word::MAX >> ((Word::BITS - partial_limb) % Word::BITS);

    for limb in zeroed_limbs.iter_mut().take(nonzero_limbs - 1) {
        rng.try_fill_bytes(&mut buffer)?;
        *limb = Limb(Word::from_le_bytes(buffer));
    }

    // This algorithm should sample the same number of random bytes, regardless of the pointer width
    // of the target platform. To this end, special attention has to be paid to the case where
    // bit_length - 1 < 32 mod 64. Bit strings of that size can be represented using `2X+1` 32-bit
    // words or `X+1` 64-bit words. Note that 64*(X+1) - 32*(2X+1) = 32. Hence, if we sample full
    // words only, a 64-bit platform will sample 32 bits more than a 32-bit platform. We prevent
    // this by forcing both platforms to only sample 4 bytes for the last word in this case.
    let slice = if partial_limb > 0 && partial_limb <= 32 {
        // Note: we do not have to zeroize the second half of the buffer, as the mask will take
        // care of this in the end.
        &mut buffer[0..4]
    } else {
        buffer.as_mut_slice()
    };

    rng.try_fill_bytes(slice)?;
    zeroed_limbs[nonzero_limbs - 1] = Limb(Word::from_le_bytes(buffer) & mask);

    Ok(())
}

impl<const LIMBS: usize> RandomBits for Uint<LIMBS> {
    fn try_random_bits<R: TryRng + ?Sized>(
        rng: &mut R,
        bit_length: u32,
    ) -> Result<Self, RandomBitsError<R::Error>> {
        Self::try_random_bits_with_precision(rng, bit_length, Self::BITS)
    }

    fn try_random_bits_with_precision<R: TryRng + ?Sized>(
        rng: &mut R,
        bit_length: u32,
        bits_precision: u32,
    ) -> Result<Self, RandomBitsError<R::Error>> {
        if bits_precision != Self::BITS {
            return Err(RandomBitsError::BitsPrecisionMismatch {
                bits_precision,
                integer_bits: Self::BITS,
            });
        }
        if bit_length > Self::BITS {
            return Err(RandomBitsError::BitLengthTooLarge {
                bit_length,
                bits_precision,
            });
        }
        let mut x = Self::ZERO;
        random_bits_core(rng, x.as_mut_limbs(), bit_length).map_err(RandomBitsError::RandCore)?;
        Ok(x)
    }
}

impl<const LIMBS: usize> RandomMod for Uint<LIMBS> {
    fn random_mod_vartime<R: Rng + ?Sized>(rng: &mut R, modulus: &NonZero<Self>) -> Self {
        let mut x = Self::ZERO;
        let Ok(()) = random_mod_vartime_core(rng, &mut x, modulus, modulus.bits_vartime());
        x
    }

    fn try_random_mod_vartime<R: TryRng + ?Sized>(
        rng: &mut R,
        modulus: &NonZero<Self>,
    ) -> Result<Self, R::Error> {
        let mut x = Self::ZERO;
        random_mod_vartime_core(rng, &mut x, modulus, modulus.bits_vartime())?;
        Ok(x)
    }
}

cpubits::cpubits! {
    32 => {
        /// Draw one limb-sized word from `rng`.
        ///
        /// `v0.7.0-rc.9` chose between `next_u32` and `next_u64` by `target_pointer_width`. On every
        /// target that revision shipped on, the pointer width equalled the limb width, so the draw
        /// was one word per limb. Since 0.7.5 the limb width comes from `cpubits`, which gives
        /// wasm32 64-bit limbs, so this fork chooses by the limb width instead: the sampler keeps
        /// drawing one word per limb, and a 64-bit-limb build on wasm32 draws exactly as native does.
        #[inline(always)]
        fn next_word<R: TryRng + ?Sized>(rng: &mut R) -> Result<Word, R::Error> {
            rng.try_next_u32()
        }
    }
    64 => {
        /// Draw one limb-sized word from `rng`.
        ///
        /// `v0.7.0-rc.9` chose between `next_u32` and `next_u64` by `target_pointer_width`. On every
        /// target that revision shipped on, the pointer width equalled the limb width, so the draw
        /// was one word per limb. Since 0.7.5 the limb width comes from `cpubits`, which gives
        /// wasm32 64-bit limbs, so this fork chooses by the limb width instead: the sampler keeps
        /// drawing one word per limb, and a 64-bit-limb build on wasm32 draws exactly as native does.
        #[inline(always)]
        fn next_word<R: TryRng + ?Sized>(rng: &mut R) -> Result<Word, R::Error> {
            rng.try_next_u64()
        }
    }
}

/// Generic implementation of `random_mod_vartime` which can be shared with `BoxedUint`.
///
/// This is the `v0.7.0-rc.9` algorithm (then named `random_mod_core`), kept in this fork so that
/// values derived from a seeded RNG stay the same as on the deployed revision: it draws one word
/// per limb, rejects on the top word first and then on the whole value. Upstream 0.7.5 instead
/// draws `n_bits` bits through `random_bits_core` and rejects on the whole value only. The one
/// departure from the rc.9 text is how the word size is chosen; see [`next_word`].
// TODO(tarcieri): obtain `n_bits` via a trait like `Integer`
pub(super) fn random_mod_vartime_core<T, R: TryRng + ?Sized>(
    rng: &mut R,
    n: &mut T,
    modulus: &NonZero<T>,
    n_bits: u32,
) -> Result<(), R::Error>
where
    T: AsMut<[Limb]> + AsRef<[Limb]> + CtLt,
{
    let n_limbs = n_bits.div_ceil(Limb::BITS) as usize;

    let hi_word_modulus = modulus.as_ref().as_ref()[n_limbs - 1].0;
    let mask = !0 >> hi_word_modulus.leading_zeros();
    let mut hi_word = next_word(rng)? & mask;

    loop {
        while hi_word > hi_word_modulus {
            hi_word = next_word(rng)? & mask;
        }
        // Set high limb
        n.as_mut()[n_limbs - 1] = Limb::from_le_bytes(hi_word.to_le_bytes());
        // Set low limbs
        for i in 0..n_limbs - 1 {
            // Need to deserialize from little-endian to make sure that two 32-bit limbs
            // deserialized sequentially are equal to one 64-bit limb produced from the same
            // byte stream.
            n.as_mut()[i] = Limb::from_le_bytes(next_word(rng)?.to_le_bytes());
        }
        // If the high limb is equal to the modulus' high limb, it's still possible
        // that the full uint is too big so we check and repeat if it is.
        if n.ct_lt(modulus).into() {
            break;
        }
        hi_word = next_word(rng)? & mask;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::uint::rand::random_bits_core;
    use crate::{Limb, NonZero, Random, RandomBits, RandomMod, U256, U1024, Uint};
    use chacha20::ChaCha8Rng;
    use rand_core::{Rng, SeedableRng};

    const RANDOM_OUTPUT: U1024 = Uint::from_be_hex(concat![
        "A484C4C693EECC47C3B919AE0D16DF2259CD1A8A9B8EA8E0862878227D4B40A3",
        "C54302F2EB1E2F69E17653A37F1BCC44277FA208E6B31E08CDC4A23A7E88E660",
        "EF781C7DD2D368BAD438539D6A2E923C8CAE14CB947EB0BDE10D666732024679",
        "0F6760A48F9B887CB2FB0D3281E2A6E67746A55FBAD8C037B585F767A79A3B6C"
    ]);

    /// Construct a 4-sequential `rng`, i.e., an `rng` such that
    /// `rng.fill_bytes(&mut buffer[..x]); rng.fill_bytes(&mut buffer[x..])` will construct the
    /// same `buffer`, for `x` any in `0..buffer.len()` that is `0 mod 4`.
    fn get_four_sequential_rng() -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(0)
    }

    /// Make sure the random value constructed is consistent across platforms
    #[test]
    fn random_platform_independence() {
        let mut rng = get_four_sequential_rng();
        assert_eq!(U1024::random_from_rng(&mut rng), RANDOM_OUTPUT);
    }

    #[test]
    fn random_mod_vartime() {
        let mut rng = ChaCha8Rng::seed_from_u64(1);

        // Ensure `random_mod_vartime` runs in a reasonable amount of time
        let modulus = NonZero::new(U256::from(42u8)).unwrap();
        let res = U256::random_mod_vartime(&mut rng, &modulus);

        // Check that the value is in range
        assert!(res < U256::from(42u8));

        // Ensure `random_mod_vartime` runs in a reasonable amount of time
        // when the modulus is larger than 1 limb
        let modulus = NonZero::new(U256::from(0x10000000000000001u128)).unwrap();
        let res = U256::random_mod_vartime(&mut rng, &modulus);

        // Check that the value is in range
        assert!(res < U256::from(0x10000000000000001u128));
    }

    #[test]
    fn random_bits() {
        let mut rng = ChaCha8Rng::seed_from_u64(1);

        let lower_bound = 16;

        // Full length of the integer
        let bit_length = U256::BITS;
        for _ in 0..10 {
            let res = U256::random_bits(&mut rng, bit_length);
            assert!(res > (U256::ONE << (bit_length - lower_bound)));
        }

        // A multiple of limb size
        let bit_length = U256::BITS - Limb::BITS;
        for _ in 0..10 {
            let res = U256::random_bits(&mut rng, bit_length);
            assert!(res > (U256::ONE << (bit_length - lower_bound)));
            assert!(res < (U256::ONE << bit_length));
        }

        // A multiple of 8
        let bit_length = U256::BITS - Limb::BITS - 8;
        for _ in 0..10 {
            let res = U256::random_bits(&mut rng, bit_length);
            assert!(res > (U256::ONE << (bit_length - lower_bound)));
            assert!(res < (U256::ONE << bit_length));
        }

        // Not a multiple of 8
        let bit_length = U256::BITS - Limb::BITS - 8 - 3;
        for _ in 0..10 {
            let res = U256::random_bits(&mut rng, bit_length);
            assert!(res > (U256::ONE << (bit_length - lower_bound)));
            assert!(res < (U256::ONE << bit_length));
        }

        // One incomplete limb
        let bit_length = 7;
        for _ in 0..10 {
            let res = U256::random_bits(&mut rng, bit_length);
            assert!(res < (U256::ONE << bit_length));
        }

        // Zero bits
        let bit_length = 0;
        for _ in 0..10 {
            let res = U256::random_bits(&mut rng, bit_length);
            assert_eq!(res, U256::ZERO);
        }
    }

    /// Make sure the `random_bits` output is consistent across platforms
    #[test]
    fn random_bits_platform_independence() {
        let mut rng = get_four_sequential_rng();

        let bit_length = 989;
        let mut val = U1024::ZERO;
        random_bits_core(&mut rng, val.as_mut_limbs(), bit_length).expect("safe");

        assert_eq!(
            val,
            RANDOM_OUTPUT.bitand(&U1024::ONE.shl(bit_length).wrapping_sub(&Uint::ONE))
        );

        // Test that the RNG is in the same state
        let mut state = [0u8; 16];
        rng.fill_bytes(&mut state);

        assert_eq!(
            state,
            [
                198, 196, 132, 164, 240, 211, 223, 12, 36, 189, 139, 48, 94, 1, 123, 253
            ]
        );
    }

    /// The `random_bits` draw rounds the last partial limb up to 4 or 8 bytes, as `v0.7.0-rc.9`
    /// did; upstream 0.7.5 draws exactly `bit_length.div_ceil(8)` bytes. 65 bits tells them
    /// apart: rc.9 leaves the RNG 12 bytes further on, 0.7.5 only 9.
    #[test]
    fn random_bits_draws_whole_words_as_rc9_did() {
        let mut rng = get_four_sequential_rng();
        let mut val = U1024::ZERO;
        random_bits_core(&mut rng, val.as_mut_limbs(), 65).expect("safe");

        let mut reference = get_four_sequential_rng();
        let mut drawn = [0u8; 12];
        reference.fill_bytes(&mut drawn);

        let mut state = [0u8; 16];
        rng.fill_bytes(&mut state);
        let mut reference_state = [0u8; 16];
        reference.fill_bytes(&mut reference_state);
        assert_eq!(state, reference_state);
    }

    /// Printed by `U256::random_mod` of crypto-bigint `v0.7.0-rc.9` (ycscaly/crypto-bigint
    /// `8aabcee5`, rand_chacha 0.9) on aarch64, for the sequence below.
    const RC9_SMALL_MODULUS_OUTPUTS: [u32; 5] = [55, 2172, 1657, 4668, 7688];
    const RC9_LARGE_MODULUS_OUTPUT: &str =
        "C54302F2EB1E2F69C3B919AE0D16DF2259CD1A8A9B8EA8E0862878227D4B40A3";
    const RC9_FINAL_RNG_STATE: [u8; 16] = [
        71, 204, 238, 147, 198, 196, 132, 164, 240, 211, 223, 12, 36, 189, 139, 48,
    ];

    cpubits::cpubits! {
        64 => {
            /// Fixed outputs of the `v0.7.0-rc.9` modulus sampler on 64-bit limbs, obtained by running
            /// this same sequence on the `v0.7.0-rc.9` crate itself. The rc.9 algorithm draws one word
            /// per limb, so its outputs depend on the limb width; this test is therefore limited to
            /// 64-bit limbs. (Upstream 0.7.5 had a `random_mod_vartime_platform_independence` test
            /// here with the outputs of its byte-based sampler, which this fork does not use.)
            #[test]
            fn random_mod_vartime_matches_rc9_on_64_bit_limbs() {
                let mut rng = get_four_sequential_rng();

                let modulus = NonZero::new(U256::from_u32(8192)).unwrap();
                let mut vals = [U256::ZERO; 5];
                for val in &mut vals {
                    *val = U256::random_mod_vartime(&mut rng, &modulus);
                }
                let expected = RC9_SMALL_MODULUS_OUTPUTS;
                for (want, got) in expected.into_iter().zip(vals.into_iter()) {
                    assert_eq!(got, U256::from_u32(want));
                }

                let modulus =
                    NonZero::new(U256::ZERO.wrapping_sub(&U256::from_u64(rng.next_u64()))).unwrap();
                let val = U256::random_mod_vartime(&mut rng, &modulus);
                assert_eq!(val, U256::from_be_hex(RC9_LARGE_MODULUS_OUTPUT));

                let mut state = [0u8; 16];
                rng.fill_bytes(&mut state);

                assert_eq!(state, RC9_FINAL_RNG_STATE);
            }
        }
    }
}
