//! Random number generator support

use super::{Uint, Word};
use crate::{Limb, NonZero, Random, RandomBits, RandomBitsError, RandomMod, Zero};
use rand_core::{RngCore, TryRngCore};
use subtle::ConstantTimeLess;

impl<const LIMBS: usize> Random for Uint<LIMBS> {
    fn try_random<R: TryRngCore + ?Sized>(rng: &mut R) -> Result<Self, R::Error> {
        let mut limbs = [Limb::ZERO; LIMBS];

        for limb in &mut limbs {
            *limb = Limb::try_random(rng)?
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
/// Note that the `TryRngCore` trait does _not_ require this behaviour from `rng`.
pub(crate) fn random_bits_core<R: TryRngCore + ?Sized>(
    rng: &mut R,
    zeroed_limbs: &mut [Limb],
    bit_length: u32,
) -> Result<(), RandomBitsError<R::Error>> {
    if bit_length == 0 {
        return Ok(());
    }

    let buffer: Word = 0;
    let mut buffer = buffer.to_be_bytes();

    let nonzero_limbs = bit_length.div_ceil(Limb::BITS) as usize;
    let partial_limb = bit_length % Limb::BITS;
    let mask = Word::MAX >> ((Word::BITS - partial_limb) % Word::BITS);

    for i in 0..nonzero_limbs - 1 {
        rng.try_fill_bytes(&mut buffer)
            .map_err(RandomBitsError::RandCore)?;
        zeroed_limbs[i] = Limb(Word::from_le_bytes(buffer));
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

    rng.try_fill_bytes(slice)
        .map_err(RandomBitsError::RandCore)?;
    zeroed_limbs[nonzero_limbs - 1] = Limb(Word::from_le_bytes(buffer) & mask);

    Ok(())
}

impl<const LIMBS: usize> RandomBits for Uint<LIMBS> {
    fn try_random_bits<R: TryRngCore + ?Sized>(
        rng: &mut R,
        bit_length: u32,
    ) -> Result<Self, RandomBitsError<R::Error>> {
        Self::try_random_bits_with_precision(rng, bit_length, Self::BITS)
    }

    fn try_random_bits_with_precision<R: TryRngCore + ?Sized>(
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
        let mut limbs = [Limb::ZERO; LIMBS];
        random_bits_core(rng, &mut limbs, bit_length)?;
        Ok(Self::from(limbs))
    }
}

impl<const LIMBS: usize> RandomMod for Uint<LIMBS> {
    fn random_mod<R: RngCore + ?Sized>(rng: &mut R, modulus: &NonZero<Self>) -> Self {
        let mut n = Self::ZERO;
        let Ok(()) = random_mod_core(rng, &mut n, modulus, modulus.bits_vartime());
        n
    }

    fn try_random_mod<R: TryRngCore + ?Sized>(
        rng: &mut R,
        modulus: &NonZero<Self>,
    ) -> Result<Self, R::Error> {
        let mut n = Self::ZERO;
        random_mod_core(rng, &mut n, modulus, modulus.bits_vartime())?;
        Ok(n)
    }
}

/// The width of the words `random_mod` draws from the RNG and places into the integer.
///
/// The deployed sampler drew one machine word per limb, so it drew 64-bit words on native
/// targets and 32-bit words on `wasm32`, and the same seed gave a different value on each.
/// `random_mod` now always draws 64-bit words, on every target, so the value is the same
/// everywhere and equals what the deployed native build drew. `ThirtyTwo` reproduces the
/// deployed 32-bit draw, for values sampled by the deployed `wasm32` build.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RandomModWordWidth {
    /// 32-bit words: the deployed draw on `wasm32`.
    ThirtyTwo,
    /// 64-bit words: the deployed draw on native targets, and the draw everywhere now.
    SixtyFour,
}

impl RandomModWordWidth {
    const fn bits(self) -> u32 {
        match self {
            Self::ThirtyTwo => 32,
            Self::SixtyFour => 64,
        }
    }
}

/// Generic implementation of `random_mod` which can be shared with `BoxedUint`.
///
/// Draws 64-bit words whatever the limb width is, so that every target computes the value the
/// deployed native build computed; see [`random_mod_core_with_word_width`].
// TODO(tarcieri): obtain `n_bits` via a trait like `Integer`
pub(super) fn random_mod_core<T, R: TryRngCore + ?Sized>(
    rng: &mut R,
    n: &mut T,
    modulus: &NonZero<T>,
    n_bits: u32,
) -> Result<(), R::Error>
where
    T: AsMut<[Limb]> + AsRef<[Limb]> + ConstantTimeLess + Zero,
{
    random_mod_core_with_word_width(rng, n, modulus, n_bits, RandomModWordWidth::SixtyFour)
}

/// Rejection sampling below `modulus`, in words of `width` bits regardless of the limb width.
///
/// The algorithm is the deployed one: the top word of the modulus (the word holding bit
/// `n_bits - 1`) fixes a mask; a top word is drawn and masked until it does not exceed the
/// modulus's top word; the lower words are drawn one by one; if the whole value is not below
/// the modulus, a new top word is drawn and everything repeats. The deployed code drew one
/// machine word per limb. Here the words are `width` bits wide and are placed into the limbs
/// little-endian, so with `width` equal to the limb width this consumes the RNG exactly as the
/// deployed code did, and with the other width it reproduces the other target's draw.
///
/// `n` must be zero on entry.
pub(super) fn random_mod_core_with_word_width<T, R: TryRngCore + ?Sized>(
    rng: &mut R,
    n: &mut T,
    modulus: &NonZero<T>,
    n_bits: u32,
    width: RandomModWordWidth,
) -> Result<(), R::Error>
where
    T: AsMut<[Limb]> + AsRef<[Limb]> + ConstantTimeLess + Zero,
{
    let word_bits = width.bits();
    let mut next_word = || -> Result<u64, R::Error> {
        match width {
            RandomModWordWidth::SixtyFour => rng.try_next_u64(),
            RandomModWordWidth::ThirtyTwo => rng.try_next_u32().map(u64::from),
        }
    };

    let n_words = n_bits.div_ceil(word_bits) as usize;

    let hi_word_modulus = word_at(modulus.as_ref().as_ref(), n_words - 1, word_bits);
    // The mask of the deployed code, `!0 >> hi_word_modulus.leading_zeros()`, computed in a
    // word of `word_bits` bits.
    let mask =
        (u64::MAX >> (64 - word_bits)) >> (hi_word_modulus.leading_zeros() - (64 - word_bits));
    let mut hi_word = next_word()? & mask;

    loop {
        while hi_word > hi_word_modulus {
            hi_word = next_word()? & mask;
        }
        // Set high word
        set_word(n.as_mut(), n_words - 1, hi_word, word_bits);
        // Set low words
        for i in 0..n_words - 1 {
            let word = next_word()?;
            set_word(n.as_mut(), i, word, word_bits);
        }
        // If the high word is equal to the modulus' high word, it's still possible
        // that the full uint is too big so we check and repeat if it is.
        if n.ct_lt(modulus).into() {
            break;
        }
        hi_word = next_word()? & mask;
    }
    Ok(())
}

/// The `index`-th little-endian word of `word_bits` bits of `limbs`. A word beyond the limbs
/// is zero.
// The casts between `Word` and `u64` are identities on 64-bit targets and conversions on
// 32-bit targets; both targets compile every branch.
#[allow(trivial_numeric_casts, clippy::unnecessary_cast, clippy::useless_conversion)]
fn word_at(limbs: &[Limb], index: usize, word_bits: u32) -> u64 {
    if word_bits == Limb::BITS {
        u64::from(limbs[index].0)
    } else if word_bits > Limb::BITS {
        // 64-bit words over 32-bit limbs: two limbs per word, low limb first.
        let low = u64::from(limbs[2 * index].0);
        let high = limbs.get(2 * index + 1).map_or(0, |limb| u64::from(limb.0));
        low | (high << 32)
    } else {
        // 32-bit words over 64-bit limbs: two words per limb, low word first.
        let shift = 32 * (index % 2) as u32;
        (u64::from(limbs[index / 2].0) >> shift) & u64::from(u32::MAX)
    }
}

/// Store `word` as the `index`-th little-endian word of `word_bits` bits of `limbs`.
// The casts between `Word` and `u64` are identities on 64-bit targets and conversions on
// 32-bit targets; both targets compile every branch.
#[allow(trivial_numeric_casts, clippy::unnecessary_cast, clippy::useless_conversion)]
fn set_word(limbs: &mut [Limb], index: usize, word: u64, word_bits: u32) {
    if word_bits == Limb::BITS {
        limbs[index] = Limb(word as Word);
    } else if word_bits > Limb::BITS {
        limbs[2 * index] = Limb(word as u32 as Word);
        match limbs.get_mut(2 * index + 1) {
            Some(high) => *high = Limb((word >> 32) as u32 as Word),
            // The top word of an odd number of 32-bit limbs: its high half is masked to zero.
            None => debug_assert_eq!(word >> 32, 0),
        }
    } else {
        let shift = 32 * (index % 2) as u32;
        let current = u64::from(limbs[index / 2].0);
        let cleared = current & !(u64::from(u32::MAX) << shift);
        limbs[index / 2] = Limb((cleared | (word << shift)) as Word);
    }
}

impl<const LIMBS: usize> Uint<LIMBS> {
    /// Sample uniformly below `modulus`, drawing words of `width` bits from `rng`.
    ///
    /// `random_mod` draws 64-bit words. This reproduces the deployed draw of the other word
    /// width as well: with [`RandomModWordWidth::ThirtyTwo`] it computes the value the deployed
    /// `wasm32` build computed from the same RNG stream, on any target.
    pub fn random_mod_with_word_width<R: RngCore + ?Sized>(
        rng: &mut R,
        modulus: &NonZero<Self>,
        width: RandomModWordWidth,
    ) -> Self {
        let mut n = Self::ZERO;
        let Ok(()) =
            random_mod_core_with_word_width(rng, &mut n, modulus, modulus.bits_vartime(), width);
        n
    }

    /// Fallible form of [`Self::random_mod_with_word_width`].
    pub fn try_random_mod_with_word_width<R: TryRngCore + ?Sized>(
        rng: &mut R,
        modulus: &NonZero<Self>,
        width: RandomModWordWidth,
    ) -> Result<Self, R::Error> {
        let mut n = Self::ZERO;
        random_mod_core_with_word_width(rng, &mut n, modulus, modulus.bits_vartime(), width)?;
        Ok(n)
    }
}

/// The deployed `random_mod_core` (revision 8aabcee5), kept verbatim as the reference the
/// tests compare the word-width-independent implementation against. It draws one machine word
/// per limb.
#[cfg(test)]
pub(super) fn random_mod_core_deployed<T, R: TryRngCore + ?Sized>(
    rng: &mut R,
    n: &mut T,
    modulus: &NonZero<T>,
    n_bits: u32,
) -> Result<(), R::Error>
where
    T: AsMut<[Limb]> + AsRef<[Limb]> + ConstantTimeLess + Zero,
{
    use crate::Encoding;

    #[cfg(target_pointer_width = "64")]
    let mut next_word = || rng.try_next_u64();
    #[cfg(target_pointer_width = "32")]
    let mut next_word = || rng.try_next_u32();

    let n_limbs = n_bits.div_ceil(Limb::BITS) as usize;

    let hi_word_modulus = modulus.as_ref().as_ref()[n_limbs - 1].0;
    let mask = !0 >> hi_word_modulus.leading_zeros();
    let mut hi_word = next_word()? & mask;

    loop {
        while hi_word > hi_word_modulus {
            hi_word = next_word()? & mask;
        }
        n.as_mut()[n_limbs - 1] = Limb::from_le_bytes(hi_word.to_le_bytes());
        for i in 0..n_limbs - 1 {
            n.as_mut()[i] = Limb::from_le_bytes(next_word()?.to_le_bytes());
        }
        if n.ct_lt(modulus).into() {
            break;
        }
        hi_word = next_word()? & mask;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::uint::rand::{random_bits_core, random_mod_core_deployed};
    use crate::{
        Limb, NonZero, Random, RandomBits, RandomMod, RandomModWordWidth, U64, U128, U256, U512,
        U1024, U2048, U4096, Uint,
    };
    use rand_chacha::ChaCha8Rng;
    use rand_core::{RngCore, SeedableRng};

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
        assert_eq!(U1024::random(&mut rng), RANDOM_OUTPUT);
    }

    #[test]
    fn random_mod() {
        let mut rng = ChaCha8Rng::seed_from_u64(1);

        // Ensure `random_mod` runs in a reasonable amount of time
        let modulus = NonZero::new(U256::from(42u8)).unwrap();
        let res = U256::random_mod(&mut rng, &modulus);

        // Check that the value is in range
        assert!(res < U256::from(42u8));

        // Ensure `random_mod` runs in a reasonable amount of time
        // when the modulus is larger than 1 limb
        let modulus = NonZero::new(U256::from(0x10000000000000001u128)).unwrap();
        let res = U256::random_mod(&mut rng, &modulus);

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

    /// Make sure the random_bits output is consistent across platforms
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

    /// Test that random bytes are sampled consecutively.
    #[test]
    fn random_bits_4_bytes_sequential() {
        // Test for multiples of 4 bytes, i.e., multiples of 32 bits.
        let bit_lengths = [0, 32, 64, 128, 192, 992];

        for bit_length in bit_lengths {
            let mut rng = get_four_sequential_rng();
            let mut first = U1024::ZERO;
            let mut second = U1024::ZERO;
            random_bits_core(&mut rng, first.as_mut_limbs(), bit_length).expect("safe");
            random_bits_core(&mut rng, second.as_mut_limbs(), U1024::BITS - bit_length)
                .expect("safe");
            assert_eq!(second.shl(bit_length).bitor(&first), RANDOM_OUTPUT);
        }
    }

    /// The word width that reproduces the deployed draw on this target.
    const DEPLOYED_WIDTH: RandomModWordWidth = if Limb::BITS == 64 {
        RandomModWordWidth::SixtyFour
    } else {
        RandomModWordWidth::ThirtyTwo
    };

    /// Draws `count` values below `modulus` with both samplers from identical RNG streams,
    /// asserting equal values and an equal RNG state after each draw.
    fn assert_matches_deployed<const LIMBS: usize>(modulus: Uint<LIMBS>, seed: u64, count: usize) {
        let modulus = NonZero::new(modulus).unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mut deployed_rng = ChaCha8Rng::seed_from_u64(seed);
        for _ in 0..count {
            let value =
                Uint::<LIMBS>::random_mod_with_word_width(&mut rng, &modulus, DEPLOYED_WIDTH);
            let mut deployed = Uint::<LIMBS>::ZERO;
            let Ok(()) = random_mod_core_deployed(
                &mut deployed_rng,
                &mut deployed,
                &modulus,
                modulus.bits_vartime(),
            );
            assert_eq!(value, deployed);
            assert!(value < *modulus);
            assert_eq!(rng.next_u64(), deployed_rng.next_u64());
        }
    }

    /// Checks every fixed modulus below and one random modulus of every bit length, so that
    /// every top-word mask width occurs, `count` draws each.
    fn check_every_modulus<const LIMBS: usize>(seed: u64, count: usize) {
        let fixed = [
            Uint::<LIMBS>::ONE,
            Uint::<LIMBS>::from(42u8),
            Uint::<LIMBS>::from(u32::MAX),
            Uint::<LIMBS>::from(u64::MAX),
            Uint::<LIMBS>::ONE.shl_vartime(32),
            Uint::<LIMBS>::ONE.shl_vartime(33).wrapping_add(&Uint::ONE),
            Uint::<LIMBS>::ONE.shl_vartime(Uint::<LIMBS>::BITS - 1),
            Uint::<LIMBS>::MAX,
        ];
        for modulus in fixed {
            assert_matches_deployed(modulus, seed, count);
        }
        let mut rng = ChaCha8Rng::seed_from_u64(seed ^ 0xA5A5);
        for bits in 1..=Uint::<LIMBS>::BITS {
            let modulus = Uint::<LIMBS>::random_bits(&mut rng, bits)
                | Uint::<LIMBS>::ONE.shl_vartime(bits - 1);
            assert_matches_deployed(modulus, seed, count);
        }
    }

    /// On this target, drawing with the deployed word width through the width-independent
    /// implementation gives the deployed values and leaves the RNG in the deployed state, for
    /// every modulus bit length at every integer size the workspace uses.
    #[test]
    fn matches_deployed_sampler_on_this_target() {
        check_every_modulus::<{ U64::LIMBS }>(1, 8);
        check_every_modulus::<{ U128::LIMBS }>(2, 8);
        check_every_modulus::<{ U256::LIMBS }>(3, 8);
        check_every_modulus::<{ U512::LIMBS }>(4, 4);
        check_every_modulus::<{ U1024::LIMBS }>(5, 2);
        check_every_modulus::<{ U2048::LIMBS }>(6, 2);
        check_every_modulus::<{ U4096::LIMBS }>(7, 1);
    }

    /// `random_mod` is the 64-bit draw on every target.
    #[test]
    fn random_mod_is_the_sixty_four_bit_draw() {
        let modulus = NonZero::new(U256::MAX.wrapping_sub(&U256::from(12345u32))).unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(9);
        let mut other = ChaCha8Rng::seed_from_u64(9);
        for _ in 0..64 {
            assert_eq!(
                U256::random_mod(&mut rng, &modulus),
                U256::random_mod_with_word_width(
                    &mut other,
                    &modulus,
                    RandomModWordWidth::SixtyFour
                )
            );
        }
    }

    /// Values fixed by the deployed native build (64-bit words) and the deployed wasm32 build
    /// (32-bit words) from the same seed and modulus, so that each width's draw is pinned on
    /// every target. Captured with revision 8aabcee5 on each target.
    const SIXTY_FOUR_BIT_DRAW: &str =
        "0000000000000000000000000000000000000000000000007746A55FBAD8C037";
    const THIRTY_TWO_BIT_DRAW: &str =
        "000000000000000000000000000000000000000000000000BAD8C037B585F767";

    #[test]
    fn both_widths_are_pinned() {
        let modulus = NonZero::new(U256::ONE.shl_vartime(65).wrapping_add(&U256::ONE)).unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let sixty_four =
            U256::random_mod_with_word_width(&mut rng, &modulus, RandomModWordWidth::SixtyFour);
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let thirty_two =
            U256::random_mod_with_word_width(&mut rng, &modulus, RandomModWordWidth::ThirtyTwo);
        assert_ne!(sixty_four, thirty_two);
        assert_eq!(sixty_four, U256::from_be_hex(SIXTY_FOUR_BIT_DRAW));
        assert_eq!(thirty_two, U256::from_be_hex(THIRTY_TWO_BIT_DRAW));
    }
}
