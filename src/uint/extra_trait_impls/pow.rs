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

                    let mut limbs = [Limb::ZERO; nlimbs!($first_bits)];
                    let mut j = 0;

                    while j < nlimbs!($second_bits % $first_bits) {
                        limbs[j] = exponent.limbs[i*nlimbs!($first_bits) + j];
                        j += 1;
                    }

                    let exponent_i: $first_type = Uint { limbs };

                    res *= shifted_base.pow_bounded_exp(
                        &exponent_i,
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
    use crate::{NonZero, Pow, U128, U192, U256, U320, U384, U4096, U4352};

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

    #[test]
    fn pow_u4096_by_u4352() {
        let n = U4096::from_be_hex("5960383b5378ad0607f0f270ce7fb6dcaba6506f9fc56deeffaf605c9128db8ccf063e2e8221a8bdf82c027741a0303b08eb71fa6225a03df18f24c473dc6d4d3d30eb9c52a233bbfe967d04011b95e8de5bc482c3c217bcfdeb4df6f57af6ba9c6d66c69fb03a70a41fe1e87975c85343ef7d572ca06a0139706b23ed2b73ad72cb1b7e2e41840115651897c8757b3da9af3a60eebb6396ffd193738b4f04aa6ece638cef1bf4e9c45cf57f8debeda8598cbef732484752f5380737ba75ee00bf1b146817b9ab336d0ce5540395377347c653d1c9d272127ff12b9a0721b8ef13ecd8a8379f1b9a358de2af2c4cd97564dbd5328c2fc13d56ee30c8a101d333f5406afb1f4417b49d7a629d5076726877df11f05c998ae365e374a0141f0b99802214532c97c1ebf9faf6e277a8f29dbd8f3eab72266e60a77784249694819e42877a5e826745c97f84a5f37002b74d83fc064cf094be0e706a6710d47d253c4532e6aa4a679a75fa1d860b39085dab03186c67248e6c92223682f58bd41b67143e299329ce3a8045f3a0124c3d0ef9f0f49374d89b37d9c3321feb2ab4117df4f68246724ce41cd765326457968d848afcc0735531e5de7fea88cf2eb35ac68710c6e79d5ad25df6c0393c0267f56e8eac90a52637abe3e606769e70b20560eaf70e0d531b11dca299104fa933f887d85fb5f72386c196e40f559baee356b9");
        let r = U4352::from_be_hex("7303C065D35BF556E4A3DB576081924C45846B3E83CE1A9178EF8D5693ECF15F670D6CDF7EDEFC390FCBEF6090021BDD3A3EEA6C47A03912AB2F85E22804DA21B523ACE4A69E2C3657AAB5343A44C2A46BD2426894B8923C9905F3B7D671DB33F9DE5AC01E525064BFD4E39D230D5B97832348D15EF5725A59A2CA5795085CB078559F822FDED48EC76E2D5B4601E1E48BF8E6BD24FBC79F9131D4D581CC2EE430E4B77F38BBBC0851C12DB448F3E9E7DC7B5ECC6CFA69110990F906E37393E3C3CE4CCDC3459071EE076F14FA123FF25B727279677E433719806F8CFB5D5EDF50BABE66791A38333F3C6AFF72866220A602C6AE0BD5CE9A22A67E46F876B977C8CB4A43B7AD16990928AA0B6145C9054C2B8A251DCF4225CBB7FCC6E6AB821338013C8F7B30EB584FA890D1C627EAA3A516E6FEF425FDBFD7061AB095BB71ED3EB88751235F96E227202B952AC1974602D626673CF261085FC6624785228A794972D8D0EA15B2B9C9C844C1C7AE3069140A858453D789F9085074D72F6D912D48DAAF976FF62C75137C653E446B2514A189E0D08DB56B15879560E77D7BD91C80327B438E2D5CE6B305F0539A51F1006A257C3565CBA76DF85ED3B361502456F59E3704BCCB9D8A6A394FA646889080A4D92F45A44D7554970C7EE21EA7762BCEA6F855794158880E54B38F35FF1201C0846FC5448D61D3166E9C0C230E49A6B3F3C825FE599E4DBDDFB4697BA5467ADBE2B50FBD0BA58A777D2528C861E430");
        let g = U4096::from_be_hex("19BB1B2E0015AA04BEE4F8321819448A2C809DF799C6627668DAA936E3A367CF87BEC43C47551221E40724FE115FF8A4E72D5D46A0E98A934C45CD6904DA0F07499D798EE611497C9493354A9A48C35ECB6318CA55B8322E4295E67F8BC0BE1E0923685E1727B7925920D4F0E9CC30C2A10135DB447EDAD3BCE87C3416252C8B4DF32C24029E0269E7103E80D02DD5A42A99B69A613C6274255DF0599B0DED35A8969463636C6D56D67A05AE11F347A5D5B81896DF5F8A52E6EA7F05359A9FEFC90297BDD298DD77714D3557325DF1C52F42470606ECBFA5E964C0A782AE19CED2E20C73F0438EB597CAE4159B5E5333C97272D8EFEDB49CEB98078E92D990076E6E4101FD97588E4BBAA9DD5D19C671424108EE7FA5F2D74F9F3DEAB4A0AC89CF9833FD9BA1F66719978D7BD13DD2ECDE2BDC9628B1AC1E0A0C44B1408E8869A8B2245DF2A877E01730500AD15466A808E6D9636EEA7A7A0A06568413408E588C52451D189774D84547FBB4171255D6E0BFC9B63C56D582E02FA0F110EEAA2B728E51BC85F529805EBA5E1D6B7323597F1647B0A3DC6D61448C1C062CADE9831DB9E3029322D79D04BB3287B7C5D857AE11802B68921FBC403E390ED693DEAD66E1A728B7F7432408EB2ED9EB9BC3B2BCD8EB2CD44D41A5EBFB32F55BAF47D3AC048F5D1F60B2CB61C0F4E3C178DC7723B8298E9D52771DCF1DABA4088EF74B");
        let g = g % NonZero::new(n).unwrap();
        let params = DynResidueParams::new(&n);
        let g = DynResidue::new(&g, params);

        assert_eq!(
            <DynResidue::<{ U4096::LIMBS }> as Pow<U4352>>::pow(&g, &r).retrieve(),
            U4096::from_be_hex("3BFC8BD71A8382135C7D8D00E5043038AD76DB4F39D94595F4D7F1A098465ADA0950BE29DB1DA00867B75F88776DF48F59B6E838F4EE10A8B42D654BFE330E91A5BE452FE9B26E8E9F720487D6E1C07BD06A7B00D7260201686EE5D877380F9F74316A609A11E5037BD82C1FC32D9F9FAEF24CC883807815BD37DEBA61A96FAE475FC5E9B1DF16471D6A5B22C57A2FB4D29270AA4566DA6A25CE1679DECE35FD2ADF27DB37CC3B1349BCBFB93D6F7E1BCC78FE056E9E1F8A23BEEF913A32EEB9EE9CDEF7595BEC881B3042842202EA6EA8F4313C457C1887B63A18E1F8D2A6D2AFFFF453E7047744323B6B27FDBB0F26840FD752395598AF28FC27FB210C93C617833C28BA4E1F74D87EF5A6BC3204E04254B95D5274FCD61A028B60B6D5D9ADD1C72B2540BFDA22BB204E0C6A44499DAA9211DB3AE32307A51F8D4D8C836A231FAADBA207043B95879CC097F2600E17F71BA8B452A91C4A5DBFA2B9B6ACB3CD4F80F2C3F0906D70ED589C25B679B04E4F255566BE8AC59101F8A8827D2BE07C094E4FC4FD11521680F7AC603E8D04C56879CA146C92CB7E24B60842719BACC931346F97392CFF97A9D1CD5EC2DF23288FDFB33445CB0CD72CAFFDF7F27C1F443DD16CD374C9C92060C1D022585264D638EFBB1C209C62FCB6FA9C29856C58E4EF1E0E20461565DFA7B15120DF5F3997435B364E2D874C4B3B8C668FEE81C2EB")
        )
    }
}
