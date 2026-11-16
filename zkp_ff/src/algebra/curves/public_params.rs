/// 椭圆曲线公共参数接口。
///
/// 每种曲线都需要实现这个 trait，
/// 例如 BN254、BLS12-381。
pub trait PublicParams: Sized {
    // ============================================================
    // Associated Types
    // ============================================================

    type Fr;
    type G1;
    type G2;

    type G1Precomp;
    type G2Precomp;

    type AffineAteG1Precomp;
    type AffineAteG2Precomp;

    type Fq;
    type Fqe;
    type Fqk;
    type GT;

    // ============================================================
    // Initialization
    // ============================================================

    fn init_public_params();

    // ============================================================
    // Final Exponentiation
    // ============================================================

    fn final_exponentiation(
        elt: &Self::Fqk,
    ) -> Self::GT;

    // ============================================================
    // G1 / G2 Precomputation
    // ============================================================

    fn precompute_g1(
        p: &Self::G1,
    ) -> Self::G1Precomp;

    fn precompute_g2(
        q: &Self::G2,
    ) -> Self::G2Precomp;

    // ============================================================
    // Miller Loop
    // ============================================================

    fn miller_loop(
        prec_p: &Self::G1Precomp,
        prec_q: &Self::G2Precomp,
    ) -> Self::Fqk;

    // ============================================================
    // Affine Ate Precomputation
    // ============================================================

    fn affine_ate_precompute_g1(
        p: &Self::G1,
    ) -> Self::AffineAteG1Precomp;

    fn affine_ate_precompute_g2(
        q: &Self::G2,
    ) -> Self::AffineAteG2Precomp;

    // ============================================================
    // Affine Ate Miller Loop
    // ============================================================

    fn affine_ate_miller_loop(
        prec_p: &Self::AffineAteG1Precomp,
        prec_q: &Self::AffineAteG2Precomp,
    ) -> Self::Fqk;

    fn affine_ate_e_over_e_miller_loop(
        prec_p1: &Self::AffineAteG1Precomp,
        prec_q1: &Self::AffineAteG2Precomp,
        prec_p2: &Self::AffineAteG1Precomp,
        prec_q2: &Self::AffineAteG2Precomp,
    ) -> Self::Fqk;

    fn affine_ate_e_times_e_over_e_miller_loop(
        prec_p1: &Self::AffineAteG1Precomp,
        prec_q1: &Self::AffineAteG2Precomp,
        prec_p2: &Self::AffineAteG1Precomp,
        prec_q2: &Self::AffineAteG2Precomp,
        prec_p3: &Self::AffineAteG1Precomp,
        prec_q3: &Self::AffineAteG2Precomp,
    ) -> Self::Fqk;

    // ============================================================
    // Double Miller Loop
    // ============================================================

    fn double_miller_loop(
        prec_p1: &Self::G1Precomp,
        prec_q1: &Self::G2Precomp,
        prec_p2: &Self::G1Precomp,
        prec_q2: &Self::G2Precomp,
    ) -> Self::Fqk;

    // ============================================================
    // Pairing
    // ============================================================

    fn pairing(
        p: &Self::G1,
        q: &Self::G2,
    ) -> Self::Fqk;

    fn reduced_pairing(
        p: &Self::G1,
        q: &Self::G2,
    ) -> Self::GT;

    fn affine_reduced_pairing(
        p: &Self::G1,
        q: &Self::G2,
    ) -> Self::GT;
}


// pub struct Bn254Params;

// impl PublicParams for Bn254Params {
//     type Fr = Fr;
//     type G1 = G1;
//     type G2 = G2;

//     type G1Precomp = G1Precomp;
//     type G2Precomp = G2Precomp;

//     type AffineAteG1Precomp = AffineAteG1Precomp;
//     type AffineAteG2Precomp = AffineAteG2Precomp;

//     type Fq = Fp;
//     type Fqe = Fp2;
//     type Fqk = Fp12;
//     type GT = Fp12;

//     // 具体方法由 BN254 配对模块实现。
// }