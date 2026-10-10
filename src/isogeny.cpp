#include "isogeny.hpp"
#include <cassert> // assertマクロを使用するために追加

Isogeny2::Isogeny2(const ProjectivePoint& K) : K(K) {
    assert(!K.X.is_zero() && "2-isogeny kernel must not be (0:1)");
}

ProjectivePoint Isogeny2::get_next_A() const{
    // (A' : C') = (2(Z^2 - 2X^2) : Z^2)
    // A' = 2(Z^2 - 2X^2)/Z^2 (Affine)
    Fp2 t0 = K.Z * K.Z;

    return ProjectivePoint(Fp2(2) * (t0 - Fp2(2) * K.X * K.X), t0);
}

ProjectivePoint Isogeny2::eval(const ProjectivePoint& P) const{
    
    Fp2 t0 = (K.X + K.Z) * (P.X - P.Z);
    Fp2 t1 = (K.X - K.Z) * (P.X + P.Z);

    Fp2 next_X = P.X * (t0 + t1);
    Fp2 next_Z = P.Z * (t0 - t1);
    
    return ProjectivePoint(next_X, next_Z);
}

Isogeny3::Isogeny3(const ProjectivePoint& K) : K(K) {}

// べルーの公式で射影座標における3-isogeny曲線変換
ProjectivePoint Isogeny3::get_next_A(const ProjectivePoint& current_A) const {
    
}

// 点Qを同種写像で写す
ProjectivePoint Isogeny3::eval(const ProjectivePoint& Q) const {
    Fp2 XQ_XK = Q.X * K.X;
    Fp2 ZQ_ZK = Q.Z * K.Z;
    Fp2 XQ_ZK = Q.X * K.Z;
    Fp2 ZQ_XK = Q.Z * K.X;

    // term1 = XQ * XK - ZQ * ZK
    Fp2 term1 = XQ_XK - ZQ_ZK;
    term1 = term1 * term1; // (XQ*XK - ZQ*ZK)^2

    // term2 = XQ * ZK - ZQ * XK
    Fp2 term2 = XQ_ZK - ZQ_XK;
    term2 = term2 * term2; // (XQ*ZK - ZQ*XK)^2

    Fp2 next_X = Q.X * term1;
    Fp2 next_Z = Q.Z * term2;

    return ProjectivePoint(next_X, next_Z);
}

IsogenyChainResult iso_chain_2e(
    const ProjectivePoint& start_A, 
    ProjectivePoint S, 
    int e, 
    ProjectivePoint P, 
    ProjectivePoint Q, 
    ProjectivePoint R
) {
    ProjectivePoint current_A = start_A;
    MontgomeryCurve curve(current_A);

    for (int i = e - 1; i >= 0; --i) {
        // 1. スカラー倍算で現在の S から "位数 2" の点 K を抽出
        ProjectivePoint K = S;
        for (int j = 0; j < i; ++j) {
            K = curve.xDBL(K);
        }

        // 2. K を核として 2-isogeny を構築
        Isogeny2 iso(K);
        current_A = iso.get_next_A();
        curve = MontgomeryCurve(current_A); // 曲線更新

        // 3. 次のステップのために点 S, P, Q, R を新しい曲線へ写像
        S = iso.eval(S);
        P = iso.eval(P);
        Q = iso.eval(Q);
        R = iso.eval(R);
    }

    return {current_A, P, Q, R};
}

IsogenyChainResult iso_chain_3e(
    const ProjectivePoint& start_A, 
    ProjectivePoint S, 
    int e, 
    ProjectivePoint P, 
    ProjectivePoint Q, 
    ProjectivePoint R
) {
    ProjectivePoint current_A = start_A;
    MontgomeryCurve curve(current_A);

    for (int i = e - 1; i >= 0; --i) {
        ProjectivePoint K = S;
        for (int j = 0; j < i; ++j) {
            K = curve.xDBL(K);
        }

        // --- 実行時アサーション: 核 K の位数3検証 ---
        //assert(!K.Z.is_zero() && "Fatal: 3-torsion kernel point K degenerated to infinity.");
        //ProjectivePoint K_times_2 = curve.xDBL(K);
        //ProjectivePoint K_times_3 = curve.xDBLADD(K_times_2, K, K); // 差分点は K
        //assert(K_times_3.Z.is_zero() && "Fatal: Kernel point K does not have strictly order 3.");
        // -------------------------------------------

        Isogeny3 iso(K);
        current_A = iso.get_next_A(current_A);
        curve = MontgomeryCurve(current_A);

        S = iso.eval(S);
        P = iso.eval(P);
        Q = iso.eval(Q);
        R = iso.eval(R);
    }

    return {current_A, P, Q, R};
}