#include "isogeny.hpp"
#include <cassert> // assertマクロを使用するために追加

Isogeny2::Isogeny2(const ProjectivePoint& K) : K(K) {}

Fp2 Isogeny2::get_next_A() const{
    Fp2 XK2 = K.X * K.X;
    Fp2 ZK2 = K.Z * K.Z;

    Fp2 num = (XK2 + ZK2) * Fp2(2,0);
    Fp2 den = ZK2 - XK2;
    return num * den.inv();
}


ProjectivePoint Isogeny2::eval(const ProjectivePoint& Q) const{
    Fp2 t0 = Q.X + Q.Z;
    Fp2 t1 = Q.X - Q.Z;
    Fp2 t2 = K.X + K.Z;
    Fp2 t3 = K.X - K.Z;
    
    t0 = t0 * t3;
    t1 = t1 * t2;
    t2 = t0 + t1;
    t3 = t0 - t1;
    
    t2 = t2 * t2; 
    t3 = t3 * t3;
    
    Fp2 next_X = Q.X * t2;
    Fp2 next_Z = Q.Z * t3;

    return ProjectivePoint(next_X, next_Z);
}

Isogeny3::Isogeny3(const ProjectivePoint& K) : K(K) {}

// べルーの公式で射影座標における3-isogeny曲線変換
Fp2 Isogeny3::get_next_A(const Fp2& current_A) const {
    // アフィン相当の式: A' = (A * xK - 6 * (xK^2 - 1)) * xK
    Fp2 xK = K.X * K.Z.inv(); // xK = XK / ZK
    Fp2 xK2 = xK * xK;

    Fp2 term1 = current_A * xK;
    Fp2 term2 = (xK2 - Fp2(1, 0)) * Fp2(6, 0);

    return (term1 - term2) * xK;
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
    const Fp2& start_A, 
    ProjectivePoint S, 
    int e, 
    ProjectivePoint P, 
    ProjectivePoint Q, 
    ProjectivePoint R
) {
    Fp2 current_A = start_A;
    MontgomeryCurve curve(current_A);

    for (int i = e - 1; i >= 0; --i) {
        // 1. スカラー倍算で現在の S から "位数 2" の点 K を抽出
        ProjectivePoint K = S;
        for (int j = 0; j < i; ++j) {
            K = curve.xDBL(K);
        }

        // --- 実行時アサーション: 核 K の位数2検証 ---
        assert(!K.Z.is_zero() && "Fatal: 2-torsion kernel point K degenerated to infinity.");
        ProjectivePoint K_times_2 = curve.xDBL(K);
        assert(K_times_2.Z.is_zero() && "Fatal: Kernel point K does not have strictly order 2.");
        // -------------------------------------------
        // iso_chain_2e のループ内


        // --- デバッグ出力 ---
        std::cout << "[DEBUG] Step i=" << i << " K.X="; K.X.print();
        std::cout << "[DEBUG] Step i=" << i << " K.Z="; K.Z.print();
        
        std::cout << "[DEBUG] K_times_2.X="; K_times_2.X.print();
        std::cout << "[DEBUG] K_times_2.Z="; K_times_2.Z.print();
        // --------------------

        assert(!K.Z.is_zero() && "Fatal: 2-torsion kernel point K degenerated to infinity.");
        assert(K_times_2.Z.is_zero() && "Fatal: Kernel point K does not have strictly order 2.");
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
    const Fp2& start_A, 
    ProjectivePoint S, 
    int e, 
    ProjectivePoint P, 
    ProjectivePoint Q, 
    ProjectivePoint R
) {
    Fp2 current_A = start_A;
    MontgomeryCurve curve(current_A);

    for (int i = e - 1; i >= 0; --i) {
        ProjectivePoint K = S;
        for (int j = 0; j < i; ++j) {
            // [3]K = [2]K + K の計算 (差分点は元の K)
            ProjectivePoint K2 = curve.xDBL(K);
            // K = curve.xDBLADD(K2, K, K);
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