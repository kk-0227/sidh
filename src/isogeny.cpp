#include "isogeny.hpp"
#include <cassert> // assertマクロを使用するために追加

Isogeny2::Isogeny2(const ProjectivePoint& K) : K(K) {
    assert(!K.X.is_zero() && "2-isogeny kernel must not be (0:1)");
}

ProjectivePoint Isogeny2::get_next_A() const{
    // (A' : C') = (2(Z^2 - 2X^2) : Z^2)
    // A' = 2(Z^2 - 2X^2)/Z^2 (Affine)
    Fp2 t0 = K.Z * K.Z;

    Fp2 next_A = Fp2(2) * (t0 - Fp2(2) * K.X * K.X);
    Fp2 next_C = t0;

    return ProjectivePoint(next_A, next_C);
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
ProjectivePoint Isogeny3::get_next_A() const {
    // (A' : C') = (Z^4 + 18 * X^2 * Z^2 - 27 * X^4 : 4 * X * Z^3)

    Fp2 t0 = K.X * K.X;
    Fp2 t1 = K.Z * K.Z;

    Fp2 next_A = t1 * t1 + Fp2(18) * t0 * t1 - Fp2(27) * t0 * t0;
    Fp2 next_C = Fp2(4) * K.X * K.Z * t1;

    return ProjectivePoint(next_A, next_C);
}

// 点Pを同種写像で写す
ProjectivePoint Isogeny3::eval(const ProjectivePoint& P) const {
    Fp2 t0 = K.X * P.X - K.Z * P.Z;
    Fp2 t1 = K.Z * P.X - K.X * P.Z;

    Fp2 next_X = P.X * t0 * t0;
    Fp2 next_Z = P.Z * t1 * t1;

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

    // 現在の曲線を生成
    ProjectivePoint current_A = start_A;
    MontgomeryCurve curve(current_A);

    for (int i = e - 1; i >= 0; --i) {
        // K = [3^i]S : 新しい同種写像の核
        ProjectivePoint K = S;

        for (int j = 0; j < i; ++j) {
            K = curve.xTPL(K);
        }
        // 写像を生成
        Isogeny3 iso(K);

        // 写像により移った先の曲線
        current_A = iso.get_next_A();
        curve = MontgomeryCurve(current_A);


        S = iso.eval(S);
        P = iso.eval(P);
        Q = iso.eval(Q);
        R = iso.eval(R);
    }

    return {current_A, P, Q, R};
}