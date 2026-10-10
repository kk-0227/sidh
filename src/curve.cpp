#include "curve.hpp"

// 2倍算
ProjectivePoint MontgomeryCurve::xDBL(const ProjectivePoint& P) const {

    Fp2 t0 = P.X - P.Z;
    Fp2 t1 = P.X + P.Z;

    t0 = t0 * t0;
    t1 = t1 * t1;

    Fp2 next_Z = C24plus * t0;
    Fp2 next_X = next_Z * t1;

    t1 = t1 - t0;
    t0 = A24plus * t1;

    next_Z = next_Z + t0;
    next_Z = next_Z * t1;

    return ProjectivePoint(next_X, next_Z);
}

// 2倍算と加算を同時に行う
DblAddResult MontgomeryCurve::xDBLADD(const ProjectivePoint& P, const ProjectivePoint& Q, const ProjectivePoint& P_minus_Q) const {
    Fp2 t0 = P.X + P.Z;
    Fp2 t1 = P.X - P.Z;

    ProjectivePoint DBL, ADD;
    
    DBL.X = t0 * t0;
    ADD.X = Q.X + Q.Z;
    Fp2 t2 = Q.X - Q.Z;

    t0 = t0 * t2;
    DBL.Z = t1 * t1;
    t1 = t1 * ADD.X;

    t2 = DBL.X - DBL.Z;

    DBL.X = DBL.X * DBL.Z;
    ADD.X = A24 * t2;
    ADD.Z = t0 - t1;
    DBL.Z = ADD.X + DBL.Z;
    ADD.X = t0 + t1;
    DBL.Z = DBL.Z * t2;
    ADD.Z = ADD.Z * ADD.Z;
    ADD.X = ADD.X * ADD.X;
    ADD.Z = ADD.Z * P_minus_Q.X;
    ADD.X = ADD.X * P_minus_Q.Z;  

    return DblAddResult(DBL, ADD);
}
// k倍算
ProjectivePoint MontgomeryCurve::xMUL(const ProjectivePoint& P, uint64_t k) const {
    // mPを計算する
    ProjectivePoint R0 = ProjectivePoint::infinity(), R1 = P;
    // ビットが立ってる時は2*R0+1,立ってない時は2*R0とする
    // R1は常にR0+1を保持する
    for (int i = 63; i >= 0; --i) {
        if ((k >> i) & 1) {
            DblAddResult result = xDBLADD(R1, R0, P);
            R0 = result.ADD;
            R1 = result.DBL;
        } else {
            DblAddResult result = xDBLADD(R0, R1, P);
            R0 = result.DBL;
            R1 = result.ADD;
        }
    }
    return R0;
}

// 3倍算
ProjectivePoint MontgomeryCurve::xTPL(const ProjectivePoint& P) const {
    
    Fp2 t0 = P.X - P.Z;
    Fp2 t1 = P.X + P.Z;

    Fp2 t2 = t0 * t0;
    Fp2 t3 = t1 * t1;
    Fp2 t4 = Fp2(2) * P.X;

    t0 = Fp2(2) * P.Z;
    t1 = t4 * t4;
    t1 = t1 - t3 - t2;

    Fp2 t5 = A24plus * t3;
    t3 = t5 * t3;

    Fp2 t6 = A24minus * t2;
    t2 = t6 * t2;
    t3 = t2 - t3;

    t2 = t5 - t6;
    t1 = t1*t2;

    t2 = t1 + t3;
    t2 = t2 * t2;

    Fp2 next_X = t2 * t4;

    t1 = t3 - t1;
    t1 = t1 * t1;

    Fp2 next_Z = t0 * t1;

    return ProjectivePoint(next_X, next_Z);
}

// P + [k]Q を計算
ProjectivePoint MontgomeryCurve::LADDER3PT(const ProjectivePoint& P, const ProjectivePoint& Q, const ProjectivePoint& P_minus_Q, uint64_t k) const {
    ProjectivePoint R0 = Q; // [2^i]Q
    ProjectivePoint R1 = P; // P + [kの下位iビット]Q
    ProjectivePoint R2 = P_minus_Q; // 差分 R1 - R0

    // R0は常に2倍したい
    // ビットが立っている時は　R1 = R0 + R1, 差分(R2)は更新なし
    // ビットが立っていない時は R1 は更新なし, R2 = R2 - R0(R0が増えた分差分が縮む)
    for( ; k; k >>= 1) {
        if(k & 1) {
            DblAddResult r = xDBLADD(R0, R1, R2);   // R1 + R0(差分は R2)
            R0 = r.DBL;
            R1 = r.ADD;
        } else {
            DblAddResult r = xDBLADD(R0, R2, R1);   // R2 + R0(差分は R1)
            R0 = r.DBL;
            R2 = r.ADD;
        }
    }
    
    return R1;
}

Fp2 calc_j_invariant(const ProjectivePoint& A) {
    // j = 256 (A^2 - 3C^2)^3 / ( C^4 (A^2 - 4C^2) )
    Fp2 t0 = A.X * A.X;
    Fp2 t1 = A.Z * A.Z;
    Fp2 t2 = (t0 - Fp2(3) * t1);

    Fp2 num = Fp2(256)  * t2 * t2 * t2;
    Fp2 den = t1 * t1 * (t0 - Fp2(4) * t1);

    return num * den.inv();
}