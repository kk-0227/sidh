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

Fp2 calc_j_invariant(const Fp2& A) {
    Fp2 A2 = A * A;
    Fp2 num = A2 - Fp2(3, 0);
    num = num * num * num * Fp2(256, 0); // 256 * (A^2 - 3)^3

    Fp2 den = A2 - Fp2(4, 0); // A^2 - 4

    return num * den.inv();
}