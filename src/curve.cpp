#include "curve.hpp"

AffinePoint MontgomeryCurve::xADD(const AffinePoint& P, const AffinePoint& Q) const {
    if(P.is_infinity)return Q;
    if(Q.is_infinity)return P;

    if(P.X == Q.X){
        if(P.Y == -Q.Y){
            return AffinePoint::Infinity();
        }
        if(P.Y == Q.Y){
            return xDBL(P);
        }
    }
    Fp2 tmp = (Q.Y - P.Y) * (Q.X - P.X).inv();
    AffinePoint next_point;
    next_point.is_infinity = false;
    next_point.X = tmp * tmp - A - P.X - Q.X;
    next_point.Y = tmp * (P.X - next_point.X) - P.Y;
    return next_point;
}

AffinePoint MontgomeryCurve::xDBL(const AffinePoint& P) const {
    Fp2 num = Fp2(3) * P.X * P.X + Fp2(2)* A * P.X + Fp2(1);
    Fp2 den = Fp2(2) * P.Y;
    Fp2 lambda = num * den.inv();
    AffinePoint next_point;
    next_point.is_infinity = false;
    next_point.X = lambda * lambda - A - Fp2(2) * P.X;
    next_point.Y = lambda * (P.X - next_point.X) - P.Y;
    return next_point;
}

AffinePoint MontgomeryCurve::xMUL(const AffinePoint& P, uint64_t k) const {
    if(k == 0)return AffinePoint::Infinity();
    AffinePoint R0 = AffinePoint::Infinity();
    AffinePoint R1 = P;

    // 最上位ビットを探す
    int msb = 63;
    while(msb >= 0 && !((k >> msb) & 1)){
        msb--;
    }
    for(int i = msb; i >= 0; i--){
        if((k >> i) & 1){
            AffinePoint next_R0 = xADD(R0, R1);
            AffinePoint next_R1 = xDBL(R1);
            R0 = next_R0;
            R1 = next_R1;
        }else{
            AffinePoint next_R1 = xADD(R0, R1);
            AffinePoint next_R0 = xDBL(R0);
            R1 = next_R1;
            R0 = next_R0;
        }
    }
    return R0;
}

ProjectivePoint MontgomeryCurve::xADD(const ProjectivePoint& P, const ProjectivePoint& Q, const ProjectivePoint& P_minus_Q) const {
    
    //if (P.Z.is_zero()) return Q;
    //if (Q.Z.is_zero()) return P;

    //if(P.X == Q.X && P.Z == Q.Z)return xDBL(P);

    Fp2 add_P = P.X + P.Z;
    Fp2 sub_P = P.X - P.Z;
    Fp2 add_Q = Q.X + Q.Z;
    Fp2 sub_Q = Q.X - Q.Z;

    Fp2 add_sub = add_P * sub_Q;
    Fp2 sub_add = sub_P * add_Q;
    

    Fp2 root_next_X = sub_add + add_sub;
    Fp2 root_next_Z = sub_add - add_sub;

    Fp2 next_X = P_minus_Q.Z * root_next_X * root_next_X;
    Fp2 next_Z = P_minus_Q.X * root_next_Z * root_next_Z;
    
    return ProjectivePoint(next_X, next_Z);
}

ProjectivePoint MontgomeryCurve::xDBL(const ProjectivePoint& P) const {

    Fp2 add = P.X + P.Z;
    Fp2 sub = P.X - P.Z;

    Fp2 sq_add = add * add;
    Fp2 sq_sub = sub * sub;

    Fp2 four_XZ = sq_add - sq_sub;

    Fp2 next_X = sq_add * sq_sub;
    // Fp2 next_Z = four_XZ * (sq_sub + A24 * four_XZ.inv());
    Fp2 next_Z = four_XZ * (sq_sub + A24 * four_XZ);

    return ProjectivePoint(next_X, next_Z);
}

ProjectivePoint MontgomeryCurve::xMUL(const ProjectivePoint& P, uint64_t k) const {

    if(k == 0) return ProjectivePoint(Fp2(1,0), Fp2(0,0));
    
    // R0:無限遠点 R1:入力点
    ProjectivePoint R0 = ProjectivePoint(Fp2(1,0), Fp2(0,0));
    ProjectivePoint R1 = P;

    // 高速化のため最上位ビットを探す
    // サイドチャネルには弱い
    int msb = 63;
    while(msb >= 0 && !((k >> msb) & 1)){
        msb--;
    }
    // モンゴメリーラダー
    // R0が無限遠点の場合、加算がうまくいかないので一手目は手書きで行う
    // k >= 1の時常に一手目は1が立っているので処理の分岐はない
    //R0 = P;
    //R1 = xDBL(P);
    for(int i = msb; i >= 0; i--){
        if((k >> i) & 1){
            R0 = xADD(R0, R1, P);
            R1 = xDBL(R1);
        }else{
            R1 = xADD(R0, R1, P);
            R0 = xDBL(R0);
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