#include "curve.hpp"

ProjectivePoint MontgomeryCurve::xDBL(const ProjectivePoint& P) const {

    Fp2 add = P.X + P.Z;
    Fp2 sub = P.X - P.Z;

    Fp2 sq_add = add * add;
    Fp2 sq_sub = sub * sub;

    Fp2 four_XZ = sq_add - sq_sub;

    Fp2 next_X = sq_add * sq_sub;
    Fp2 next_Z = four_XZ * (sq_sub + A24 * four_XZ.inv());

    return ProjectivePoint(next_X, next_Z);
}

ProjectivePoint MontgomeryCurve::xADD(const ProjectivePoint& P, const ProjectivePoint& Q, const ProjectivePoint& P_minus_Q) const {
    
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

ProjectivePoint MontgomeryCurve::xMUL(const ProjectivePoint& P, uint64_t k) const {

    if(k == 0) return ProjectivePoint(Fp2(0,0), Fp2(1,0));
    
    // R0:無限遠点 R1:入力点
    ProjectivePoint R0 = ProjectivePoint(Fp2(0,0), Fp2(1,0));
    ProjectivePoint R1 = P;

    // 高速化のため最上位ビットを探す
    // サイドチャネルには弱い
    int msb = 63;
    while(msb >= 0 && (!(k >> msb)) & 1){
        msb--;
    }
    // モンゴメリーラダー
    // R0が無限遠点の場合、加算がうまくいかないので一手目は手書きで行う
    // k >= 1の時常に一手目は1が立っているので処理の分岐はない
    R0 = P;
    R1 = xDBL(P);
    for(int i = msb-1; i >= 0; i--){
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