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