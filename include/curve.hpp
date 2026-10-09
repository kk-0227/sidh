#ifndef CURVE_HPP
#define CURVE_HPP

#include "fp.hpp"
#include "fp2.hpp"
#include <iostream>

// 射影座標上の点を表す構造体
struct ProjectivePoint {
    Fp2 X;
    Fp2 Z;
    ProjectivePoint() : X(0), Z(1) {};
    ProjectivePoint(const Fp2& x, const Fp2& z) : X(x), Z(z) {}

    void print() const {
        std::cout << "X : ";
        X.print();
        std::cout << "Z : ";
        Z.print();
    }

    static ProjectivePoint infinity() {
        return ProjectivePoint(Fp2(1), Fp2(0));
    }

    bool is_infinity() const {
        return Z == Fp2(0);
    }
};

struct DblAddResult {
    ProjectivePoint DBL;
    ProjectivePoint ADD;

    DblAddResult(const ProjectivePoint& dbl, const ProjectivePoint& add) : DBL(dbl), ADD(add) {}
};

class MontgomeryCurve {
private:

    // 4の逆元を一度だけ計算
    static const Fp2& inv4(){
        static const Fp2 inv4_value = Fp2(4).inv();
        return inv4_value;
    }
    // モンゴメリ係数(Cy^2 = x^3 + Ax^2 + x)
    // A = (A : C)
    ProjectivePoint A;

    // 何度も使うので前計算しておく
    Fp2 A24plus; // A + 2C
    Fp2 C24plus; // 4C
    Fp2 A24minus; // A - 2C
    Fp2 A24; // (A + 2C)/(4C)

public:

    // 曲線の生成
    // 初期値は(A : C) = (A : 1)と考える
    MontgomeryCurve(const Fp2& a) : A(ProjectivePoint(a, Fp2(1))), A24plus(a + Fp2(2)), C24plus(Fp2(4)), A24minus(a - Fp2(2)), A24(A24plus * inv4()) {}
    MontgomeryCurve(const ProjectivePoint& a) : A(a), A24plus(a.X + Fp2(2) * a.Z), C24plus(a.Z * Fp2(4)), A24minus(a.X - Fp2(2) * a.Z), A24(A24plus * C24plus.inv()) {}
    
    ProjectivePoint get_A() const { return A; }

    // 射影座標上の演算
    // P_minus_Qは最初に入力された点を入れる
    ProjectivePoint xDBL(const ProjectivePoint& P) const ;
    DblAddResult xDBLADD(const ProjectivePoint& P, const ProjectivePoint& Q, const ProjectivePoint& P_minus_Q) const ;
    ProjectivePoint xMUL(const ProjectivePoint& P, uint64_t k) const ;
    ProjectivePoint xTPL(const ProjectivePoint& P) const ;
    ProjectivePoint LADDER3PT(const ProjectivePoint& P, const ProjectivePoint& Q, const ProjectivePoint& P_minus_Q, uint64_t k) const ;
};

// j不変量の計算
Fp2 calc_j_invariant(const Fp2& A);
#endif // CURVE_HPP