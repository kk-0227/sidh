#include "fp.hpp"
#include "fp2.hpp"
#include <iostream>

// 射影座標上の点を表す構造体
struct ProjectivePoint {
    Fp2 X;
    Fp2 Z;
    ProjectivePoint() : X(0), Z(1) {};
    ProjectivePoint(const Fp2& x, const Fp2& z) : X(x), Z(z) {}
};

class MontgomeryCurve {
private:
    // モンゴメリ係数(By^2 = x^3 + Ax^2 + x)
    Fp2 A;
    // (A+2)/4
    // 何度も使うので前計算しておく
    Fp2 A24;
public:
    MontgomeryCurve(const Fp2& a) : A(a), A24((a + 2)*Fp2(4).inv()) {}

    Fp2 get_A() const { return A; }
    // 2倍算
    ProjectivePoint xDBL(const ProjectivePoint& P) const ; 
    // 加算
    // P_minus_Qは最初に入力された点を入れる
    ProjectivePoint xADD(const ProjectivePoint& P, const ProjectivePoint& Q, const ProjectivePoint& P_minus_Q) const ;
    // スカラー倍算
    ProjectivePoint xMUL(const ProjectivePoint& P, uint64_t k) const ;
};