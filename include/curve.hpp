#ifndef CURVE_HPP
#define CURVE_HPP

#include "fp.hpp"
#include "fp2.hpp"
#include <iostream>

// アフィン座標上の点を表す構造体
struct AffinePoint {
    Fp2 X;
    Fp2 Y;
    bool is_infinity;

    // コンストラクタ
    AffinePoint(const Fp2& x, const Fp2& y)
        : X(x), Y(y), is_infinity(false) {}
    // デフォルトの場合は無限遠点を返す
    AffinePoint() : X(0), Y(0), is_infinity(true) {}

    // 無限遠点を返す関数
    static AffinePoint Infinity() {
        return AffinePoint();
    }

    void print() const {
        std::cout << "x : ";
        X.print();
        std::cout << "y : ";
        Y.print();
    }
};

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


    // アフィン座標上の演算
    // xADD : 加算, xDBL : 2倍算, xMUL : k倍算
    AffinePoint xADD(const AffinePoint& P, const AffinePoint& Q) const ;
    AffinePoint xDBL(const AffinePoint& P) const ;
    AffinePoint xMUL(const AffinePoint& P, uint64_t k) const ;


    // 射影座標上の演算
    // P_minus_Qは最初に入力された点を入れる
    ProjectivePoint xADD(const ProjectivePoint& P, const ProjectivePoint& Q, const ProjectivePoint& P_minus_Q) const ;
    ProjectivePoint xDBL(const ProjectivePoint& P) const ;
    ProjectivePoint xMUL(const ProjectivePoint& P, uint64_t k) const ;
};

// j不変量の計算
Fp2 calc_j_invariant(const Fp2& A);
#endif // CURVE_HPP