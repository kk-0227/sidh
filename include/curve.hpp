#include "fp.hpp"
#include "fp2.hpp"
#include <iostream>

struct ProjectivePoint {
    Fp2 X;
    Fp2 Z;
    ProjectivePoint() : X(0), Z(1) {};
};

class MontgomeryCurve {
private:
    //モンゴメリ曲線(By^2 = x^3 + Ax^2 + x)で演算に使う係数
    Fp2 A;
public:

};