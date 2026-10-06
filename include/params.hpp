#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <cstdint>
#include "fp2.hpp"
#include "curve.hpp"

// 曲座標のパラメータ
namespace ToyParams {
    // 基本パラメータ
    const uint64_t P_PRIME = 71;
    const int E_A = 3;
    const int E_B = 2;
    const Fp2 INITIAL_A(6, 0);

    // Alice の基底点 (位数 2^eA)
    const ProjectivePoint P_A(Fp2(69, 22), Fp2(1, 0));
    const ProjectivePoint Q_A(Fp2(25, 21), Fp2(1, 0));
    const ProjectivePoint R_A(Fp2(48, 59), Fp2(1, 0));

    // Bob の基底点 (位数 3^eB)
    const ProjectivePoint P_B(Fp2(4, 40), Fp2(1, 0));
    const ProjectivePoint Q_B(Fp2(12, 58), Fp2(1, 0));
    const ProjectivePoint R_B(Fp2(0, 39), Fp2(1, 0));
}

#endif // PARAMS_HPP