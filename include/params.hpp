#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <cstdint>
#include "fp2.hpp"
#include "curve.hpp"

namespace ToyParams {
    const uint64_t P_PRIME = 71;
    const int E_A = 3;
    const int E_B = 2;
    const Fp2 INITIAL_A(6, 0);

    // Alice 基底点 (位数 2^3)
    const ProjectivePoint P_A(Fp2(19, 66), Fp2(1, 0));
    const ProjectivePoint Q_A(Fp2(19, 62), Fp2(1, 0));
    const ProjectivePoint R_A(Fp2(62, 26), Fp2(1, 0));

    // Bob 基底点 (位数 3^2)
    const ProjectivePoint P_B(Fp2(45, 0), Fp2(1, 0));
    const ProjectivePoint Q_B(Fp2(17, 0), Fp2(1, 0));
    const ProjectivePoint R_B(Fp2(52, 16), Fp2(1, 0));
}

#endif // PARAMS_HPP