#ifndef PARAMS_AFFINE_HPP
#define PARAMS_AFFINE_HPP

#include <cstdint>
#include "fp2.hpp"
#include "curve.hpp"

// アフィン座標のパラメータ
const uint64_t P_PRIME = 71;
const int E_A = 3;
const int E_B = 2;
const Fp2 INITIAL_A(6, 0);

// Alice's Basis Points (Order 2^eA)
const AffinePoint P_A(Fp2(69, 22), Fp2(70, 49));
const AffinePoint Q_A(Fp2(25, 21), Fp2(56, 28));
const AffinePoint R_A(Fp2(48, 59), Fp2(63, 24));

// Bob's Basis Points (Order 3^eB)
const AffinePoint P_B(Fp2(4, 40), Fp2(45, 42));
const AffinePoint Q_B(Fp2(12, 58), Fp2(68, 32));
const AffinePoint R_B(Fp2(0, 39), Fp2(34, 22));

#endif // PARAMS_AFFINE_HPP
