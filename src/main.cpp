#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <string>
#include "params.hpp"
// #include "params_affine.hpp"
#include "fp.hpp"
#include "fp2.hpp"
#include "curve.hpp"
#include "isogeny.hpp"
using namespace ToyParams;

void print_affine(ProjectivePoint P) {
    Fp2 P_af = P.X*P.Z.inv();
    P_af.print();
}

int main() {
    //AffinePoint P = P_A, Q = Q_A;
    ProjectivePoint P = P_A, Q = Q_A, R = R_A;
    MontgomeryCurve curve(INITIAL_A);
    // アリスの秘密点(位数8)の計算(P+mQ)
    // m = 1とした
    ProjectivePoint S_A = curve.xADD(P, Q, R);
    // 2-同種写像用の位数2の点を生成
    ProjectivePoint K = curve.xMUL(S_A, 4); // S_Aの4倍点を計算してKとする
    Isogeny2 iso(K);
    //IsogenyChainResult iso_result = iso_chain_2e(INITIAL_A, S_A, E_A, P, Q, R);
    

    //std::cout << "S_A'" << std::endl;
    
    std::cout << "K" << std::endl;
    print_affine(K);
    std::cout << "A' = ";
    iso.get_next_A().print();
    std::cout << "P_B'" << std::endl;
    iso.eval(P).print();
    std::cout << "Q_B'" << std::endl;
    iso.eval(Q).print();
    std::cout << "R_B'" << std::endl;
    iso.eval(R).print();

    return 0;
}