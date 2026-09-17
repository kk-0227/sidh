#include <iostream>
#include "fp.hpp"
#include "fp2.hpp"
#include "curve.hpp"

void print_affine(const ProjectivePoint& P){
    if(P.Z == Fp2(0,0)){
        std::cout << "無限遠点" << std::endl;
        return ;
    } 
    //　射影座標をアフィン座標に変換 x = X/Z
    Fp2 affine_x = P.X * P.Z.inv();
    affine_x.print(); 
}

int main(){
    // モンゴメリ係数の決定
    Fp2 A(6,0);
    // モンゴメリ曲線の生成
    MontgomeryCurve curve(A);
    // 始点を生成
    // x = 2, z = 1
    ProjectivePoint P(Fp2(2,0), Fp2(1,0));
    // 2倍点の計算
    ProjectivePoint P2_dbl = curve.xDBL(P);
    print_affine(P2_dbl);
    // 無限遠点の生成
    ProjectivePoint O(Fp2(1,0), Fp2(0,0));
    print_affine(O);
    
    ProjectivePoint P3_add = curve.xADD(P, P2_dbl, P);
    ProjectivePoint P3_dbl = curve.xADD(P2_dbl, P, P);
    print_affine(P3_add);
    print_affine(P3_dbl);
    ProjectivePoint P3_mongomer = curve.xMUL(P,3);
    print_affine(P3_mongomer);
    
}