#include <cassert>
#include <iostream>
#include "params.hpp"
#include "fp.hpp"
#include "fp2.hpp"
#include "curve.hpp"
#include "isogeny.hpp"

using namespace ToyParams;

// X/Z の比が等しいかで射影点の同値を判定する
static bool proj_eq(const ProjectivePoint& P, const ProjectivePoint& Q) {
    return P.X * Q.Z == Q.X * P.Z;
}

// xDBL / xDBLADD / xMUL が互いに整合しているか確認する
static void test_montgomery_ladder() {
    MontgomeryCurve curve(INITIAL_A);
    const ProjectivePoint P = P_A;

    ProjectivePoint P2 = curve.xDBL(P);
    ProjectivePoint P3 = curve.xDBLADD(P, P2, P).ADD;
    ProjectivePoint P4 = curve.xDBL(P2);

    assert(proj_eq(curve.xMUL(P, 1), P));
    assert(proj_eq(curve.xMUL(P, 2), P2));
    assert(proj_eq(curve.xMUL(P, 3), P3));
    assert(proj_eq(curve.xMUL(P, 4), P4));
    assert(curve.xMUL(P, 0).is_infinity());

    std::cout << "[OK] montgomery ladder\n";
}

int main() {
    test_montgomery_ladder();
    return 0;
}