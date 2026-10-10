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

// xTPL が xMUL(P, 3) と一致し、位数 3^2 の点は 2 回で無限遠点になるか確認する
static void test_tripling() {
    MontgomeryCurve curve(INITIAL_A);

    for (const ProjectivePoint& P : {P_A, Q_A, R_A, P_B, Q_B, R_B}) {
        assert(proj_eq(curve.xTPL(P), curve.xMUL(P, 3)));
    }
    for (const ProjectivePoint& P : {P_B, Q_B, R_B}) {
        assert(curve.xTPL(curve.xTPL(P)).is_infinity());
    }

    std::cout << "[OK] tripling\n";
}

static void test_ladder3pt() {
    MontgomeryCurve curve(INITIAL_A);
    const ProjectivePoint G = P_B;

    for (uint64_t k = 0; k < 20; ++k) {
        ProjectivePoint P = curve.xMUL(G, 5);
        ProjectivePoint Q = curve.xMUL(G, 2);
        ProjectivePoint D = curve.xMUL(G, 3);   // P - Q
        assert(proj_eq(curve.LADDER3PT(P, Q, D, k), curve.xMUL(G, 5 + 2 * k)));
    }
    std::cout << "[OK] ladder3pt\n";
}

static void test_isogeny2() {
    MontgomeryCurve E(INITIAL_A);

    for (const ProjectivePoint& G : {P_A, R_A}) {      // Q_A は核が (0:1) になるので使わない
        ProjectivePoint K = E.xMUL(G, 4);               // 位数 2 の核
        Isogeny2 iso(K);
        MontgomeryCurve E2(iso.get_next_A());
        ProjectivePoint phiG = iso.eval(G);

        assert(iso.eval(K).is_infinity());
        assert(E2.xMUL(phiG, 4).is_infinity());
        assert(!E2.xMUL(phiG, 2).is_infinity());
        assert(proj_eq(iso.eval(E.xDBL(G)), E2.xDBL(phiG)));
        assert(proj_eq(iso.eval(E.xMUL(G, 3)), E2.xMUL(phiG, 3)));
    }
    std::cout << "[OK] isogeny2\n";
}

int main() {
    test_montgomery_ladder();
    test_tripling();
    test_ladder3pt();
    test_isogeny2();
    return 0;
}