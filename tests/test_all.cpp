#include <cassert>
#include <iostream>
#include "params.hpp"
#include "fp.hpp"
#include "fp2.hpp"
#include "curve.hpp"
#include "isogeny.hpp"
#include "key_exchange.hpp"

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

static void test_isogeny3() {
    MontgomeryCurve E(INITIAL_A);

    for (const ProjectivePoint& G : {P_B, Q_B, R_B}) {
        ProjectivePoint K = E.xMUL(G, 3);
        Isogeny3 iso(K);
        MontgomeryCurve E2(iso.get_next_A());
        ProjectivePoint phiG = iso.eval(G);

        assert(iso.eval(K).is_infinity());                          // 核は無限遠点に写る
        assert(E2.xMUL(phiG, 3).is_infinity());                     // phi(G) の位数は 3
        assert(!phiG.is_infinity());                                // 1 ではない
        assert(proj_eq(iso.eval(E.xDBL(G)), E2.xDBL(phiG)));        // phi(2G) = 2 phi(G)
        assert(proj_eq(iso.eval(E.xTPL(G)), E2.xTPL(phiG)));        // phi(3G) = 3 phi(G)
        assert(proj_eq(iso.eval(E.xMUL(G, 4)), E2.xMUL(phiG, 4)));  // phi(4G) = 4 phi(G)
    }
    std::cout << "[OK] isogeny3\n";
}

static void test_iso_chain_3e() {
    MontgomeryCurve E(INITIAL_A);
    ProjectivePoint A0(INITIAL_A, Fp2(1));

    for (const ProjectivePoint& S : {P_B, Q_B, R_B}) {
        IsogenyChainResult r = iso_chain_3e(A0, S, E_B, P_A, Q_A, R_A);
        MontgomeryCurve F(r.final_A);

        assert(F.xMUL(r.phi_P, 8).is_infinity());    // 位数 8 のまま
        assert(!F.xMUL(r.phi_P, 4).is_infinity());

        for (uint64_t k = 0; k < 8; ++k) {            // 準同型性
            ProjectivePoint X = E.LADDER3PT(P_A, Q_A, R_A, k);   // P + kQ
            IsogenyChainResult rx = iso_chain_3e(A0, S, E_B, X, Q_A, R_A);
            assert(proj_eq(rx.phi_P, F.LADDER3PT(r.phi_P, r.phi_Q, r.phi_R, k)));
        }
    }
    std::cout << "[OK] iso_chain_3e\n";
}

// Alice 側: 核 S = P_A + [m]Q_A (位数 2^3) で 2^e 同種写像を作り、Bob の点を写す
static void test_iso_chain_2e() {
    MontgomeryCurve E(INITIAL_A);
    ProjectivePoint A0(INITIAL_A, Fp2(1));

    for (uint64_t m = 0; m < 8; ++m) {
        ProjectivePoint S = E.LADDER3PT(P_A, Q_A, R_A, m);   // P_A + [m]Q_A
        IsogenyChainResult r = iso_chain_2e(A0, S, E_A, P_B, Q_B, R_B);
        MontgomeryCurve F(r.final_A);

        assert(F.xMUL(r.phi_P, 9).is_infinity());             // 位数 9 のまま
        assert(!F.xMUL(r.phi_P, 3).is_infinity());

        for (uint64_t k = 0; k < 9; ++k) {                    // 準同型性
            ProjectivePoint X = E.LADDER3PT(P_B, Q_B, R_B, k);   // P_B + kQ_B
            IsogenyChainResult rx = iso_chain_2e(A0, S, E_A, X, Q_B, R_B);
            assert(proj_eq(rx.phi_P, F.LADDER3PT(r.phi_P, r.phi_Q, r.phi_R, k)));
        }
    }
    std::cout << "[OK] iso_chain_2e\n";
}

static void test_j_invariant() {
    ProjectivePoint a(INITIAL_A, Fp2(1));
    ProjectivePoint a5(INITIAL_A * Fp2(5), Fp2(5));    // (5A : 5)、同じ曲線

    assert(calc_j_invariant(ProjectivePoint(Fp2(0), Fp2(1))) == Fp2(24));   // 1728 mod 71
    assert(calc_j_invariant(a) == calc_j_invariant(a5));                    // 射影のスケールに依らない
    assert(calc_j_invariant(a) == calc_j_invariant(ProjectivePoint(-INITIAL_A, Fp2(1))));   // A と -A は同じ j
    std::cout << "[OK] j-invariant\n";
}

static void test_key_exchange() {
    for (uint64_t m_A = 0; m_A < 8; ++m_A) {
        for (uint64_t m_B = 0; m_B < 9; ++m_B) {
            PublicKey pkA = alice_keygen(m_A);
            PublicKey pkB = bob_keygen(m_B);
            assert(alice_shared(m_A, pkB) == bob_shared(m_B, pkA));
        }
    }
    std::cout << "[OK] key exchange\n";
}

int main() {
    test_montgomery_ladder();
    test_tripling();
    test_ladder3pt();
    test_isogeny2();
    test_isogeny3();
    test_iso_chain_3e();
    test_iso_chain_2e();
    test_j_invariant();
    test_key_exchange();
    return 0;
}