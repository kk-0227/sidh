#include "key_exchange.hpp"
#include "params.hpp"

using namespace ToyParams;

PublicKey alice_keygen(uint64_t m_A) {

    // 曲線の生成
    MontgomeryCurve curve(INITIAL_A);

    // aliceの秘密点 S = P + [m]Q を生成
    ProjectivePoint S = curve.LADDER3PT(P_A, Q_A, R_A, m_A);

    // bobの基底点をSで移したものを公開鍵として公開する
    return iso_chain_2e(curve.get_A(), S, E_A, P_B, Q_B, R_B);
}

PublicKey bob_keygen(uint64_t m_B) {

    MontgomeryCurve curve(INITIAL_A);

    ProjectivePoint S = curve.LADDER3PT(P_B, Q_B, R_B, m_B);
    
    return iso_chain_3e(curve.get_A(), S, E_B, P_A, Q_A, R_A);
}


Fp2 alice_shared(uint64_t m_A, const PublicKey& bob_pk) {

    // bobの公開鍵から曲線を生成
    MontgomeryCurve curve(bob_pk.final_A);
    
    // bobが写したphi_B(P_A), phi_B(Q_A), phi_B(R_A)から核Sを生成
    ProjectivePoint S = curve.LADDER3PT(bob_pk.phi_P, bob_pk.phi_Q, bob_pk.phi_R, m_A);

    // 写す点は使わないので、ダミーを渡す
    PublicKey ans = iso_chain_2e(bob_pk.final_A, S, E_A, ProjectivePoint(), ProjectivePoint(), ProjectivePoint());

    return calc_j_invariant(ans.final_A);
}

Fp2 bob_shared(uint64_t m_B, const PublicKey& alice_pk) {

    MontgomeryCurve curve(alice_pk.final_A);
    
    ProjectivePoint S = curve.LADDER3PT(alice_pk.phi_P, alice_pk.phi_Q, alice_pk.phi_R, m_B);

    PublicKey ans = iso_chain_3e(alice_pk.final_A, S, E_B, ProjectivePoint(), ProjectivePoint(), ProjectivePoint());

    return calc_j_invariant(ans.final_A);
}