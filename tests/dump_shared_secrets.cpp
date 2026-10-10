// 全秘密鍵の組み合わせについて、共有秘密(j不変量)を "mA mB j" の形式で出力する。
// tests/sage/verify_key_exchange.sage の出力と diff で照合するために使う。
#include <iostream>
#include "key_exchange.hpp"

int main() {
    for (uint64_t m_A = 0; m_A < 8; ++m_A) {
        for (uint64_t m_B = 0; m_B < 9; ++m_B) {
            PublicKey pk_B = bob_keygen(m_B);
            std::cout << m_A << " " << m_B << " ";
            alice_shared(m_A, pk_B).print();
        }
    }
    return 0;
}