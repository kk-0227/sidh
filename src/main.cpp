#include <iostream>
#include "key_exchange.hpp"

// 鍵交換のデモ: 秘密鍵 m_A, m_B から、双方が同じ共有秘密を得ることを示す。
int main() {
    const uint64_t m_A = 3;   // Alice の秘密鍵 (0〜7)
    const uint64_t m_B = 5;   // Bob の秘密鍵 (0〜8)

    PublicKey pk_A = alice_keygen(m_A);
    PublicKey pk_B = bob_keygen(m_B);

    Fp2 secret_A = alice_shared(m_A, pk_B);
    Fp2 secret_B = bob_shared(m_B, pk_A);

    std::cout << "Alice の共有秘密: "; secret_A.print();
    std::cout << "Bob の共有秘密:   "; secret_B.print();
    std::cout << (secret_A == secret_B ? "一致しました" : "一致しませんでした") << std::endl;

    return secret_A == secret_B ? 0 : 1;
}