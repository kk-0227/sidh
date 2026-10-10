#ifndef KEY_EXCHANGE_HPP
#define KEY_EXCHANGE_HPP

#include "isogeny.hpp"

using PublicKey = IsogenyChainResult;

// Alice: 秘密鍵 m_A から公開鍵を作る
PublicKey alice_keygen(uint64_t m_A);

// Bob: 秘密鍵 m_B から公開鍵を作る
PublicKey bob_keygen(uint64_t m_B);

// Alice: 自分の秘密鍵 m_A と、Bob の公開鍵から、共有秘密(j 不変量)を作る
Fp2 alice_shared(uint64_t m_A, const PublicKey& bob_pk);

// Bob: 同様
Fp2 bob_shared(uint64_t m_B, const PublicKey& alice_pk);

#endif // KEY_EXCHANGE_HPP