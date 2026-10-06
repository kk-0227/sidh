#include "fp.hpp"

// デフォルトコンストラクタ
Fp::Fp() : value(0) {}

// 値を指定するコンストラクタ
Fp::Fp(int64_t v) {
    int64_t rem = v % static_cast<int64_t>(MODULUS);
    if (rem < 0) {
        rem += MODULUS;
    }
    value = static_cast<uint64_t>(rem);
}

uint64_t Fp::get_value() const {
    return value;
}

// 加算
Fp Fp::add(const Fp& rhs) const {
    return Fp((value + rhs.value) % MODULUS);
}

// 減算
Fp Fp::sub(const Fp& rhs) const {
    return Fp((value + MODULUS - rhs.value) % MODULUS);
}

// 乗算
Fp Fp::mul(const Fp& rhs) const {
    return Fp((value * rhs.value) % MODULUS);
}

// 繰り返し二乗法
Fp Fp::binpow(uint64_t exp) const {
    Fp ans(1);
    // thisは自身を指すポインタなので*thisで自分自身が取得できる
    Fp base = *this;
    while(exp){
        if(exp&1)ans = ans.mul(base);
        base = base.mul(base);
        exp >>= 1;
    }
    return ans;
}

// フェルマーの小定理
Fp Fp::inv() const {
    return binpow(MODULUS-2);
}

// 除算
Fp Fp::div(const Fp& rhs) const {
    return mul(rhs.inv());
}

void Fp::print() const {
    std::cout << value << std::endl;
}
