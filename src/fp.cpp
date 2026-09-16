#include "fp.hpp"

// デフォルトコンストラクタ
Fp::Fp() : value(0) {}

// 値を指定するコンストラクタ
Fp::Fp(uint64_t v) : value(v % MODULUS) {}

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
