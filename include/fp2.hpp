#ifndef FP2_HPP
#define FP2_HPP

#include "fp.hpp"
#include <iostream>

// Fp2 の要素: a + b*i (a, b ∈ Fp, i^2 = -1)
class Fp2 {
private:
    Fp a; // 実部
    Fp b; // 虚部

public:
    // コンストラクタ、引数の渡し方によって初期化方法を変える
    Fp2() : a(0), b(0) {}
    Fp2(const Fp& real, const Fp& imag) : a(real), b(imag) {}
    Fp2(uint64_t real, uint64_t imag = 0) : a(real), b(imag) {}

    Fp get_real() const { return a; }
    Fp get_imag() const { return b; }

    // 加算: (a + bi) + (c + di) = (a + c) + (b + d)i
    Fp2 add(const Fp2& rhs) const {
        return Fp2(a + rhs.a, b + rhs.b);
    }

    // 減算: (a + bi) - (c + di) = (a - c) + (b - d)i
    Fp2 sub(const Fp2& rhs) const {
        return Fp2(a - rhs.a, b - rhs.b);
    }

    // 乗算: (a + bi)(c + di) = (ac - bd) + (ad + bc)i
    Fp2 mul(const Fp2& rhs) const{ return Fp2(a*rhs.a-b*rhs.b, a*rhs.b+b*rhs.a); }

    // 逆元: (a + bi)^(-1) = (a - bi) / (a^2 + b^2)
    Fp2 inv() const{
        Fp norm_inv = (a*a + b*b).inv();
        return Fp2(a*norm_inv, -b*norm_inv);
    }

    // 除算
    Fp2 div(const Fp2& rhs) const { return mul(rhs.inv()); }

    void print() const{
        std::cout << a.get_value() << "+" << b.get_value() << "i" << std::endl;
    }
    Fp2 operator+(const Fp2& rhs) const { return add(rhs); }
    Fp2 operator-(const Fp2& rhs) const { return sub(rhs); }
    Fp2 operator*(const Fp2& rhs) const { return mul(rhs); }
    Fp2 operator/(const Fp2& rhs) const { return div(rhs); }

    Fp2 operator+=(const Fp2& rhs) { *this = add(rhs); return *this; }
    Fp2 operator-=(const Fp2& rhs) { *this = sub(rhs); return *this; }
    Fp2 operator*=(const Fp2& rhs) { *this = mul(rhs); return *this; }
    Fp2 operator/=(const Fp2& rhs) { *this = div(rhs); return *this; }

    Fp2 operator-() const {return Fp2(-a, -b); }
    // operator==の中ではFp2の==まだ定義されていない
    bool operator==(const Fp2& rhs) const {return (a == rhs.a) && (b == rhs.b); }
    bool operator!=(const Fp2& rhs) const {return *this != rhs; }
};

#endif // FP2_HPP