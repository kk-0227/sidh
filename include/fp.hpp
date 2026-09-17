#ifndef FP_HPP
#define FP_HPP

#include <iostream>
#include <array>
#include <cstdint>

// SIDHで用いる素体 F_p の要素を扱うクラス
class Fp {
private:
    uint64_t value;

    // SIDHで使用する素数 p (例として小さな素数、またはパラメータを定義)
    static constexpr uint64_t MODULUS = 97; // 初期テスト用の小さな素数

public:
    // コンストラクタ
    Fp();
    Fp(uint64_t v);

    // 値の取得
    uint64_t get_value() const;

    // 有限体 F_p 上の四則演算（メソッド）
    Fp add(const Fp& rhs) const;
    Fp sub(const Fp& rhs) const;
    Fp mul(const Fp& rhs) const;
    Fp div(const Fp& rhs) const;
    Fp binpow(uint64_t exp) const;
    Fp inv() const;

    // 演算子オーバーロード
    Fp operator+(const Fp& rhs) const { return add(rhs); }
    Fp operator-(const Fp& rhs) const { return sub(rhs); }
    Fp operator*(const Fp& rhs) const { return mul(rhs); }
    Fp operator/(const Fp& rhs) const { return div(rhs); }
    // +=はFpのvalueを和に書き換えたものを返す
    Fp operator+=(const Fp& rhs)  {
        value = add(rhs).value;
        return *this;
    }
    Fp operator-=(const Fp& rhs)  {
        value = sub(rhs).value;
        return *this;
    }
    Fp operator*=(const Fp& rhs)  {
        value = mul(rhs).value;
        return *this;
    }
    Fp operator/=(const Fp& rhs)  {
        value = div(rhs).value;
        return *this;
    }
    Fp operator-() const { return Fp(0) - *this; }
    bool operator==(const Fp& rhs) const {
        return value == rhs.value;
    }
    bool operator!=(const Fp& rhs) const {
        return value != rhs.value;
    }
    // デバッグ出力用
    void print() const;
};

#endif // FP_HPP