#ifndef ISOGENY_HPP
#define ISOGENY_HPP

#include "fp2.hpp"
#include "curve.hpp"

// 2-同種写像の計算を行うクラス
class Isogeny2 {
private:
    ProjectivePoint K; // 核(Kernel)となる位数2の点 (XK : ZK)

public:
    // コンストラクタ: 核となる点 K を渡して初期化
    explicit Isogeny2(const ProjectivePoint& K);

    // 新しい曲線 E' のパラメータ A' を計算する (BuildIsogeny)
    Fp2 get_next_A() const;

    // 点 Q を新しい曲線 E' 上の点 phi(Q) に写像する (EvalIsogeny)
    ProjectivePoint eval(const ProjectivePoint& Q) const;
};


// 3-同種写像の計算を行うクラス
class Isogeny3 {
private:
    ProjectivePoint K; // 核(Kernel)となる位数3の点 (XK : ZK)

public:
    // コンストラクタ: 核となる点 K を渡して初期化
    explicit Isogeny3(const ProjectivePoint& K);

    // 新しい曲線 E' のパラメータ A' を計算する (BuildIsogeny)
    Fp2 get_next_A(const Fp2& current_A) const;

    // 点 Q を新しい曲線 E' 上の点 phi(Q) に写像する (EvalIsogeny)
    ProjectivePoint eval(const ProjectivePoint& Q) const;
};

struct IsogenyChainResult {
    Fp2 final_A;              // 最終的に到達した曲線のパラメータ A'
    ProjectivePoint phi_P;    // 写像された点 phi(P)
    ProjectivePoint phi_Q;    // 写像された点 phi(Q)
    ProjectivePoint phi_R;    // phi(R) ここでR=P-Qでこれは楕円曲線上の加算P+Qで必要
};

// 2^e同種写像
// 最初のラウンドではstart_Aは共通の楕円曲線E_0,Sは位数2^eの点(ここから位数2の点をe個取り出し、veluの公式で使う)
// eは同種写像のステップ数(最初に指定される),P,Qはボブの基底点でこれをアリスの写像に通したあと計算用のRと一緒にボブに渡す
// 二度目のラウンドではstart_Aはボブから送られてきた楕円曲線,Sはボブから来たP,Qで作った位数2^eの点,PQRは使わない
IsogenyChainResult iso_chain_2e(
    const Fp2& start_A, 
    const ProjectivePoint S,
    int e, 
    const ProjectivePoint P, 
    const ProjectivePoint Q, 
    const ProjectivePoint R
);

// 3^e同種写像
IsogenyChainResult iso_chain_3e(
    const Fp2& start_A, 
    const ProjectivePoint S,
    int e, 
    const ProjectivePoint P, 
    const ProjectivePoint Q, 
    const ProjectivePoint R
);

#endif // ISOGENY_HPP