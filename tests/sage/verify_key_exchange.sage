# Sage 組み込みの isogeny を使って共有秘密(j 不変量)を計算し、C++ 実装の結果と照合する。
# 実行: sage tests/sage/verify_key_exchange.sage
# 基底点は params.hpp の x 座標から作り、R = P - Q の関係に符号を合わせている。

p = 71
Fp2.<i> = GF(p^2, modulus=x^2 + 1)
E0 = EllipticCurve(Fp2, [0, 6, 0, 1, 0])


# params.hpp の x 座標から、Sage の点を作る
P_A = E0.lift_x(Fp2(19 + 66*i))
Q_A = E0.lift_x(Fp2(19 + 62*i))
R_A_x = Fp2(62 + 26*i)
if (P_A - Q_A).xy()[0] != R_A_x:    # 符号の向きを、R_A = P_A - Q_A に合わせる
    Q_A = -Q_A
assert (P_A - Q_A).xy()[0] == R_A_x

P_B = E0.lift_x(Fp2(45))
Q_B = E0.lift_x(Fp2(17))
R_B_x = Fp2(52 + 16*i)
if (P_B - Q_B).xy()[0] != R_B_x:
    Q_B = -Q_B
assert (P_B - Q_B).xy()[0] == R_B_x

def fmt(j):                 # C++ の出力 "a+bi" に合わせる
    c = j.polynomial().list() + [0, 0]
    return f"{int(c[0])}+{int(c[1])}i"

for mA in range(8):
    for mB in range(9):
        SA, SB = P_A + mA * Q_A, P_B + mB * Q_B
        phiA = E0.isogeny(SA)
        EA = phiA.codomain()
        EAB = EA.isogeny(phiA(SB)).codomain()       # 共有秘密の曲線
        print(mA, mB, fmt(EAB.j_invariant()))