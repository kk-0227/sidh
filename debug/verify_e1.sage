# ==========================================
# verify_e1.sage
# ==========================================
eA = 3
eB = 2
p = 1 * (2**eA) * (3**eB) - 1

Fp = GF(p)
R = PolynomialRing(Fp, 'x')
x = R.gen()
Fp2 = GF(p**2, name='i', modulus=x**2 + 1)
i = Fp2.gen()

A0 = Fp2(6)
E0 = EllipticCurve(Fp2, [0, A0, 0, 1, 0])

# generate_params と同じ方法で完全に同一の点を生成
#max_order = p + 1
#G = E0.abelian_group()
#gens = G.gens()
#G1 = gens[0].element()
#G2 = gens[1].element()

#cofactor_A = max_order // (2**eA)
#PA = cofactor_A * G1
#QA = cofactor_A * G2

#cofactor_B = max_order // (3**eB)
#PB = cofactor_B * G1
#QB = cofactor_B * G2
#RB = PB - QB
PA = E0(
    Fp2(69 + 22*i),
    Fp2(70 + 49*i)
)

QA = E0(
    Fp2(25 + 21*i),
    Fp2(56 + 28*i)
)

RA = E0(
    Fp2(48 + 59*i),
    Fp2(63 + 24*i)
)

PB = E0(
    Fp2(4 + 40*i),
    Fp2(45 + 42*i)
)

QB = E0(
    Fp2(12 + 58*i),
    Fp2(68 + 32*i)
)

RB = E0(
    Fp2(0 + 39*i),
    Fp2(34 + 22*i)
)

# ==========================================
# 1ステップの同種写像の検証 (e = 1)
# ==========================================

# 1. 秘密鍵 m_A = 1 の秘密点 S_A = PA + QA
SA = PA + QA

# 2. 2倍算を2回(4倍)して、位数 2 の核 K を抽出
K = 4 * SA
SA = PA + QA
K = 4 * SA

print("S_A =", SA.x())
print("K =", K.x())

print("=== SageMath Verification (e = 1) ===")
print(f"Kernel K (Affine x): {K.x()}")

# 3. K を核とする 2-同種写像を計算
phi = E0.isogeny(K)
E1 = phi.codomain()

print("\n[2] Next Curve Parameter (SageMath):")
print(f"A' (a2 parameter) : {E1.a2()}")

print("\n[3] Pushed-forward Points (Bob's Basis in SageMath):")
print(f"phi(P_B) (Affine x): {phi(PB).x()}")
print(f"phi(Q_B) (Affine x): {phi(QB).x()}")
print(f"phi(R_B) (Affine x): {phi(RB).x()}")