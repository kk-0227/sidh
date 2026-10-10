# ==========================================
# verify_e1_affine.sage
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
max_order = p + 1
G = E0.abelian_group()
gens = G.gens()
G1 = gens[0].element()
G2 = gens[1].element()

cofactor_A = max_order // (2**eA)
PA = cofactor_A * G1
QA = cofactor_A * G2

cofactor_B = max_order // (3**eB)
PB = cofactor_B * G1
QB = cofactor_B * G2
RB = PB - QB

# ==========================================
# 1ステップの同種写像の検証 (e = 1) - アフィン完全版
# ==========================================

# 1. 秘密鍵 m_A = 1 の秘密点 S_A = PA + QA
SA = PA + QA

# 2. 2倍算を2回(4倍)して、位数 2 の核 K を抽出
K = 4 * SA

print("=== SageMath Verification (e = 1) [Affine Coordinates] ===")

# アフィン座標表示用ヘルパー関数
def format_affine_point(Pt):
    if Pt.is_zero():
        return "Point at Infinity"
    
    x_val = Pt.x()
    x_coeffs = x_val.polynomial().list()
    xr = x_coeffs[0] if len(x_coeffs) > 0 else 0
    xi = x_coeffs[1] if len(x_coeffs) > 1 else 0

    y_val = Pt.y()
    y_coeffs = y_val.polynomial().list()
    yr = y_coeffs[0] if len(y_coeffs) > 0 else 0
    yi = y_coeffs[1] if len(y_coeffs) > 1 else 0

    return f"AffinePoint(Fp2({xr}, {xi}), Fp2({yr}, {yi}))"

print(f"Kernel K : {format_affine_point(K)}")

# 3. K を核とする 2-同種写像を計算
phi = E0.isogeny(K)
E1 = phi.codomain()

print("\n[2] Next Curve Parameter (SageMath):")
print(f"A' (a2 parameter) : {E1.a2()}")

print("\n[3] Pushed-forward Points (Bob's Basis in SageMath):")
print(f"phi(P_B) : {format_affine_point(phi(PB))}")
print(f"phi(Q_B) : {format_affine_point(phi(QB))}")
print(f"phi(R_B) : {format_affine_point(phi(RB))}")