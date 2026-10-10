# ==========================================
# verify_e3.sage
# ==========================================
p = 71
Fp = GF(p)
R = PolynomialRing(Fp, 'x')
x = R.gen()
Fp2 = GF(p**2, name='i', modulus=x**2 + 1)
i = Fp2.gen()

E0 = EllipticCurve(Fp2, [0, 6, 0, 1, 0])

# --- 基底点の生成 ---
max_order = p + 1
G = E0.abelian_group().gens()
G1, G2 = G[0].element(), G[1].element()
PA, QA = (max_order // 8) * G1, (max_order // 8) * G2
PB, QB = (max_order // 9) * G1, (max_order // 9) * G2
RB = PB - QB

# --- SIKE Montgomery Formulas in Python ---
def get_next_A(XK, ZK):
    C24 = (XK**2 - ZK**2)**2
    A24plus = 2 * (XK**2 + ZK**2)**2
    return (4 * A24plus - 2 * C24) / C24

def isogeny_eval(XQ, ZQ, XK, ZK):
    t0 = (XQ + ZQ) * (XK - ZK)
    t1 = (XQ - ZQ) * (XK + ZK)
    t2 = t0 + t1
    t3 = t0 - t1
    t2 = t2**2
    t3 = t3**2
    next_X = XQ * t2
    next_Z = ZQ * t3
    return next_X, next_Z

def xDBL(X, Z, A):
    t0 = (X + Z)**2
    t1 = (X - Z)**2
    C = t0 - t1
    next_X = t0 * t1
    A24 = (A + Fp2(2)) / Fp2(4)
    t0 = t0 + A24 * C
    next_Z = C * t0
    return next_X, next_Z

# --- ループの初期設定 ---
S = PA + QA
SX, SZ = S.x(), Fp2(1)

PX, PZ = PB.x(), Fp2(1)
QX, QZ = QB.x(), Fp2(1)
RX, RZ = RB.x(), Fp2(1)
A = Fp2(6)

print("=== SIKE e=3 Chain Verification ===")
print(f"[1] Initial Secret Point S (Affine x): {SX}")

# --- e=3 の同種写像チェーンループ ---
for step in range(3):
    # 核 K の計算: 2^(2 - step) 倍
    KX, KZ = SX, SZ
    for _ in range(2 - step):
        KX, KZ = xDBL(KX, KZ, A)
    
    # 次の曲線を計算
    A_next = get_next_A(KX, KZ)
    
    # 点の移送
    if step < 2:
        SX, SZ = isogeny_eval(SX, SZ, KX, KZ)
    PX, PZ = isogeny_eval(PX, PZ, KX, KZ)
    QX, QZ = isogeny_eval(QX, QZ, KX, KZ)
    RX, RZ = isogeny_eval(RX, RZ, KX, KZ)
    
    # 曲線パラメータの更新
    A = A_next

# --- アフィン座標に戻して出力 ---
print(f"\n[2] Expected Next Curve Parameter A': {A}")
print("\n[3] Expected Pushed-forward Points (Affine x):")
print(f"phi(P_B) : {PX / PZ}")
print(f"phi(Q_B) : {QX / QZ}")
print(f"phi(R_B) : {RX / RZ}")