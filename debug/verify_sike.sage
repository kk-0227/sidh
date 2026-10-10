# ==========================================
# verify_sike.sage
# ==========================================
p = 71
Fp = GF(p)
R = PolynomialRing(Fp, 'x')
x = R.gen()
Fp2 = GF(p**2, name='i', modulus=x**2 + 1)
i = Fp2.gen()

E0 = EllipticCurve(Fp2, [0, 6, 0, 1, 0])

# --- Generate Basis ---
max_order = p + 1
G = E0.abelian_group().gens()
G1, G2 = G[0].element(), G[1].element()

PA, QA = (max_order // 8) * G1, (max_order // 8) * G2
PB, QB = (max_order // 9) * G1, (max_order // 9) * G2
RB = PB - QB

# --- Kernel K ---
K = 4 * (PA + QA)
xK = K.x()

print("=== SIKE Specification Verification ===")
print(f"Kernel K (Affine x) : {xK}")

# --- SIKE 2-Isogeny Formulas ---
XK2 = xK^2
ZK2 = Fp2(1)

C24 = (XK2 - ZK2)^2
A24plus = 2 * (XK2 + ZK2)^2
next_A = (4 * A24plus - 2 * C24) / C24

print(f"\n[2] Expected Next A'  : {next_A}")

def sike_eval(Q):
    XQ, ZQ = Q.x(), Fp2(1)
    t0 = (XQ + ZQ) * (xK - Fp2(1))
    t1 = (XQ - ZQ) * (xK + Fp2(1))
    next_X = XQ * (t0 + t1)
    next_Z = ZQ * (t0 - t1)
    return next_X / next_Z

print("\n[3] Expected Pushed-forward Points (Affine x):")
print(f"phi(P_B) : {sike_eval(PB)}")
print(f"phi(Q_B) : {sike_eval(QB)}")
print(f"phi(R_B) : {sike_eval(RB)}")