# ==============================================================================
# SIDH Toy Parameter Generator (p = 71) - Stable & Correct Version
# ==============================================================================

eA = 3
eB = 2
f = 1
p = f * (2**eA) * (3**eB) - 1 # 1 * 8 * 9 - 1 = 71

print(f"/* ===== Generating Parameters for p = {p} ===== */")

# 拡大体 F_{p^2} の定義
Fp = GF(p)
R = PolynomialRing(Fp, 'x')
x = R.gen()
Fp2 = GF(p**2, name='i', modulus=x**2 + 1)
i = Fp2.gen()

# 初期モンゴメリ曲線 E_0: y^2 = x^3 + 6x^2 + x
A0 = Fp2(6)
E0 = EllipticCurve(Fp2, [0, A0, 0, 1, 0])

# 超特異楕円曲線の各生成元の最大位数は p+1
max_order = p + 1

# 群の生成元を2つ取得
G = E0.abelian_group()
gens = G.gens()
G1 = gens[0].element()
G2 = gens[1].element()

# Alice の基底点 (位数 2^eA)
cofactor_A = max_order // (2**eA) # 72 // 8 = 9
PA = cofactor_A * G1
QA = cofactor_A * G2
RA = PA - QA

# Bob の基底点 (位数 3^eB)
cofactor_B = max_order // (3**eB) # 72 // 9 = 8
PB = cofactor_B * G1
QB = cofactor_B * G2
RB = PB - QB

# C++ 定数出力関数
def print_cpp_point(name, Pt):
    if Pt.is_zero():
        print(f"// ERROR: {name} is the Point at Infinity")
        return
    
    x_val = Pt.x()
    # Fp2要素から実部と虚部を安全に抽出
    coeffs = x_val.polynomial().list()
    real_part = coeffs[0] if len(coeffs) > 0 else 0
    imag_part = coeffs[1] if len(coeffs) > 1 else 0
    
    print(f"const ProjectivePoint {name}(Fp2({real_part}, {imag_part}), Fp2(1, 0));")

print("// ===== C++ Global Constants (Copy & Paste to your code) =====")
print(f"const uint64_t P_PRIME = {p};")
print(f"const int E_A = {eA};")
print(f"const int E_B = {eB};")
print(f"const Fp2 INITIAL_A(6, 0);\n")

print("// Alice's Basis Points (Order 2^eA)")
print_cpp_point("P_A", PA)
print_cpp_point("Q_A", QA)
print_cpp_point("R_A", RA)

print("\n// Bob's Basis Points (Order 3^eB)")
print_cpp_point("P_B", PB)
print_cpp_point("Q_B", QB)
print_cpp_point("R_B", RB)