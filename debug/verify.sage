import re

def parse_fp2_str(s, p):
    """ '9+0i' のような文字列を SageMath の Fp2 エレメントに変換 """
    s = s.strip()
    match = re.match(r"([+-]?\d+)\+([+-]?\d+)i", s)
    if match:
        re_val = int(match.group(1)) % p
        im_val = int(match.group(2)) % p
        return F2(re_val + im_val * i)
    return F2(int(s) % p)

# --- 1. パラメータ初期化 ---
p = 71
F2.<i> = GF(p^2, modulus=x^2+1)

print(f"=== SIDH Verification (p = {p}) ===")

# --- 2. ログファイルの読み込み ---
log_path = "debug/debug_log.txt"

try:
    with open(log_path, "r") as f:
        lines = f.readlines()
except FileNotFoundError:
    print(f"Error: {log_path} が見つかりません。")
    exit(1)

log_data = {}
current_step = None

for line in lines:
    line = line.strip()
    if line.startswith("[") and line.endswith("]"):
        current_step = line[1:-1]
        log_data[current_step] = {}
    elif "=" in line and current_step:
        key, val = line.split("=")
        log_data[current_step][key.strip()] = parse_fp2_str(val, p)

# --- 3. SageMathでの曲線・基準点構築 ---
A_cpp = log_data["Step0_Init"]["A"]
X0_cpp = log_data["Step0_Init"]["X"]
Z0_cpp = log_data["Step0_Init"]["Z"]

E = EllipticCurve(F2, [0, A_cpp, 0, 1, 0])

x0_affine = X0_cpp / Z0_cpp
try:
    P = E.lift_x(x0_affine)
except ValueError:
    print(f"\n❌ エラー: x = {x0_affine} に対応する点が曲線上に存在しません。")
    exit(1)


# --- [Step 1: xDBL 検証 (Q = 2 * P)] ---
if "Step1_After_xDBL" in log_data:
    P_dbl_sage = P * 2
    
    X1_cpp = log_data["Step1_After_xDBL"]["X"]
    Z1_cpp = log_data["Step1_After_xDBL"]["Z"]
    x1_cpp_affine = X1_cpp / Z1_cpp
    x1_sage_affine = P_dbl_sage[0]

    print("\n[Step1_After_xDBL (2P)]")
    print(f"  C++ 側 x 座標     : {x1_cpp_affine}")
    print(f"  SageMath 側 x 座標: {x1_sage_affine}")

    assert x1_cpp_affine == x1_sage_affine, "❌ Step1 (xDBL) で不一致が発生しました！"
    print("  ✅ SUCCESS: 完全一致")


# --- [Step 2: xADD 検証 (R = P + Q = 3 * P)] ---
if "Step2_After_xADD" in log_data:
    # 2P (Step1の結果) に P を足す (3P)
    P_add_sage = P * 3
    
    X2_cpp = log_data["Step2_After_xADD"]["X"]
    Z2_cpp = log_data["Step2_After_xADD"]["Z"]
    x2_cpp_affine = X2_cpp / Z2_cpp
    x2_sage_affine = P_add_sage[0]

    print("\n[Step2_After_xADD (3P)]")
    print(f"  C++ 側 x 座標     : {x2_cpp_affine}")
    print(f"  SageMath 側 x 座標: {x2_sage_affine}")

    assert x2_cpp_affine == x2_sage_affine, "❌ Step2 (xADD) で不一致が発生しました！"
    print("  ✅ SUCCESS: 完全一致")


# --- [Step 3: xMUL 検証 (S = 5 * P)] ---
if "Step3_After_xMUL" in log_data:
    # 5 * P
    P_mul_sage = P * 5
    
    X3_cpp = log_data["Step3_After_xMUL"]["X"]
    Z3_cpp = log_data["Step3_After_xMUL"]["Z"]
    x3_cpp_affine = X3_cpp / Z3_cpp
    x3_sage_affine = P_mul_sage[0]

    print("\n[Step3_After_xMUL (5P)]")
    print(f"  C++ 側 x 座標     : {x3_cpp_affine}")
    print(f"  SageMath 側 x 座標: {x3_sage_affine}")

    assert x3_cpp_affine == x3_sage_affine, "❌ Step3 (xMUL) で不一致が発生しました！"
    print("  ✅ SUCCESS: 完全一致")

print("\n全ステップの検証が正常に完了しました。")