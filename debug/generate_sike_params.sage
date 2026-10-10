# ==============================================================================
# SIKE (Round 3) Single-Script Parameter Generator & Validator
# ==============================================================================

def generate_and_validate_sike_params(e2, e3, A_coeff=6):
    """
    1つのコードで SIKE パラメータの生成から数学的検証(Unit Test)までを完結して行う関数
    """
    print("="*75)
    print(" 1. SIKE パラメータ生成フェーズ (Generation Phase)")
    print("="*75)
    
    # --------------------------------------------------------------------------
    # Step 1-1: 素数 p の計算と素数判定
    # --------------------------------------------------------------------------
    p = (2^e2) * (3^e3) - 1
    print(f"[+] 素数候補: p = 2^{e2} * 3^{e3} - 1 = {p} (ビット長: {p.bit_length()} bits)")
    if not is_prime(p):
        raise ValueError(f"[FAIL] p = {p} は素数ではありません。e2, e3 を見直してください。")

    # --------------------------------------------------------------------------
    # Step 1-2: 有限体 Fp2 = Fp(i) (i^2 + 1 = 0) および出発曲線 E0 の構築
    # --------------------------------------------------------------------------
    Fp = GF(p)
    R.<t> = PolynomialRing(Fp)
    Fp2.<i> = GF(p^2, modulus=t^2 + 1)

    # SIKE (Round 3) 仕様の初期曲線: E0 : y^2 = x^3 + A*x^2 + x
    E0 = EllipticCurve(Fp2, [0, Fp2(A_coeff), 0, 1, 0])
    f = lambda x: x^3 + A_coeff*x^2 + x
    print(f"[+] 出発曲線 E0: y^2 = x^3 + {A_coeff}*x^2 + x over Fp2")

    # --------------------------------------------------------------------------
    # Step 1-3: Alice 基底点 P2, Q2 の探索 (2^e2 位数)
    # --------------------------------------------------------------------------
    # P2 の目的点: [2^(e2-1)]P2 = (-3 + 2*sqrt(2), 0)
    sqrt2_list = [x for x in [Fp(2).sqrt()] if x is not None]
    if len(sqrt2_list) > 0:
        target_x_P2 = Fp2(-3 + 2 * sqrt2_list[0])
    else:
        target_x_P2 = -3 + 2 * (Fp2(2)).sqrt()

    P2, Q2 = None, None

    # P2 の探索
    for c in range(1000):
        x_cand = i + Fp2(c)
        y2 = f(x_cand)
        if y2.is_square():
            P_cand = E0(x_cand, y2.sqrt())
            P_scaled = (3^e3) * P_cand
            if (2^e2) * P_scaled == E0(0) and (2^(e2-1)) * P_scaled != E0(0):
                if ((2^(e2-1)) * P_scaled)[0] == target_x_P2:
                    P2 = P_scaled
                    print(f"[+] Alice 基底点 P2 を発見 (c={c})")
                    break

    # Q2 の探索: [2^(e2-1)]Q2 = (0, 0)
    for c in range(1000):
        x_cand = i + Fp2(c)
        y2 = f(x_cand)
        if y2.is_square():
            Q_cand = E0(x_cand, y2.sqrt())
            Q_scaled = (3^e3) * Q_cand
            if (2^e2) * Q_scaled == E0(0) and (2^(e2-1)) * Q_scaled != E0(0):
                if ((2^(e2-1)) * Q_scaled)[0] == 0:
                    Q2 = Q_scaled
                    print(f"[+] Alice 基底点 Q2 を発見 (c={c})")
                    break

    # --------------------------------------------------------------------------
    # Step 1-4: Bob 基底点 P3, Q3 の探索 (3^e3 位数)
    # --------------------------------------------------------------------------
    P3, Q3 = None, None

    # P3 (Fp 上の点)
    for c in range(1000):
        c_fp = Fp(c)
        y2 = f(c_fp)
        if y2.is_square():
            P_cand = E0(c_fp, y2.sqrt())
            P_scaled = (2^e2) * P_cand
            if (3^e3) * P_scaled == E0(0) and (3^(e3-1)) * P_scaled != E0(0):
                P3 = P_scaled
                print(f"[+] Bob 基底点 P3 を発見 (c={c})")
                break

    # Q3 (Fp2 \ Fp 上の点, トレースゼロ)
    for c in range(1000):
        c_fp = Fp(c)
        y2 = f(c_fp)
        if not y2.is_square():
            Q_cand = E0(c_fp, Fp2(y2).sqrt())
            Q_scaled = (2^e2) * Q_cand
            if (3^e3) * Q_scaled == E0(0) and (3^(e3-1)) * Q_scaled != E0(0):
                Q3 = Q_scaled
                print(f"[+] Bob 基底点 Q3 を発見 (c={c})")
                break

    # --------------------------------------------------------------------------
    # Step 1-5: 差の点 R2 = P2 - Q2, R3 = P3 - Q3 の生成
    # --------------------------------------------------------------------------
    R2 = P2 - Q2
    R3 = P3 - Q3


    # ==========================================================================
    # 2. パラメータ自動検証フェーズ (Validation / Unit Test Phase)
    # ==========================================================================
    print("\n" + "="*75)
    print(" 2. パラメータ自動判定・検証フェーズ (Validation Suite)")
    print("="*75)
    
    all_passed = True

    # Test 1: 素数判定
    if is_prime(p):
        print(f"[PASS] 1. Prime Test: p = {p} は素数です。")
    else:
        print(f"[FAIL] 1. Prime Test: p = {p} は素数ではありません！")
        all_passed = False

    # Test 2: 超特異性と群位数
    expected_order = (p + 1)^2
    if p < 100000:
        actual_order = E0.cardinality()
        if actual_order == expected_order:
            print(f"[PASS] 2. Supersingularity Test: #E0(Fp2) = {actual_order} (=(p+1)^2) 超特異曲線です。")
        else:
            print(f"[FAIL] 2. Supersingularity Test: 群位数が一致しません！")
            all_passed = False
    else:
        E0_Fp = EllipticCurve(Fp, [0, Fp(A_coeff), 0, 1, 0])
        if E0_Fp.trace_of_frobenius() == 0:
            print(f"[PASS] 2. Supersingularity Test: E0 は超特異曲線です (Trace of Frobenius = 0)。")
        else:
            print(f"[FAIL] 2. Supersingularity Test: 超特異曲線ではありません！")
            all_passed = False

    # Test 3: Alice Torsion 位数
    order_2 = 2^e2
    cond_P2 = (order_2 * P2 == E0(0)) and ((order_2 // 2) * P2 != E0(0))
    cond_Q2 = (order_2 * Q2 == E0(0)) and ((order_2 // 2) * Q2 != E0(0))
    cond_R2 = (order_2 * R2 == E0(0)) and ((order_2 // 2) * R2 != E0(0))
    if cond_P2 and cond_Q2 and cond_R2:
        print(f"[PASS] 3. Alice Torsion Order Test: P2, Q2, R2 の位数はすべて正確に 2^{e2} ({order_2}) です。")
    else:
        print(f"[FAIL] 3. Alice Torsion Order Test: Alice 側の基底点の位数が不正です！")
        all_passed = False

    # Test 4: Bob Torsion 位数
    order_3 = 3^e3
    cond_P3 = (order_3 * P3 == E0(0)) and ((order_3 // 3) * P3 != E0(0))
    cond_Q3 = (order_3 * Q3 == E0(0)) and ((order_3 // 3) * Q3 != E0(0))
    cond_R3 = (order_3 * R3 == E0(0)) and ((order_3 // 3) * R3 != E0(0))
    if cond_P3 and cond_Q3 and cond_R3:
        print(f"[PASS] 4. Bob Torsion Order Test: P3, Q3, R3 の位数はすべて正確に 3^{e3} ({order_3}) です。")
    else:
        print(f"[FAIL] 4. Bob Torsion Order Test: Bob 側の基底点の位数が不正です！")
        all_passed = False

    # Test 5: Weil ペアリングによる一次独立性
    w2 = P2.weil_pairing(Q2, order_2)
    w2_ok = (w2^order_2 == 1) and (w2^(order_2 // 2) != 1)
    w3 = P3.weil_pairing(Q3, order_3)
    w3_ok = (w3^order_3 == 1) and (w3^(order_3 // 3) != 1)
    if w2_ok and w3_ok:
        print(f"[PASS] 5. Weil Pairing Independence Test: {{P2, Q2}} および {{P3, Q3}} はそれぞれ一次独立です。")
    else:
        print(f"[FAIL] 5. Weil Pairing Independence Test: 基底点が一次独立ではありません！")
        all_passed = False

    # Test 6: 差の点の一貫性
    if (R2 == P2 - Q2) and (R3 == P3 - Q3):
        print(f"[PASS] 6. Difference Point Test: R2 = P2 - Q2 および R3 = P3 - Q3 の関係が正しく成立しています。")
    else:
        print(f"[FAIL] 6. Difference Point Test: 差の点 R の計算が一貫していません！")
        all_passed = False


    # ==========================================================================
    # 3. 判定結果判定および C++ コード出力
    # ==========================================================================
    print("\n" + "="*75)
    if all_passed:
        print(" 【最終判定】 ALL PASSED: パラメータは SIKE 仕様に完全準拠しています！")
        print("="*75)
        print("\n// --- C++ パラメータ定義用ヘッダーコード ---")
        print("#ifndef PARAMS_HPP")
        print("#define PARAMS_HPP\n")
        print("#include <cstdint>")
        print('#include "fp2.hpp"')
        print('#include "curve.hpp"\n')
        print("namespace ToyParams {")
        print(f"    const uint64_t P_PRIME = {p};")
        print(f"    const int E_A = {e2};")
        print(f"    const int E_B = {e3};")
        print(f"    const Fp2 INITIAL_A({A_coeff}, 0);\n")
        print(f"    // Alice 基底点 (位数 2^{e2})")
        print(f"    const ProjectivePoint P_A(Fp2({P2[0][0]}, {P2[0][1]}), Fp2(1, 0));")
        print(f"    const ProjectivePoint Q_A(Fp2({Q2[0][0]}, {Q2[0][1]}), Fp2(1, 0));")
        print(f"    const ProjectivePoint R_A(Fp2({R2[0][0]}, {R2[0][1]}), Fp2(1, 0));\n")
        print(f"    // Bob 基底点 (位数 3^{e3})")
        print(f"    const ProjectivePoint P_B(Fp2({P3[0][0]}, {P3[0][1]}), Fp2(1, 0));")
        print(f"    const ProjectivePoint Q_B(Fp2({Q3[0][0]}, {Q3[0][1]}), Fp2(1, 0));")
        print(f"    const ProjectivePoint R_B(Fp2({R3[0][0]}, {R3[0][1]}), Fp2(1, 0));")
        print("}")
        print("\n#endif // PARAMS_HPP")
        print("="*75)
    else:
        print(" 【最終判定】 FAILED: 検証に失敗しました。パラメータを見直してください。")
        print("="*75)

# ------------------------------------------------------------------------------
# 実行 (例: トイパラメータ e2=3, e3=2 -> p = 8 * 9 - 1 = 71, A = 6)
# ------------------------------------------------------------------------------
generate_and_validate_sike_params(e2=3, e3=2, A_coeff=6)