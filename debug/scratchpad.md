sageの起動方法
'conda activate sage'

- `Fp` : 有限体 ($\mathbb{F}_p$)
  - `get_value()` : 値の取得

- `Fp2` : 二次の拡大体 ($\mathbb{F}_{p^2} = \mathbb{F}_p[i] / (i^2 + 1)$)
  - `a + bi` の形で表現
  - 初期化
    - `Fp2 P;` : `0+0i`で初期化
    - `Fp2 P(5)` : `5+0i`で初期化
    - `Fp2 P(5,3)` : `5+3i`で初期化
  - `get_real()` : 実部を`Fp`で返す
  - `get_imag()` : 虚部を`Fp`で返す
    - 値が欲しい時は`get_real().get_value()`みたいに書く
  - `is_zero()` : `0+0i`かどうかを`bool`で判定
    - 射影座標のZが`0+0i`ならその点は無限遠点を表す

- `ProjectivePoint` : 射影座標上の点$(X:Z)$
  - 初期化
    - `ProjectivePoint P;` : `X = 0+0i, Z = 1+0i`で初期化
    - `ProjectivePoint P(x,z);` : `X = x, Z = z`で初期化

- `MontgomeryCurve` : モンゴメリー曲線($y^2=x^3+Ax^2+x$)の生成
    - 初期化
      - `MontgomeryCurve curve(A)` : (`A`は`Fp2`)
    - `get_A()` : Aを`Fp2`で返す
    - 楕円曲線上の演算($P,Q$は楕円曲線上の射影座標表示の点)
      - `xDBL(P)` : 2倍算( $2P$ を返す)
      - `xADD(P, Q, P_minus_Q)` : 加算(P+Q)
        - `P_minus_Q`はモンゴメリーラダーでは常に最初の点
      - `xMUL(P,n)` : n倍算(nは検証用の`uint64_t`)
      - `calc_j_invariant(A)` : パラメータが $A$ のj-不変量を返す(型は`Fp2`)

- `Isogey2` : 2-同種写像の計算
  - 初期化
    - `Isogeny2 iso2(K);` : (`K` : 写像の核となる位数2の`ProjectivePoint`)
  - `get_next_A()` : 行き先の楕円曲線 $E'$ の係数 $A'$ を返す
  - `eval(Q)` : 点 $Q$ を新しい曲線上へ写像した点 $\phi(Q)$ (`ProjectivePoint`) を返す

- `Isogeny3` : 3-同種写像の計算クラス
  - コンストラクタ: `Isogeny3 iso3(K);`
    - 引数 `K`: 核となる**位数3の点**
  - `get_next_A(current_A)` : 現在の曲線 $A$ を受け取り、新しい曲線 $A'$ を算出
  - `eval(Q)` : 点 $Q$ を新しい曲線上へ写像した点 $\phi(Q)$ を算出

- `IsogenyChainResult` : 同種写像チェーンの計算結果を保持する構造体
  - メンバ変数:
    - `final_A` (`Fp2`) : 最終的に到達した曲線のパラメータ $A'$
    - `phi_P`, `phi_Q`, `phi_R` (`ProjectivePoint`) : 写像後の各点 $\phi(P), \phi(Q), \phi(R)$
      - $P$,$Q$の二つの点とその差分$R = P - Q$
        - $R$は加算$P + Q$で使用

- `iso_chain_2e(start_A, S, e, P, Q, R)` : $2^e$-同種写像の連鎖計算
  - 引数:
    - `start_A` : 開始曲線のパラメータ $A$
    - `S` : 核を生成するための**位数 $2^e$ の点**
      - これは与えられた位数 $2^e$ の点と自分のパラメータの演算で作ったもの
    - `e` : ステップ数（同種写像を繰り返す回数）
    - `P`, `Q`, `R` : 相手に渡すために写像（移送）させる基底点群（$R = P - Q$）
      -  $P,Q,R$ は位数 $3^e$ の側の基底点であることに注意
  - 内部処理の概要:
    1. ループ内で $S$ から2倍算（`xDBL`）を繰り返し、位数2の点 $K$ を取り出す
    2. `Isogeny2(K)` を構築し、曲線 $A'$ の更新および $S, P, Q, R$ の評価（`eval`）を行う
    3. これを $e$ 回繰り返す
  - 戻り値: 最終曲線 $A'$ と写像された点群を含む `IsogenyChainResult`

- `iso_chain_3e(start_A, S, e, P, Q, R)` : $3^e$-同種写像の連鎖計算