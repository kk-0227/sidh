# SIDH

耐量子計算機暗号として提唱された暗号の一つである、SIDH(Supersingular Isogeny Diffie-Hellman)のC++実装です。

## 目的
独学で学んだ SIDH(同種写像を用いた鍵交換)を、実際に実装して動かす経験をすること、
そして、理論と実装の違いやギャップを理解することが目的です。
実装を通して、理論の理解を深めることにもつながると考えています。

## できること
Alice の秘密鍵 `m_A` と Bob の秘密鍵 `m_B` を用いて、二人だけが共有できる秘密の数値(j-invariant)を生成できます。

## 使用技術
- 言語: C++17(標準ライブラリのみ)
- 検証: SageMath(同種写像の照合)
- バージョン管理: Git/GitHub

## 工夫した点
- 機能ごとに段階的に実装し、都度テストを行い正当性を担保しながら進めました。(有限体→楕円曲線上の演算→同種写像→連鎖→鍵交換)
- 理論の理解を深め可読性の高いコードにするため、リファレンスでの速さや安全性を優先した書き方ではなく、演算子オーバーロードや値を管理する構造体を用いた書き方で実装しました。

## ビルドと実行
```bash
git clone https://github.com/kk-0227/sidh.git
cd sidh

# ビルドして、テストを実行
./compile.sh

# 鍵交換のデモ
./build/main
```
C++17 対応の `g++` が必要です。

## 検証
- **単体テスト**(`tests/test_all.cpp`): 曲線演算、同種写像、鍵交換を、性質(準同型性など)で確認します。
- **SageMath との照合**(`tests/sage/`): Sage 組み込みの同種写像で、共有秘密を異なる方法で計算し、
  C++ の出力(`./build/dump_shared_secrets`)と、秘密鍵のペア全 72 通りで一致することを確認しています。
  詳しくは [tests/sage/README.md](tests/sage/README.md) を参照してください。

## 今後の課題
- 本番パラメータ(P434 など)への対応(多倍長整数が必要)。

## 注意
- 学習目的です。実運用には使わないでください。
- トイパラメータのため、安全性はまったくありません。
- 定数時間実装などのサイドチャネル対策はしていません。
- SIDH/SIKEは2022年に[効率的な鍵回復攻撃](https://eprint.iacr.org/2022/975)が発表され、現在は安全でないと分かっています。

## 参考
- C. Costello, P. Longa, M. Naehrig,
  [Efficient algorithms for supersingular isogeny Diffie-Hellman](https://eprint.iacr.org/2016/413)
  (CRYPTO 2016)
- [microsoft/PQCrypto-SIDH](https://github.com/microsoft/PQCrypto-SIDH)
