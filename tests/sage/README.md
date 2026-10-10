# Sage による検証

SageMath 組み込みの同種写像(`EllipticCurve.isogeny`)で、鍵交換の共有秘密(j 不変量)を
独立に計算し、C++ 実装の出力と比較します。

## 実行方法
```bash
sage tests/sage/verify_key_exchange.sage
```
`mA mB j` の 72 行が出力されます(`mA` は 0〜7、`mB` は 0〜8)。
C++ 側の共有秘密(`alice_shared` の出力)と、全 72 通りで一致することを確認しています。