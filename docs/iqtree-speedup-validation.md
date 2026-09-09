# IQ-TREE／NWKITの追加最適化（2026-09-09）

[性能レビュー](iqtree-performance-review.md)の優先項目1〜3を実装した。

- IQ-TREEの `SCORE` 要求は尤度と全枝の一次微分のみを計算・返却する。二階微分が必要な `EVAL` は維持する。通常の `value_gradient()` は `SCORE` を使い、同じ点で後から対角を要求した場合はキャッシュを一度更新する。
- 固定したモデル・アラインメント・計算表現に対して、末端尤度表とパターン情報を再利用する。枝長に依存する部分尤度は毎要求で無効化し、spectral／original-stateの表現切替時には末端表も作り直す。
- NWKITのmarginalプロファイルで、共分散、観測尤度の逆Hessian、root contrastなどを読み取り専用の共通構造として共有する。年代制約は各点で作り直し、rhoや尤度近似が変わった場合は共通構造も再構築する。求積点も再利用する。

モデル、最適化の開始点数、求積精度、近似の検証条件、95%区間の定義は変更していない。IQ-TREEの計算変更はIQ-TREE、推定処理の再利用はNWKITに置いた。GeneGalleonのcore実装には追加変更を加えていない。

## 測定条件

Apple M2 Max／macOS 26.6.2上のDockerのGeneGalleon runtime、Linux arm64（Docker割当12 CPU・約7.65 GiB）、IQ-TREEは1スレッド、`OPENBLAS_NUM_THREADS=1`、`OMP_NUM_THREADS=1`。比較元は追加最適化直前の常駐実装であり、毎回プロセスを起動する旧方式ではない。

カーネル比較は4配列×150コドン、16／64配列×1,500コドンの平衡木、固定 `GY{0.5,2}+FQ+G4{1}`、各10回の異なる枝長要求、各条件3反復。各試行は新しいコンテナ／Python／IQ-TREE workerで開始し、1回のウォームアップ後に評価時間を測る。条件の実行順を反復ごとに入れ替え、初期化時間と子プロセスの最大RSSも記録する。RSSはPythonとの合計常駐量ではない。

```bash
python3 /Users/kf/repos/nwkit/tools/benchmark_radte_iqtree_compare.py \
  --baseline-image local/genegalleon:iqtree-before-speedup \
  --candidate-image local/genegalleon:iqtree-speedup-dev \
  --output /tmp/iqtree-speedup-validation/ablation.json
```

各要求の尤度と全枝勾配を保存して照合する。許容差は尤度絶対値 `2e-6`、勾配 `rtol=2e-7, atol=2e-5`。対角の有無以外は同じ要求を使用する。

## 評価カーネルの結果

各セルは10評価の中央値［最小〜最大］、単位は秒。

| 入力 | 変更前 | 固定キャッシュのみ | キャッシュ＋SCORE |
| --- | ---: | ---: | ---: |
| 4配列 × 150コドン | 0.147 [0.143–0.161] | 0.142 [0.137–0.153] | 0.137 [0.135–0.140] |
| 16配列 × 1,500コドン | 9.251 [9.099–10.298] | 9.230 [9.160–9.277] | 9.267 [9.204–9.298] |
| 64配列 × 1,500コドン | 41.409 [40.957–42.846] | 42.357 [41.843–43.571] | 42.647 [42.470–42.996] |

4配列ではSCOREまで含めて中央値で6.5%短縮した。16配列はほぼ同じ、64配列は変更後の中央値が3.0%長かった。64配列では範囲が重なっており、3反復から小さな差の原因を特定はできないが、通常の大規模評価を高速化できたとは結論しない。

27試行の全要求で尤度・全枝勾配は完全一致した（最大絶対差0）。変更後の全180要求は固定キャッシュを再利用し、安定計算への切替は0回。64配列の初期化中央値は43.406秒→43.382秒、IQ-TREE子プロセス最大RSS中央値は226.73→226.79 MiBで、こちらも実質的な削減はなかった。

## marginalプロファイル構築の結果

同じ平衡木と合成された正定値の二次尤度を使い、親問題を構築した後、年代を1点固定した問題を40個構築した。シーケンス尤度のfitと目的関数の評価時間は除く。変更前→変更後、変更後→変更前の順で各3反復、計6観測ずつ。親問題の初回構築は短縮対象外である。

| 配列数 | ρ | 変更前：40構築 | 変更後：40構築 | 倍率 |
| --- | ---: | ---: | ---: | ---: |
| 4 | 0 | 0.0144秒 | 0.0095秒 | 1.5× |
| 4 | 0.5 | 0.0214秒 | 0.0166秒 | 1.3× |
| 64 | 0 | 0.0562秒 | 0.0110秒 | 5.1× |
| 64 | 0.5 | 0.0680秒 | 0.0246秒 | 2.8× |
| 256 | 0 | 1.0834秒 | 0.0201秒 | 54.0× |
| 256 | 0.5 | 1.0696秒 | 0.0581秒 | 18.4× |

全ケースで構築後の尤度と勾配が絶対差 `1e-12` 以内で一致した。これは共通行列を含む問題構築部分の改善であり、年代推定全体やexact-onlyの速度倍率ではない。目的関数評価中の分散依存の分解計算は各問題に残る。

再現用スクリプトはNWKITの `tools/benchmark_radte_marginal.py`。`benchmark_radte_iqtree.py` と同じディレクトリに置いて両イメージへマウントし、`--output` を指定する。両スクリプトは実行先イメージにインストールされたNWKITを使用する。

## 残る計算負荷

16配列×1,500コドン、10回のSCORE要求の独立した分解計測では、約9.188秒のうち `computeLikelihood()` が2.160秒、全枝の微分ループが7.025秒、残りは約0.004秒だった。微分ループには方向付き部分尤度の更新も含まれる。

別途 `-pg` を付けたバイナリのgprofでは、`productVecMat<Vec2d, double, false>` と `dotProductDualVec<Vec2d, double, false>` が主要な計算箇所だった。これは計測用ビルドのサンプルであり、通常ビルドの厳密な時間比率には使わない。ARMではVec2dによるベクトル経路を実行していた。

現在も方向付き部分尤度を枝間で共有している。全枝APIを追加するだけで木全体の再計算を除去できる、という構造ではない。次の大きな改善には61状態の行列・ベクトル計算、メモリ配置、走査とスコア計算の統合を対象にした比較が必要になる。今回、全枝APIの大規模変更や浮動小数点演算順序の変更は実施していない。

## ワークフローの一致

GeneGalleonの実際の年代推定ステージを、変更前後それぞれ既定モデルと `GY+F3X4+R4` で実行した。両ケースとも固定モデル、推定方式、診断、区間の状態が一致し、ノードの年代・95%区間・推定速度の最大絶対差は0だった。species表も一致した。両ケースの最終推定方式は `sequence-empirical-bayes-map`、区間は `conditional-profile` である。

この比較は4配列・150コドンの機能検証であり、各1回のワークフロー所要時間から全体の速度倍率は算出しない。marginal区間については、既存の独立再最適化による尤度比endpoint検証と、共通構造を共有／再構築した値・勾配・事後速度の比較を実行する。

最終runtime `local/genegalleon:iqtree-speedup-final` のIQ-TREEバイナリSHA-256は、比較に用いたcandidateと同じ `a0d22a6fe1d93a0c742522e9ec6617150d918b222280d47ef320f183f14ea09d`。最終検証対象にはNWKIT側の型ガードとドキュメントの更新も含めた。

## 検証成果物

- [コンパクトな測定記録](benchmarks/iqtree-speedup-2026-09-09.json)
- 比較要求・生の結果・image ID: `/tmp/iqtree-speedup-validation/ablation.json`
- 関数プロファイル: `/tmp/iqtree-speedup-validation/gprof.txt`
- 実装元: `/Users/kf/repos/iqtree`、`/Users/kf/repos/nwkit`

## 最終チェック

GeneGalleon Docker runtime内のNWKIT検証用スナップショットで、次を実行した。

- Ruff lint／format、mypy（178ソース）、依存関係整合性、Bandit、pip-audit: 成功。
- 全テスト: **3,269成功、48スキップ、1失敗**（620.48秒）。IQ-TREE／marginalの対象テストは成功。
- 唯一の失敗は `test_interface_conventions.py::test_every_visible_long_option_has_a_canonical_kebab_case_spelling`。別機能の `regress-select --predictor-file` に対してテストがunderscore aliasを要求する。変更前の保存済みソースと変更前イメージでも同じ失敗を再現した。今回そのCLIや規約テストは変更していない。
- 全テストのカバレッジ: **85%**。保守性の上限チェック: 成功。
- wheel／sdistの内容、メタデータ、再現性: 成功。テスト失敗でreleaseチェックが中断したため、カバレッジ以降は個別に実行した。
- IQ-TREE producerのSCORE／EVAL／STATSプロトコルテスト: 成功。同じバイナリを最終runtimeでも使用する。
- GeneGalleonの既定／FreeRateステージを変更前後で実行し、PDF、樹形、manifest、provenanceと固定されたspecies root age=10を確認した。

今回変更したNWKITの実装・テスト・ツール・RADTE文書11ファイルが検証用スナップショットと一致することを確認した。全体チェックを全件成功とは扱わない。SIFおよびx86 runtimeは実行しておらず、確認範囲はDocker arm64である。

開発用 `local/genegalleon:iqtree-session-dev` は最終イメージ
`sha256:4c9840f86183d73e939132ef0e0b3d5cd3ec4761122c6d5bb76f1d32cedd9e08`
へ更新した。比較元の専用タグは保持した。

詳細ログは `/tmp/iqtree-speedup-validation/` の `full-check.log`、
`preexisting-cli-failure.log`、`coverage-report.log`、`maintainability.log`、
`dist-check.log`、`stage-comparison.json` に保存している。
