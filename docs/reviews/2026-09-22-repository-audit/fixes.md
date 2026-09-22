# 調査指摘への対応 — 0.8.0

[修正前の調査](review.md)で挙げた6項目に対応した。既存coreの構成・実行順序と、既存SQLiteの読み取り互換性を維持する。

| 項目 | 対応 |
|---|---|
| 統合fixtureとscratch契約 | core単独のfixtureにworkspaceモードを明示。supervisorを通る外部scratchテストは独立して保持 |
| fullの3Di skip矛盾 | strict manifestで実3Diを有効化。Dockerの実行UIDから書き込めるモデルキャッシュを設定し、未許可skipの失敗条件を維持 |
| GOの観測後選別 | `cafe_branch_flags_all_go` を追加。対象イベントと独立に背景GO集合を確定し、1ファミリー1観測、ゼロ該当GOを含むBH補正、bothでは増減をまとめた補正を実装 |
| 誤ったDBフィルタ | 不正書式、未知列、NaN/Infinityの閾値を拒否。空文字でのみ明示的にフィルタ解除。欠損値は数値閾値を通過させない |
| DB再実行の重複 | 既定のcreateは既存DBを拒否。replaceは既存の `--overwrite 1` を保持。明示的appendは重複ファミリーを拒否し、PGLS既存行を保持して補正を再計算。すべて一時DBで処理し、公開更新を直列化 |
| macOS lint | Bash 4+の選択・要件診断を追加。Git列挙失敗を検出し、Bashラッパーの標準入力読み取りで検査対象が欠落しないよう修正 |

GOの従来 `event` / `cafe_branch_flags` の数値・出力パスは維持する。event利用時には探索的解釈を警告する。新方式の最小例はGO10項目を保持し、調整P値0.0822487644を再現する。これは対象観測後のGO除外を解消する修正であり、CAFEによる選別・注釈依存性を含む解析全体の校正を完了したという主張ではない。結果metadataにも探索的解釈を記録する。

CLI移行上、従来の暗黙のDB追加を利用していた場合は `--mode append` を明示し、追加する新規ファミリーだけを入力する。既存の標準workflowは `--overwrite 1` を使うため呼び出し変更は不要。

## 検証環境

他の作業によるテスト整理が同じ作業ツリーで進行していたため、HEAD `716261c` から隔離したチェックアウトを作り、今回の変更のみを適用した。別作業のテスト削除・変更は今回のcommitに含めない。検証中にその整理が `ec539e5` としてmainへ公開されたため、今回の修正はその上へ適用する。隔離検証では整理前のより広いテスト群を維持している。

2026-09-22にビルドした `local/genegalleon:dev` を使用：

- image ID: `sha256:509ccb36468841a78f963568534437922c67b75be058eab8c20024907414ec1b`
- runtime-input: `443934b758122d5aa44a4177d02d11d79d40c339773d1296ed1f87d901231337`
- Docker Linux arm64。SIF/HPCの実行結果ではない。

隔離チェックアウトでも同じ日次上流スナップショットを比較するため、`GG_RUNTIME_FRESHNESS_CACHE_DIR` に元チェックアウトの当日キャッシュを指定した。freshness検査は有効。上流が調査中に進んだため別の再ビルドも開始したが中止し、最初に正常ビルドした上記スナップショットを検証対象として固定した。SHAは診断記録であり、上流の既定ブランチ設定を変更していない。

## 最終検証結果

| 検証 | 結果 |
|---|---|
| `dev check full --tb=short` | **2,823 passed、skipなし**、11件の依存側/既存fork使用に関するDeprecationWarning |
| fullに含まれるRチェック | **16コマンドすべて成功**、treevisの `R CMD check` はStatus: OK |
| 隔離環境のDB・PGLS・開発ツール追加確認 | 98 passed |
| 新GO方式の実CAFE Base/Gamma連携 | 2 passed |
| 最新mainとの重複箇所（開発ツール） | 33 passed |
| host lint（コンテナBashパーサ使用） | Bash構文・Ruff・8 entrypoints/23 common parametersの設定検査が成功 |
| 変更shellのShellCheck | warning severityで成功 |
| 差分・報告書リンク | 成功 |

[fullログ](evidence/fixed-full.txt)、[追加確認](evidence/fixed-focused.txt)、
[実CAFE](evidence/fixed-native-go.txt)、[main重複箇所](evidence/fixed-main-overlap.txt)、
[lint](evidence/fixed-lint.txt)、[検証ランタイムの上流revision記録](evidence/runtime-sources.tsv)。

新しいコンテナを使いfreshness検査を有効にして検証した。SIF/HPCおよび大規模生物データによる統計的校正は実施していない。

保存ログは行末の空白と末尾の空行のみ正規化した。
