# 項目11：GBIF観測指標の実装記録

対象は[採用計画](item-11-gbif-plan.md)の出力契約、取得状態、再現可能な取得と記述的感度分析。
利用法・列移行・設定は [GBIF observation traits](../../gbif-observation-traits.md) を参照。

## 実装した変更

- 形質表は10種類の `gbif_observed_*` 観測指標とし、取得件数・照合品質・取得状態を別TSV/JSONに分離した。
  旧来の平面近似凸包は廃止し、球面セル面積を真の分布面積・IUCN AOOとは表示しない。
- table hashとroleを持つmetadataを生成し、取得query、taxon照合、UTC、正規化レコードsnapshot/hash、
  採用・除外理由、未知の品質項目、dataset/地域/観測種別/由来/測地系の内訳を保存する。
  出力一式は一時ファイルから公開し、公開途中の失敗時には以前のファイルへ戻す。
- searchの上限、ページ位置・件数変動・短いページ・終了フラグ・ID重複・種key不整合を検査する。
  部分取得/未解決照合は記述用sidecarに残し、本表は欠測にする。完全取得必須モードは旧出力を置換せず失敗する。
- 保存済みGBIF SIMPLE_CSV/TSV、gzip、単一テーブルZIPを読める。既存download metadataと全行数を照合し、
  未検証ファイルを完全なdownloadとは扱わない。取込はネットワークを使わず、downloadの新規申請はしない。
- 日付・地域・basis・自然由来・不確実性を明示的に選べる。
  旧centroid距離オプションを、既知の地理参照中心点への**最小**距離に改めた。未知を在来や高精度へ変換しない。
- 日付変更線、極、端数格子、定義不能な円平均を処理する。
  同一gbifIDは重複除去し、別イベントの同一座標はrecord平均に残す。
- コピー数PGLS・RSC/species PGLSは `all` で観測/品質列を選ばない。
  明示指定した観測列は別名でも意味と取得品質を保持し、不適格な種を欠測にする。
  コピー数結果には選択表/metadata/選択監査JSONと、観測値の関連であることを示す図captionを含める。
  codeml/HyPhy/CSUBSTのforegroundには観測指標を自動変換しない。
- オフライン感度分析は格子、期間、地域、日付不明、自然由来、dataset除外、record/位置/セル単位の抽出、
  固定標本数を比較する。snapshot、seed、条件、採用ID集合hashを記録し、標本不足は欠測にする。

## 初回実装の検証範囲

取得境界/不正応答、filter/unknown、全球セル面積の独立解析解、日付変更線・極・円平均、
cache再利用/更新/破損、圧縮download、出力rollback、別名role・部分取得の除外、
固定seed再解析を自動試験に含めた。runtime manifestへGBIF生成→コピー数PGLS→RSC入力の統合試験を追加した。

コンテナ検証結果（2026-09-10）：

| 検証 | 結果 |
|---|---|
| `dev check full` の全Python試験 | 2,348件中2,344件成功、初回4件失敗。20分11秒 |
| GBIF・生成器・入力生成coreの最終コード | 95件成功。新規統合試験の期待列名を修正し、同試験の再実行も成功 |
| CSUBSTの形質/色付け/候補出力 | 69件成功 |
| 実CSUBSTコマンド/結果取込 | 4件成功 |
| 修正後のGBIF統合試験＋shell static safety | 237件成功 |
| 修正後のASR core cache＋上流source契約 | 3件成功 |
| `dev check r` | manifestの必須15コマンドすべて成功 |
| hostでの変更ファイルRuff・Bash構文・ShellCheck・設定schema・差分空白検査 | 成功 |

初回4件は、新規統合試験の `analysis_method` / `species_nwkit` 参照、変更したentrypoint説明文、
古い次ステージ見出し `l1ou`、上流ソース集合に欠けていた `iqtree` という試験期待値の問題だった。
後者2件は変更前のコードにもある不整合で、現在の見出し・必須ソース集合へ追従させた。
検査範囲は削減せず、4件すべて修正後の対象再試験で成功を確認した。全Python試験の一括再実行はしていない。
コピー数のSVGは描画して、観測値に関するcaptionが収まることも確認した。

検証コンテナは `local/genegalleon:gbif-observations-dev`（linux/arm64）。
イメージIDは `sha256:226b336b916ad91a271f51d458c24cd0b1d57ba2d6198592eef99ebc13a4c9b3`、
runtime-inputは `b4e8ea38741783f045c4f0b8e1b6ee001b37f963f669fc63b6e388c9daae6ab4`。
検証wrapperで現行入力とdaily owned upstream snapshotの一致を確認し、鮮度検査の無効化は行っていない。

2026-09-10 01:33 UTCの実APIプローブでは、Arabidopsis thaliana（speciesKey 3052436）を
取得上限1件で照合・検索した。当該queryの報告件数303,775、取得1件、`capped_partial`、
`analysis_eligible=false`、本表の観測指標すべて欠測を確認した。
これはAPI接続・上限状態の確認であり、分布推定や代表標本の検証ではない。

全体のhost lintには、今回変更していない `repair_coge_transcripts.py`、`test_real_aha_dataset.py`、
`test_repair_coge_transcripts.py`、`test_synteny_cutoff_metadata.py` の既存19件の指摘がある。
今回変更したPythonファイルのRuff、設定schema検査、差分空白検査は合格した。

## 科学的な限界と別段階の作業

`complete` は保存したquery/downloadの取得完了を指す。GBIF全体や自然分布の完全性、ライブsearchの
厳密な同時点snapshot、不偏性を意味しない。保存したdownload metadataは利用者提供の証跡であり、
同一行数のファイルが公式配布物と同一であることをリモート照合したとは主張しない。

観測努力だけで平均・境界・面積が変わる負の対照を含むが、これは推論法の校正ではない。
計画の項目10/12と共通化する2,000反復以上の系統的交絡シミュレーション、Type I error/被覆率/FDRの校正、
ゲノム品質共変量・測定誤差モデルは本変更で実装・検証済みとはしない。
コピー数同時選択への `species_trait_contract` 接続は、下記の本体統合レビューで追加した。
観測努力・ゲノム品質の共有交絡を含む推論校正は、引き続き別段階の作業である。
自然分布推定には研究対象・地域/時期・独立調査/努力量データと妥当性検証が別途必要であり、新しい推定器は追加していない。

実GBIFの大規模取得、download申請、Linux/HPCのSIF実行は未実施である。
この作業環境にはApptainer/Singularity実行系と `genegalleon.sif` がなく、Docker結果をSIF互換性とは呼ばない。


## 本体統合前の再レビュー（2026-09-10）

本体 `main` の `db7e6ba` を基に隔離チェックアウトで統合し、並行作業完了後の `18c134f`
（MCMCtree並列chain・形質型schema追加を含む）へ再統合した。
GIFT取得改善、形質欠測保持、コピー数応答分布/同時選択、表現量PGLSの多重補正を保持した。
本体の作業中の `orthogroup_statistics.py` と対応試験、未追跡レビュー資料は変更対象に含めない。
同じcore/entrypointで進行中のnative OU追加も、行単位で分離して未ステージのまま保持する。

追加で修正した問題：

- 形質表・品質/記述sidecar・metadata・統計JSONの公開を一つにまとめ、途中失敗で旧出力一式へ戻す。
  GBIF原データ、入力設定、cache内、別出力への上書きも事前に拒否する。
- 列数不足/超過、正規化後の重複ヘッダー、hard linkで同一実体を指す出力を拒否する。
  cache manifest自体と集計値の破損も検知して再取得する。
- 未使用GBIF設定の事前取得を避け、array準備時には間接指定されたGBIF入力ファイルもhash化する。
- コピー数同時選択でも `all` 除外・明示指定・部分取得種の欠測化・監査付き選択表を接続する。
  各形質の選択結果metadataにも観測値の意味と除外理由を保存する。
- RSC準備の全出力をrollback可能にし、既定のRSC-only経路でも最終method-auditへ入力監査を残す。
  一時準備ディレクトリを削除しても科学的な意味・除外理由が失われない。
- 前景用二値化が普通の形質の欠測/未知値を背景0へ置換していた処理を修正する。
- 前回記録した既存lint 19件も、import整理・文分割・長さ確認済みzipのstrict指定で解消する。

GBIF取得/生成器、GIFT、破損とrollback、RSC、コピー数同時選択、core二値化の回帰試験を追加・再実行した。
形質型schemaとの統合では、コピー数PGLS・同時選択の観測制約と数値型検証を順に適用する。
型による除外と観測制約による除外を `trait_selection.tsv` に保存し、派生形質表には
metadataとschemaの両方を付ける。既存のtext/categorical除外・不正コピー数拒否を保持した。

再レビュー中の全体実行は2,455 passed / 4 failed（662.76秒）。失敗4件は次の理由を確認して修正した。

- 現行NWKITの推定分散ゼロ境界ではprofile区間が定義できないため、年代推定の旧期待値3件を更新した。
  `unavailable-strict-clock-limit`、境界診断、全区間の欠測を検査し、正分散では従来のprofile検査を残す。
- 既存GIFT改善コミット `2024612` が更新した入力設定に対し、実データmanifestの保存hashが古かった。
  対応する設定のhashのみ更新した。

両試験ファイルの再実行は **12 passed**（223.71秒）。全Python suiteを修正後に一括再実行したとは主張しない。
Rパッケージ検査で新規gene-structure処理の裸の `head()` がNOTEとなったため、`utils::head()` に修正した。
R必須15コマンドはすべて合格。型schema統合後のコピー数R試験も別途合格した。

統合後のGBIF・生成器・GIFT・型schema・コピー数同時選択重点試験は **178 passed**（58.03秒）。
最後のmetadata公開順・派生schema確認を含む再試験は **2 passed**（9.54秒）。
追加されたMCMCtree・genome core・shell静的回帰試験は **295 passed**（471.89秒）。
変更対象の全シェルもリポジトリ指定除外付きShellCheckに合格した。

host `dev lint`（全Python Ruff、8 entrypoints / 23 common parameters、Bash構文）は合格。
本体に追加されたMCMCtree関連のimport順と例外連鎖指定のlint指摘も最小変更で解消した。

最終Docker image: `local/genegalleon:gbif-main-review-dev`

- image: `sha256:6aac044b761c9f888f7a0996055ee0383602ea324ec4717f7f0e80867011a30d`
- runtime-input: `0ffbdedd69bbf00c99826efaa38d3d065e827ce10d8f78fde275351ae6ebcd3b`

各wrapperは現行container入力とdaily snapshotの一致を確認した。NWKITの同時選択を含む検証では、
所有repositoryのコミット済み `e2465b048f7529b7dad84cb6722278406bb49bbe` を `git archive` で取り出し、
read-only bindと `PYTHONPATH` で使用した。未コミット変更は含めていない。
これはローカル依存統合検証であり、公開NWKIT・公開コンテナだけで同時選択が動作するとの主張ではない。
公開側へ同機能を配布する際は所有NWKIT側の公開更新が必要であり、本変更に旧版fallbackや固定SHAの既定値は追加していない。
