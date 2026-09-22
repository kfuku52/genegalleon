# リポジトリ調査・改善提案 — 2026-09-22

> この文書は修正前の調査記録です。修正内容と最終検証は [対応記録](fixes.md) を参照してください。

対象：GeneGalleon 0.7.133、HEAD `716261c07f467ff4a4eea91d635cf019cf84cb2e`。調査開始時の追跡ファイルに未コミット変更なし。実装は変更せず、本報告と証拠だけを追加した。

優先すべきなのは、検証基盤の不整合の修正と、解析条件を黙って無視する挙動の解消である。単体テストや保存形式の防御は充実しているが、単体テストの成功だけでは統合動作や統計的解釈を保証できない。以下の6項目は、今回の再現結果に基づく。GOの問題は新しい回帰ではなく、互換性のため残されている既知の制約である。

## 調査範囲と限界

- README、開発・実行・ランタイム文書、10本のcoreスクリプトの構成を確認し、起動、一時領域、provenance、ロック、ZIP公開、DB生成、検定補正、テスト選択を重点的に追跡した。全行を精査したという意味ではない。
- 標準検証入口、静的チェック、fast・static・download・workflowレーン、宣言されたRチェック16コマンド、問題別の最小例を実行した。
- macOSの標準Bashは3.2.57。Docker `local/genegalleon:dev` は2026-09-20作成、image IDは `sha256:794ec414717b243f72a39e88d4ae8698f568cd06e5327e8bbef7a484bad8c251`。
- 標準 `bash ./dev check fast` はfreshnessチェックで停止した。観測runtime-inputは `4d0ce409a58ef0936880c78b84bbb448f6d1d5f3f2ffe1b4935cb647611c3861`、期待値は `554748f7a80a7285fe93ed75866c48c8b9621562d23a00b13f6ca4fa5afe333c`。これは古い環境の利用を検出した正常な防御動作であり、それ自体をコード不具合とは判定しない。
- 以後のDocker検証は、原因切り分けのため明示的に `GG_RUNTIME_FRESHNESS=off` とした。チェックアウトをマウントして現在のコードを実行したが、現在の依存ブランチ先端との互換性を保証する結果ではない。
- 新規コンテナビルド、全runtime Pythonレーン、実3Di予測器、SIF、HPC実ジョブ、大規模生物データでの一貫実行・統計的校正は未実施。最新GitHub Actionsの成否も照会していない。

## 実行結果

| 検証 | 結果 | 解釈 |
|---|---|---|
| `bash ./dev lint` | 失敗 | macOSのBash 3.2によるcore構文検査で停止 |
| Ruff | 全件成功 | `ruff check workflow container/scripts --output-format concise` |
| config schema | 成功 | 8 entrypoints、23 common parameters |
| cache guard audit | 検出0 | 自動検出パターン内の結果。全キャッシュの正しさの証明ではない |
| fast | 1,934 passed、8 warnings | [ログ](evidence/fast.txt) |
| static | 385 passed | [ログ](evidence/static.txt) |
| integration-download | 87 passed、2 skipped | 実3Diテスト2件がopt-in未指定でskip。[ログ](evidence/download.txt) |
| integration-workflow | 7 passed、47 failed | scratch契約に追随していないfixture。[ログ](evidence/workflow.txt) |
| 同上、診断用workspace指定 | 43 passed、11 failed | genome系の36件が改善。入力生成側は環境変数を独自生成するため指定を継承しない。[ログ](evidence/workflow-workspace.txt) |
| R manifest | 16コマンドすべて終了コード0 | treevisパッケージチェックを含む。[ログ](evidence/r.txt) |
| strict + 実3Diファイルのみ | 2 skipped、終了コード1 | skip禁止とopt-in既定値の矛盾を再現。[ログ](evidence/strict-3di.txt) |

Rログには異常入力を拒否するテストの意図したエラー出力も含まれる。コマンド終了コードとテストの成功メッセージで判断した。Rテストがrepo-rootに生成した `Rplots.pdf` は今回の一時成果物として退避済み。

## 1. P1：scratch既定値変更で統合テストが起動直後に止まる

根拠：`workflow/gg_common_params.sh:8`、`workflow/support/gg_util/05_workspace_io.sh:808`、`workflow/tests/test_genome_evolution_protein_mode.py:807`、`workflow/tests/test_gg_input_generation_end_to_end.py:396`。

既定の `GG_COMMON_TMP_ROOT=/tmp` はentrypoint supervisorが設定する `GG_TMP_TASK_ROOT` を必要とする。しかし既存の統合fixtureはcoreを直接起動し、その契約を用意していない。結果は `External scratch requires the GeneGalleon entrypoint supervisor.`。解析本体に到達する前にテストが終了する。

診断としてDocker内に `GG_COMMON_TMP_ROOT=workspace` を与えると、47失敗から11失敗に減った。残る入力生成fixtureは許可したキーだけで環境を組み直すため、この指定を継承しない。一方、supervisorを実際に通す `test_task_tmp_integration.py` の2件は成功した。したがって、通常のentrypointもすべて起動不能だという指摘ではない。

改善案：core単独のテストにはworkspaceモードを明示し、supervisorとの結合を検証するテストは本来の起動経路を通す。fixtureの環境構築を揃え、通常モードと外部scratchモードの両方を検証する。安全ガードを削除して通すべきではない。

完了条件：標準環境のintegration-workflow 54件が成功し、外部scratchの独立2件も成功すること。新しいコンテナでsmoke、関連runtime、fullを再確認する。

## 2. P1：`dev check full` の必須検証と3Diのopt-in設定が矛盾する

根拠：`workflow/tests/run_checks.py:36`、`workflow/tests/conftest.py:69`、`workflow/tests/test_csubst_3di_runtime.py:22`、`workflow/tests/validation_manifest.json:3`。

fullは全Pythonテストを `--gg-strict-runtime` 付きで実行する。strictでは未許可のskipが失敗になる。一方、実3Diテストは `GG_TEST_CSUBST_3DI=1` がなければ必ずskipし、共通manifestはその変数を設定せず、allowed_skipsも空である。該当ファイルだけのstrict実行でも、2 skipと終了コード1を再現した。full全体は未実行だが、この組合せはそのPythonフェーズに含まれる。

CIのSIF actionには3Di専用の別実行があり、そこでのopt-inは明示されている。したがって「CI全体で3Diが未検証」とは言えない。問題はローカルfullと専用実行の契約が一致しない点である。

改善案：fullが実3Diを含むなら専用設定・キャッシュを宣言して確実に実行する。別レーンとするならfullの仕様・選択規則を明文化し、集約した必須検証で実3Diの成功を要求する。skipを無差別に許容する修正は避ける。

完了条件：未宣言のskipで失敗する防御を保ちつつ、文書どおりのfullコマンドと実3Diチェックが両立すること。

## 3. P1・既知の制約：GOの観測後選別で補正対象が縮み、有意性が変わる

根拠：`workflow/support/cafe_go_enrichment.r:165` と `:209`。GO解析は `run_go_enrichment=0` が既定で、ユーザーが有効化した場合の既定方式が `event`。影響を全ワークフローに一般化しない。

`event` は対象群に出現しないGOを補正前に落とす。以前のレビューに付属する再現コードを現在のソースに対して再実行し、次を確認した。

| 指標 | 現在の結果 |
|---|---:|
| 入力中のGO項目 | 10 |
| 補正に残るGO項目 | 2 |
| 対象GOの未補正P値 | 0.008224876 |
| 現行のBH調整P値 | 0.01644975 |
| 他の8項目のP=1を含めたBH調整P値 | 0.08224876 |

0.05での結論が逆転する。[再現ログ](evidence/go.txt)。RのBH実装は既定で渡されたP値の個数を比較数として使う。[R公式文書](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)。この最小例は補正母集団による違いを示すもので、実データの偽陽性率を測定したものではない。

`docs/go-enrichment.md` は互換性維持と探索的解釈を明記しており、追加の `cafe_branch_flags` も既存の対象GO集合を保持する。この既知の設計を未認識の新規バグとして扱うべきではないが、確証的な解析に使用する際には優先して解決する必要がある。同一ファミリーの複数イベントを数える方式と、1ファミリー1観測の方式も区別が必要。

改善案：既存出力の互換性を保持したまま、事前に検定集合を定義する別の解析モードを設計する。対象群ゼロのGOも含め、観測単位・背景・補正範囲・探索的解釈を結果自身に記録する。BH計算の変更だけで解析全体のFDRが保証されるとはしない。

完了条件：固定した検定集合の回帰例と、帰無データでの選択・依存構造を含む校正結果を整備する。既存方式の数値を無断で変更しない。

## 4. P2：DBの絞り込み条件の誤記が正常終了のまま結果を変える

根拠：`workflow/support/generate_orthogroup_database.py:591` と `:616`。

数値として解釈できない閾値は `continue` で無視され、存在しない列名も無視される。NaNは数値として受理される。OCN値が0.1と0.9の2行で再現した。

| 条件 | 残る値 |
|---|---|
| `OCNany2spe,0.8` | 0.9 |
| `OCNany2spe,typo` | 0.1、0.9 |
| `OCNany2sp,0.8`（列名誤記） | 0.1、0.9 |
| `OCNany2spe,nan` | なし |

`test_parse_cutoff_stat_parses_valid_tokens_and_ignores_invalid` は無視する挙動を現行仕様として固定している。テスト成功が望ましい利用者体験を保証しない例である。

改善案：CLI境界で書式と有限値を検査し、各対象テーブルで指定列の存在を確認する。明示的な無フィルタ指定と不正入力を区別し、診断付きで停止する。既存の寛容な公開挙動を変更するなら、厳格モードの追加と移行方針を用意する。

影響範囲：直接CLI・独自設定で特に問題になる。標準gene-summary呼び出しは固定の正しい閾値を渡しており、通常の全実行で発生するという意味ではない。

完了条件：誤記、未知列、NaN、Infinityを検出し、エラー時に既存DBを保持する。[最小例](evidence/reproduce_database.py)、[結果](evidence/database.txt)。

## 5. P2：DB生成CLIの既定の再実行で同じ行が重複登録される

根拠：`workflow/support/generate_orthogroup_database.py:160`、`:632`、`:976`。

`--overwrite` の既定値0では既存DBをそのまま開き、入力行をappendする。同じ最小入力で同じコマンドを2回実行すると、両方が終了コード0で、tree・branchテーブルが各1行から各2行になった。既存DBを保護する「上書きしない」と、重複しない「再開」は別の意味だが、CLIヘルプでは明確に区別されていない。

標準coreは `--overwrite 1` を指定しており、この重複の直接対象は主に補助CLIの利用者である。appendが意図された利用法を否定するものではないが、同一入力の再投入を検出しないことは集計の信頼性上のリスクになる。

改善案：create・replace・appendの意味を明確にし、appendには入力識別と重複検出を持たせる。既存引数を互換レイヤーとして残し、既存DBへの意図しない追加に明瞭な診断を出す。新規作成にも一時DBからの公開を適用し、失敗した途中DBを完成品と混同させない。

完了条件：同一入力再投入で行数が増えないか、明示的な重複エラーになること。追加・置換の失敗時に公開DBを保持すること。再現は項目4と同じ証拠ファイルに含む。

## 6. P2：macOS向けhost lintがコンテナ用BashスクリプトをBash 3.2で検査する

根拠：`dev:53`。73本の追跡shellファイルをmacOS標準Bashで調べると、`gg_transcriptome_generation_core.sh` のみ構文検査が失敗する。GeneGalleon Docker内のBashでは同じファイルが成功する。

`dev lint` は最初の失敗で終了するため、後続Ruffとconfig schemaにも到達しない。それらを独立実行すると成功した。したがって、表示されたEOFエラーをそのままLinuxの本体スクリプト破損と解釈してはいけない。

改善案：hostで実行するラッパーと、コンテナ内で実行するcoreの構文検査対象・Bash要件を分ける。host lintで新しいBashを必要とするなら実行前に明示し、検証用Bashを選択できるようにする。coreのパーサ互換性をmacOSの旧Bashに合わせるためだけに変更しない。

完了条件：macOSとLinuxで、同じ実行環境を想定したlint結果が得られること。Bash不足時はファイル末尾の構文エラーではなく、必要なバージョンと対処を表示すること。

## 実施順序と継続的改善

1. **検証の回復**：項目1・2・6。正しいscratch契約、レーン選択、検証用shellを整える。最新コンテナを構築し、freshnessチェックを有効にした状態で関連レーンとfullを通す。
2. **入力ミス・再投入への対処**：項目4・5。明示的な失敗と公開DB保護を先に実装し、既存CLIの移行を文書化する。
3. **科学的解釈の明確化**：項目3。探索的方式を維持しつつ、確証的な追加方式の仕様と校正を別途設計する。
4. **変更影響の追跡**：provenanceは入力とパラメータを中心に比較し、ツール・producerバージョンだけの変更では失効させない設計である。これを一律バージョン失効に変えるより、結果の意味を変える修正ごとに対象成果物の再生成条件・移行手順を記載する運用を整える。
5. **保守性**：coreは合計24,170行。承認なしにステージ別ファイルへ分割せず、既存core内の責務・入力・出力・更新条件とテストの対応表を整備する。文字列検査に加え、途中失敗→再開、入力変更→失効、ZIP/liveの混在、並列実行などの観測可能な挙動を増やす。

既存の強みとして、provenance、ZIP検査と原子的公開、共有ロック、入力スキーマ検査、Rのパッケージ検査、上流のmoving branchを検査する仕組みがある。これらを弱める回避策ではなく、今回見つかった検証契約・CLI境界の穴を埋める方針を推奨する。性能の優劣は測定しておらず、速度向上やメモリ削減の数値的主張は行わない。

## 再現コマンド

以下のfreshness無効化は今回の診断条件を再現するためであり、リリース検証にそのまま使わない。

```bash
bash ./dev lint
GG_RUNTIME_FRESHNESS=off bash ./dev check fast
GG_RUNTIME_FRESHNESS=off bash ./dev check static
GG_RUNTIME_FRESHNESS=off bash ./dev check integration-download
GG_RUNTIME_FRESHNESS=off bash ./dev check integration-workflow
GG_RUNTIME_FRESHNESS=off bash ./dev check r
GG_RUNTIME_FRESHNESS=off bash workflow/tests/run_in_runtime.sh \
  env GG_COMMON_TMP_ROOT=workspace python workflow/tests/run_checks.py integration-workflow --workers 2
GG_RUNTIME_FRESHNESS=off bash workflow/tests/run_in_runtime.sh \
  python -m pytest -q --gg-strict-runtime workflow/tests/test_csubst_3di_runtime.py
GG_RUNTIME_FRESHNESS=off bash workflow/tests/run_in_runtime.sh \
  Rscript docs/reviews/2026-09-10-scientific-validity/reproduce_go_filter.R "$PWD"
docker run --rm -v "$PWD:/repo:ro" -w /tmp local/genegalleon:dev \
  python /repo/docs/reviews/2026-09-22-repository-audit/evidence/reproduce_database.py
```
