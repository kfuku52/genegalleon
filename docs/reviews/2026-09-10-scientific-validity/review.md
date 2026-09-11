# GeneGalleon 科学的妥当性レビュー — 2026-09-10

> 履歴資料：以下は2026-09-10のレビュー時点の記録です。その後の修正を反映した現在の不具合一覧ではありません。各項目の対応記録と現在の実装を参照してください。`validation.txt` は当時の実行結果で、`reproduce_go_filter.R` の同じ数値を再現するにはレビュー対象時点のソースが必要です。

科学的な結論に直接使うには修正または追加検証が必要な機能がある。特に、GO富化の多重検定補正、MCMCの収束確認と証跡保存、実験的な遺伝子年代推定、イントロン消失からの機構判定、探索後の有意性解釈が重要である。一方、ワークフロー全体が科学的に不適切という結論ではない。

対象は HEAD `15481838663118b1bb0ff4e456a0be46ab3c34eb` に未コミット変更を含む現在の作業ツリー。実装・設定・関連テスト・文書を静的に横断調査し、重点箇所を原論文・公式資料と照合した。NWKITの説明はローカル `/Users/kf/repos/nwkit` の文書も参照した。このチェックアウトと実行コンテナの実装がすべて同一とは仮定しない。解析コードや既存のユーザー変更は変更していない。

「確定」はコード上の挙動を確認したという意味で、実データの偽陽性率まで測定したという意味ではない。「要検証」はモデルの校正や入力条件への頑健性が未確認という意味。P1は論文の主要結論に用いる前に対応すべき問題、P2は使用条件・解釈・検証の整備が必要な問題とする。

## 1. P1・再現済み：GO項目を観測結果で選別してからBH補正する

場所：[cafe_go_enrichment.r](/Users/kf/repos/genegalleon/workflow/support/cafe_go_enrichment.r:165)。対象枝に該当イベントがないGOを `next` で除外し、残ったP値だけにBH補正を行う。同じ片側検定で対象群の出現数がゼロならP値は1だが、それを検定数から除くと補正が軽くなる。フィルタが対象群の結果そのものに依存している。

実際の関数をGeneGalleon Docker内で呼び出した最小例：対象10イベント、その他90イベント。あるGOは対象4、その他6に存在する。全イベントに共通のGOを1項目、対象0のGOを8項目用意した。全イベントにGO注釈があるため、上流の注釈済みファミリー限定にも整合する。

| 指標 | 結果 |
|---|---:|
| 片側Fisher P値 | 0.008224876 |
| 現実装の補正対象数 | 2 |
| 現実装のBH調整P値 | 0.01644975 |
| 全10項目を含むBH調整P値 | 0.08224876 |
| 0.05での判定 | 現実装は有意、全項目補正では非有意 |

修正案：対象群の観測結果を見る前にGO検定集合を確定し、ゼロ出現の項目も補正母集団に含める。背景での注釈数などによる事前フィルタを使う場合も、帰無仮説下での妥当性を検証する。[独立フィルタの原論文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2906865/)。再現コードは [reproduce_go_filter.R](/Users/kf/repos/genegalleon/docs/reviews/2026-09-10-scientific-validity/reproduce_go_filter.R)。

## 2. P1・構造を確認：GO検定で同一ファミリーの複数枝イベントを独立に数える

場所：[cafe_go_enrichment.r](/Users/kf/repos/genegalleon/workflow/support/cafe_go_enrichment.r:361)。FamilyID×枝の表を展開し、変化した各イベントをFisher検定の1観測にする。同じファミリーが複数枝で変化すれば、そのGO注釈を繰り返し数える。対象枝と背景枝の双方に同じファミリーが入る場合もある。

これは「対象枝の変化ファミリー対その他ファミリー」という通常の富化検定とは異なる。ファミリー固有の変化しやすさや、同じ木から復元した祖先数の依存性が検定に入らず、反復イベントの多いファミリーが背景に強く寄与する。方向・大きさはデータ依存で、今回偽陽性率の実測はしていない。

修正案：科学的な帰無仮説をまず明示する。ファミリー富化なら1ファミリー1観測の重複しない比較集合にする。枝間のイベント富化なら、ファミリーごとのイベント構造を保存する置換や階層モデルで校正する。項目1の修正だけではこの問題は解消しない。

## 3. P1・確定：種年代MCMCが収束診断なしで採用され、生サンプルも通常削除される

場所：[gg_genome_evolution_core.sh](/Users/kf/repos/genegalleon/workflow/core/gg_genome_evolution_core.sh:3998)。MCMCtreeを1回呼び、FigTreeブロックを抽出できれば公開結果に進む。通常は作業ディレクトリを削除し、公開する「raw」出力には換算後のFigTreeブロックだけを残す。独立チェーン間比較、ESS、R-hatを採否条件にする処理は該当経路にない。

有限の年代と95%区間が出ることは収束の証拠ではない。後段のRADTE、CAFE、PGLSにもこの年代が伝わる。公式PAMLは少なくとも2回の独立実行による確認を求めている。[MCMCtree公式説明](https://github.com/abacus-gene/paml/wiki/MCMCtree)。

修正案：独立seedの複数チェーン、生サンプル・制御ファイル・ログの永続保存、収束診断と明示的な未収束状態を追加する。暫定的に `GG_KEEP_MCMCTREE_RAW_DEBUG=1` と `delete_tmp_dir=0` で証跡を保全できるが、それだけで収束は確認できない。今回実データのMCMC再実行はしていない。

## 4. P1・依存側も実験的と明記：遺伝子年代推定の精度・95%区間を一般化できない

場所：[gene-tree-dating.md](/Users/kf/repos/genegalleon/docs/gene-tree-dating.md:1)、[NWKIT RADTE.md](/Users/kf/repos/nwkit/RADTE.md:8)、[検証記録](/Users/kf/repos/nwkit/RADTE_VALIDATION.md:1)。`run_tree_dating=1` の既定エンジンは独自のnative推定器。NWKIT自身が実験的で、広範な精度・区間被覆率検証が必要と明記している。既存の数値検証や小規模シミュレーションは価値があるが、幅広い遺伝子ファミリーに対する95%被覆率の証明にはならない。

種年代・reconciliation・一度推定した置換モデルを固定するため、通常の出力区間はそれらに条件付けたもの。種年代CIのTSVを追加しても表示用であり、不確実性は伝播されない。遺伝子系統のcoalescence時間と種分化時間が違う状況、誤った重複／種分化の対応、系統間速度変化も検証対象となる。IQ-TREEを選んでも年代制約・時計推定・区間構成はNWKITに残るので、この懸念は消えない。

必要な検証：既知真値の多数の独立ファミリーで、配列長・種数・重複位置・喪失・速度変化・誤reconciliationを振り、バイアス、RMSE、区間の利用可能率、失敗を含めた被覆率を評価する。主張は当面「固定した前提下の探索的年代推定」に限定する。実装修正が必要なら所有元はNWKITである。

## 5. P1・確定：イントロン有無の確率差から「Retrotransposition」と表示する

場所：[ASRアダプター](/Users/kf/repos/genegalleon/workflow/support/asr_intron_evolution.py:65)、[枝分類](/Users/kf/repos/genegalleon/workflow/support/treevis/R/04_branch_synteny.R:18)、[凡例](/Users/kf/repos/genegalleon/workflow/support/treevis/R/00_core.R:386)。イントロン数を有無に縮約し、既定ではgain=0.0001、loss=0.001、根事前分布=(0.5,0.5)の固定CTMCで復元する。重複の片側で `delta_intron_present <= -0.5` ならRとし、凡例は「Retrotransposition」。

確率差は枝上の機構の事後確率ではない。全イントロン消失、DNA重複後の欠失、注釈漏れでも同じ観測を作れる。1イントロンと多数イントロンも同じ状態となる。固定Qの数学が正しくても、そのrateと二状態近似の生物学的妥当性は別問題である。

修正案：表示を「イントロン消失を伴う重複候補」等にし、機構推定と分離する。位置を対応付けたイントロン、親コピー、挿入部位、利用可能なpoly(A)/target-site duplicationの証拠を別々に記録する。rate・根事前分布・年代への感度分析が必要。古いretrocopiesでは挿入痕跡が失われ得るため、痕跡の不在だけで否定もしない。[retrocopyの証拠と限界](https://pmc.ncbi.nlm.nih.gov/articles/PMC5470649/)。

## 6. P1・既定値を確認：CSUBSTの探索補正を計算しても、候補選択は解析的P値由来

場所：[候補選択既定値](/Users/kf/repos/genegalleon/workflow/gg_gene_summary_entrypoint.sh:112)、[支持数別の再補正](/Users/kf/repos/genegalleon/workflow/support/plot_csubst_aa_change_summary.py:466)。scan側は既定で `full_scan` の経験的補正を行うが、候補サイトレポートは `q_rate_enrichment_global`、すなわち解析的P値のBH補正を既定で使用する。さらに観測されたsupportで集合を絞り直してBHを再計算し、複数の閾値でそれぞれ候補を採る。

確認したコンテナのCSUBST実装では、解析的検定は推定イベント量と推定exposureをPoisson LRTに入れる。full-scan補正は置換ごとに探索を再実行する別経路である。したがって「full_scanを指定済みだから、後段の選択にもその補正が効いている」とは言えない。解析的P値の校正、supportフィルタ下の校正、複数閾値探索の全体誤検出率を区別する必要がある。

修正案：主解析用の検定集合・閾値・P値列を事前に固定する。経験的maxT列を用いる場合も、置換の交換可能性、置換失敗数、解像度、ファミリー間の多重性を検証する。支持数別出力は感度分析と明示し、最も有利な閾値の採用を確証的FDR制御と呼ばない。[選別と検定の独立性](https://pmc.ncbi.nlm.nih.gov/articles/PMC2906865/)。今回実データのFDR過大を実測したわけではない。

## 7. P2・既定値を確認：OUシフトはAICc選択・bootstrapなし

場所：[設定](/Users/kf/repos/genegalleon/workflow/gg_gene_evolution_entrypoint.sh:274)、[モデル実行](/Users/kf/repos/genegalleon/workflow/support/detect_OU_shift_kfl1ou.r:817)。ワークフローは `AICc`、`nbootstrap=0`、収斂レジーム推定ありを既定とする。シフト場所を多数探索した最良モデルの採択は、そのシフトの統計的確証や適応の証明ではない。

l1ou原論文は、シフト配置の組合せ数と系統相関を考慮するpBICを提案している。AIC系で多くのシフトを拾う懸念があるが、改変エンジンkfl1ouの現状に対する誤検出率は別途測定が必要。[原論文](https://doi.org/10.1111/2041-210x.12534)。なお現実装は共有measurement errorを推定し、replicate由来の誤差も渡す。「測定誤差を全く無視」は誤った指摘になる。それでも組織・batch・発現尺度の非比較性はそのままでは解消されない。[OUモデルと測定誤差](https://pubmed.ncbi.nlm.nih.gov/27478249/)。

修正案：pBIC等との感度比較、帰無シミュレーション、シフト選択を含むbootstrap、組織・条件の対応監査を用意し、図は「選択されたOUレジーム」と解釈する。

## 8. P2・確定した補正範囲：発現PGLSの補正はファミリー内

場所：[RSC summary](/Users/kf/repos/genegalleon/workflow/support/reconciled_speciation_contrast.py:1180)、[species PGLS summary](/Users/kf/repos/genegalleon/workflow/support/species_tree_pgls.py:1006)。Holm/BHの範囲は明示的に `all_usable_family_associations` 等。DBの全体BH追加はCSUBSTのaa_change用で、発現PGLSの全ファミリー横断補正はこの経路で確認できない。

1ファミリーを事前指定した検定としては問題とは限らない。数千ファミリーを走査して各summaryの調整P<0.05を並べる場合は、実験全体の多重性が残る。複数手法のうち有意なものを採用する場合も追加の選択がある。

修正案：全ファミリー×事前指定した応答×説明変数の検定表と、適切な全体補正を追加する。ファミリー単位のomnibusを対象にするなら、その仮説単位を明示する。既存のファミリー内補正を「何も補正していない」と扱ってはいけない。

## 9. P2・科学的解釈の制約：TimeTreeのCIを自動的に較正分布の境界にする

場所：[GeneGalleon呼出し](/Users/kf/repos/genegalleon/workflow/core/gg_genome_evolution_core.sh:3718)、[NWKIT変換](/Users/kf/repos/nwkit/nwkit/mcmctree.py:343)。TimeTreeの `precomputed_ci_low/high` がMCMCtreeの `B(lower, upper, tails...)` に変換される。soft boundsであり、ここをhard boundsと呼ぶのは正確ではない。

TimeTreeの集約区間と、対象研究の較正事前分布は意味が同一ではない。元研究の重複、化石較正・配列の共有、複数ノード間の相関があれば、独立の情報としての再利用は過剰な精度を生み得る。TimeTreeは研究間の時刻分布を用いる説明をしている。[公式FAQ](https://timetree.temple.edu/faqs)。二次較正の誤差についてはシミュレーション研究もある。[Schenk 2016](https://pubmed.ncbi.nlm.nih.gov/PMC4732660)。

修正案：自動値は較正候補と位置付け、元研究・対象ノード・較正根拠・データ重複を確認する。独立に正当化した較正、prior-only解析、較正を外す感度分析と比較する。二次較正を一律禁止すべきという結論ではない。

## 10. P2・確定した未考慮項目：CAFEとコピー数回帰に注釈・検出誤差が入らない

場所：[CAFE呼出し](/Users/kf/repos/genegalleon/workflow/core/gg_genome_evolution_core.sh:5650)、[copy-number regression](/Users/kf/repos/genegalleon/workflow/support/orthogroup_copy_number_trait_pgls.r:306)。CAFEはgammaカテゴリを指定するがerror modelを渡さず、コピー数回帰も観測カウントを確実な説明変数として使う。gammaによるファミリー間速度差と、観測誤差は別物である。

ゲノム品質・組織依存のトランスクリプトーム・isoformやhaplotig・断片注釈の違いでカウントが変わると、増減や形質関連に混入する。BUSCOの図や整数検証だけでは個別ファミリーの欠損観測を訂正できない。CAFEの公称0.05の星を集計した「significant」図も、全ファミリー×枝のFDR保証ではない。[CAFE公式説明](https://github.com/hahnlab/CAFE5)。

修正案：入力の観測単位と品質を揃え、CAFEの適切なerror model、品質層別の感度分析、コピー数の独立検証を行う。PGLSでは必要に応じ観測過程／交絡共変量を扱う。完全なゲノム注釈を前提にした既存PGLSの数式そのものが誤りという指摘ではない。CAFEの根での存在条件は依存側のフィルタも関係するので、今回独立した確定バグとは判定していない。

## 11. P2・確定した記述統計：GBIF由来の値は生物学的分布域の不偏推定ではない

場所：[分布指標生成](/Users/kf/repos/genegalleon/workflow/support/generate_species_trait.py:1318)、[取得上限](/Users/kf/repos/genegalleon/workflow/support/generate_species_trait.py:1454)。緯度極値・平均座標・占有グリッド面積等を取得済みpresence座標から直接計算する。ページ取得には上限があり、超過時の採用点は無作為／地域層別標本ではない。truncatedの記録はある。

観測密度の違いは面積・限界・重心を変える。既定1度グリッドの占有面積はIUCN等の標準AOOと同一ではない。地理エラーフラグ除外もサンプリング努力の差を直さない。[GBIFの観測・サンプリングデータの区別](https://techdocs.gbif.org/en/data-publishing/data-quality-recommendations)、[空間的偏りの実証研究](https://doi.org/10.1016/j.ecoinf.2013.11.002)。

修正案：「取得した観測記録の範囲」と表示し、標本数・切詰め・空間解像度・時期を必須メタデータにする。自然分布を問うなら移入・栽培等の扱いも定義し、努力量／空間偏りを考慮した推定と比較する。これらの値をそのまま遺伝子コピー数の説明対象とすると、研究努力量との相関を拾う可能性がある。

## 12. P2・追加検証が必要：RSC階層モデルと非Gaussian回帰の少数標本推論

場所：[RSC設定](/Users/kf/repos/genegalleon/workflow/gg_gene_evolution_entrypoint.sh:220)、[コピー数モデル説明](/Users/kf/repos/genegalleon/docs/copy-number-trait-models.md:1)。RSCは独自のreconciled contrastを階層モデルに入れ、既定はWald推論・最小species events=2。コピー数モデルはGaussianだけでなく、Laplace MLのbinomial／Poisson／NBを提供する。

同じ種分化に属する複数paralogを扱うための重み・階層化・replicate処理は既にあり、単純な擬似反復との断定は不適切。ただし最小件数の通過やoptimizer成功は、少数の独立イベント、偏った二値応答、分散境界、モデル誤指定下のP値校正を保証しない。新しい同時選択にはnested CVとbaselineがあるため、選択後の通常P値を出しているという指摘も当たらない。

必要な検証：独立イベント数・不均衡・重複数・共変量相関・測定誤差・欠損を振ったシミュレーションで、Type I error・被覆率・失敗率を確認する。RSCでは利用可能な分散パラメータ再推定付きbootstrapとの比較も行う。今回の調査はこれらの統計的妥当性を否定も確証もしていない。

## その他の機能と、誤って問題視しなかった点

| 機能 | 今回の評価・使用上の境界 |
|---|---|
| MAFFT、IQ-TREE、OrthoFinder、ASTRAL等 | 確立した方法を組み合わせている。ただしalignment、遺伝子木誤差、ILS、重複、サンプリングに対する頑健性は実データで別途評価が必要。各依存ソース全体の監査はしていない。 |
| GeneRax/HGT | DTLに基づく候補抽出として根拠がある。`score_hgt_candidates.py` は現在、主に証拠列を集約しており、未校正の合成スコアを勝手にHGT確率へ変換しているわけではない。HGT確証には別のroot／taxon sampling、支持、コンタミ・ゲノム文脈を確認する。[GeneRax原論文](https://pubmed.ncbi.nlm.nih.gov/32502238/) |
| コンタミ除去 | 分類学的不一致の除去はQCとして有用だが、外来由来の真の遺伝子も同じ条件を満たし得る。HGT用途では除外配列を追跡し、陰性を「HGTなし」と解釈しない。 |
| dN/dS、CodeML two-ratio、RELAX | 率推定や選択圧変化の解析として理解する。two-ratioの率差だけでは有意な正の選択とは言えず、RELAXの強化も常に正の選択を意味しない。今回はこの解釈を誤って自動確定している実装は確認していない。 |
| gene-family presence/absence | 「入力から検出できない」と「生物学的喪失」を区別する。UFBootをreconciliation事象の確率としない説明は既にあり、これは適切。 |
| self fractionation | 文書はextant annotated genesに条件付けたself-synteny retentionで、祖先喪失数ではないと明示している。これを根拠なく祖先喪失推定のバグとはしない。 |
| イントロン欠損 | unknownを観測0に変換しない実装があり、HGTのイントロン支持もimputed値と区別されている。この点は適切。 |
| コピー数PGLS・同時選択 | 前者にはglobal／per-trait BHがある。後者はnested group CV・学習fold内標準化を実施し、通常の選択後P値を出さない。この点は適切。 |

## 実行した検証と限界

GeneGalleon Dockerイメージ `local/genegalleon:trait-selection-pilot-20260910-precision` にリポジトリをread-only mountして実行した。

```sh
docker run --rm --entrypoint /bin/bash \
  -v /Users/kf/repos/genegalleon:/review:ro -w /review \
  local/genegalleon:trait-selection-pilot-20260910-precision -lc \
  'python -m pytest -q -p no:cacheprovider workflow/tests/test_asr_intron_evolution.py workflow/tests/test_plot_csubst_aa_change_summary.py workflow/tests/test_reconciled_speciation_contrast.py workflow/tests/test_generate_species_trait.py'

docker run --rm --entrypoint Rscript \
  -v /Users/kf/repos/genegalleon:/review:ro \
  local/genegalleon:trait-selection-pilot-20260910-precision \
  /review/docs/reviews/2026-09-10-scientific-validity/reproduce_go_filter.R /review
```

既存の選択した4テストファイルは **63 passed in 11.37s**。GOの最小例も前述の値を再現した。これらは局所的な実装検証であり、解析全体の科学的検証ではない。実データの全工程再実行、大規模帰無シミュレーション、全R／全Python suite、SIF実行、新規コンテナビルドは行っていない。今回のレビュー資料だけを追加したため、解析コード変更用の全検証は実施対象にしていない。

優先順は、(1) GO補正と検定単位、(2) MCMC証跡・収束判定、(3) 実験的年代と機構ラベルの適切な表示、(4) CSUBST・発現検定の全体多重性、(5) モデルごとの帰無校正・入力品質の感度分析。数値計算の正しさ、統計的な誤差率、生物学的な因果・機構の主張をそれぞれ別に検証することが必要である。
