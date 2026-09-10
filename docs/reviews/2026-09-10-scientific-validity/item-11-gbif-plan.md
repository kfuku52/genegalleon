# 項目11：GBIF観測範囲と分布推定を区別する対応計画

2026-09-10。調査と計画のみ。解析コード・設定の変更、コミット、push、他タスク作成、GBIF download申請は実施していない。

## 結論

元レビューの中心的な指摘は維持する。現在の値は「指定条件で取得できたGBIF presence記録の座標要約」であり、真の生物学的分布域の不偏推定とは位置付けられない。ただし、記述統計としての使用を禁止する根拠にはならず、GeneGalleonの実データで遺伝子関連解析の偽陽性が増えたことも今回は実証していない。

推奨は、①観測指標・品質情報の契約を明確にする、②取得完了状態と座標処理を直す、③再現可能な全件取得と感度分析を整える、④必要な研究だけ自然分布推定に進む、の順である。全件取得だけでもGBIFに未収録の分布は回復しない。標本数を揃えるだけでも観測努力の偏りは除去できない。

## 調査対象と後続修正の照合

- 元レビュー：[review.md](/Users/kf/repos/genegalleon/docs/reviews/2026-09-10-scientific-validity/review.md:97)。記載された基点は `15481838663118b1bb0ff4e456a0be46ab3c34eb`＋当時の未コミット変更。
- 調査時点では元チェックアウトと本作業ツリーのHEADはいずれも `021c5b62c551d073016ce14676a56ec9768f75c8`。本作業ツリーは調査開始時にclean。元チェックアウトにはレビュー以外にも多数の未コミット変更があるため、読み取り専用で参照した。
- 基点から現在HEADまで `generate_species_trait.py` の差分はない。元チェックアウトの未コミット差分は二値欠測の保持、strict時の全欠測列検出、形質別観測数statsの追加。GBIF取得・幾何・既定の形質一覧には差分がない。対応する追加テストも欠測処理が対象である。
- 元チェックアウトのコピー数解析にはNWKITへのアダプター化、応答分布の拡張、nested CVによる同時選択がある。古い作業ツリーだけで下流を判定せず、これらを含む現在の実装を確認した。GBIF品質による入力選別は確認できない。
- 上記は現在のファイルとGit基点の照合であり、元レビュー実施時の全未コミット状態を完全に復元したものではない。

## コードで確定したこと、既存対策、残る不確実性

主要な根拠は [集計・幾何](/Users/kf/repos/genegalleon/workflow/support/generate_species_trait.py:1236)、[取得・フィルタ](/Users/kf/repos/genegalleon/workflow/support/generate_species_trait.py:1398)、[キャッシュ](/Users/kf/repos/genegalleon/workflow/support/generate_species_trait.py:1504)、[既定形質一覧](/Users/kf/repos/genegalleon/workflow/support/generate_species_trait.py:137) にある。

| 論点 | 確認した現在の挙動 | 判断 |
|---|---|---|
| 座標品質 | APIでhasCoordinate、hasGeospatialIssue=false、PRESENTを指定。ローカルで座標範囲、任意のbasisOfRecord、不確実性閾値を扱う | 対策なしという批判は誤り。ただし不確実性欠測は閾値指定時も通過する |
| 分類群照合 | taxon key、match type、confidenceをキャッシュに記録。既定confidence閾値90 | 既存対策。confidence欠測は棄却されず、rankを厳密に要求する処理もないため高位分類群への照合などを試験する |
| 取得上限 | page size最大300、取得最大100,000。offset順に取得し、無作為抽出・地域層別抽出はしていない | 無作為性は保証されない。APIの実際の並びと偏りの方向・大きさは未検証 |
| 完了状態 | truncatedは初回countと上限だけで確定。途中の空ページ、短いページ、不正なpayloadでbreakしても更新しない | 中途終了が完全取得に見える経路がある。通信例外の全てが黙って成功するという意味ではない |
| 件数 | occurrence_countはAPI条件適用後、ローカル除外前。usedは採用座標数 | used/countは取得率と品質除外率を混同し、真の分布の網羅率にはならない |
| 時期・自然由来・努力 | 採用点をlat/lon/countryのタプルに縮約。eventDate、dataset、event、effort、establishment情報を集計に残さない | 事後に時期・栽培・努力量を精査できない。basisOfRecordの除外だけで野生・在来を保証できない |
| 重心 | 緯度の算術平均＋経度の円平均。記録を同じ重みで数える | 面積重心・個体数重心とは異なる。繰返し記録で動く |
| 経度幅 | 最大経度ギャップの補集合による最短区間。日付変更線の既存テストあり | 単純max−minではなく、既存対策を維持する |
| 占有面積 | 原点(-90,-180)の緯経度格子。占有セルを重複除去し、球面帯の式で面積加算。既定1度 | セル面積の緯度依存は考慮済み。ただし等面積格子、IUCN AOO、実生息地面積ではない |
| 凸包 | 平均緯度を使う正距円筒近似の平面上で凸包・面積計算 | 球面/測地線の凸包ではない。広域・極域での精度は未検証。IUCN EOOと同一視しない |
| provenance | coreには入力・パラメータのartifact契約が既にある。GBIFキャッシュも種集合と設定をhash化 | metadata皆無ではない。しかしレコードのsnapshot、取得時刻、時期・除外内訳・格子定義などの科学的証跡が不足 |
| 下流混入 | 件数・used・truncatedもspecies_trait.tsvに出る。コピー数解析のtrait=allはspecies以外を選ぶ | 数値化・変動条件を満たせば品質列も応答になる。生成もコピー数解析も既定offなので、全実行に自動混入するとの断定は誤り |

下流の根拠：[列選択](/Users/kf/repos/genegalleon/workflow/support/orthogroup_copy_number_trait_pgls.r:220)、[既定off・trait指定](/Users/kf/repos/genegalleon/workflow/gg_genome_evolution_entrypoint.sh:141)、[入力生成provenance](/Users/kf/repos/genegalleon/workflow/core/gg_input_generation_core.sh:1655)。コピー数解析にはglobal補正も既存実装にあり、多重補正が全くないわけではない。

**追加で確認した修正候補：** `gbif_max_distance_from_centroid_m` は、[説明](/Users/kf/repos/genegalleon/workflow/gg_input_generation_entrypoint.sh:155)では「種の観測重心からの距離」だが、実装はGBIFの `distanceFromCentroidInMeters` が閾値より大きいレコードを除く。このフィールドは地理参照に使われる既知の中心点からの距離であり、種の観測重心ではない。[GBIF形式仕様](https://techdocs.gbif.org/en/data-use/download-formats)。中心点への誤配置を除く意図なら近すぎる点を識別すべきで、現実装と説明は整合しない。既定無効なので全実行への影響ではない。GeneGalleon側のフィールド解釈・命名の修正対象とする。

コードから導かれる幾何の境界事例も試験対象にする。+180と−180が異なるセル番号になり得る、360を割り切らない格子幅で経度端セルの幅が切り詰められない、円平均の合成ベクトルがほぼゼロで経度重心が不安定になる、という経路がある。今回は実関数を動かした再現結果ではなく、静的に導いた試験仮説として扱う。

## 科学的な目的と出力の意味

### A. 推奨既定：観測記録の記述

推定対象は「分類群照合・地域・期間・品質条件と取得snapshotを固定した採用レコード集合」。記述そのものには帰無仮説やP値を付けない。

- 極値：採用記録の最北・最南緯度と観測経度の最短区間。
- 重心：record-weighted mean latitude / circular mean longitude。位置重複除去後の平均や占有セル重心を追加する場合は別指標とする。
- 面積：指定格子の少なくとも1採用記録を含むセルの全面積。未観測セルを不在としない。標本数0の分布指標はNAを保つ。
- 名前は `gbif_observed_*` など観測値と判るものに統一し、`distribution`、`limit`、`centroid`の図ラベルにも意味を添える。移行表で旧列と新列・単位・式・変更理由を対応付ける。

### B. 条件付きの別機能：自然分布・現在の占有

対象地域、期間、季節、在来/移入/栽培の扱い、解析単位を先に定義する。現在の自然分布を問うなら化石・過去記録を無制限に混ぜない。移入後の野生定着と栽培は別の状態であり、在来限定とも区別する。

presence-onlyの観測密度は、真の生物学的過程と観測・収録過程の双方で変わる。独立調査または観測過程に関する仮定なしに両者を一般に分離できない、という判断に立つ。GBIFはsampling-eventデータに方法・努力量等を記録する枠組みを持つ。[GBIF公式資料](https://techdocs.gbif.org/en/data-publishing/data-quality-recommendations)。presence-onlyと調査データの統合を行う場合の設計根拠は [Fithianらの原論文](https://www.stat.berkeley.edu/~wfithian/biasCorrection.pdf) とする。適合度だけで絶対的な在確率や真の面積を保証しない。

### C. 遺伝子との関連

記述指標を応答にする主解析のH0は「事前指定した観測範囲指標と遺伝子コピー数の回帰係数が、指定共変量・系統共分散の下で0」。生物学的分布との因果関係・適応の検定と読み替えない。

検証用の科学的H0は「潜在的な真の分布形質とコピー数の関連が0」。この状態で研究努力とゲノム品質が相関するシミュレーションを用い、観測値の関連を生物学的関連として誤解するリスクを測る。努力量補正後も未測定交絡は残り得る。観測件数は真の分布サイズの結果でもあり、単純なlog(count)調整を万能な補正にしない。

## 推奨案と代替案

| 選択肢 | 用途・判断 |
|---|---|
| 観測指標＋品質sidecar＋明示的な解析採用 | 最初に実装する推奨案。現在の軽量presetを維持できる |
| GBIF非同期downloadまたは既存downloadファイルの取込 | 論文用・上限超過・多種の本解析に推奨。query、download key/DOI、ファイルhashを記録。認証情報は出力しない |
| search APIの上限内取得 | 小規模確認・探索用。上限超過はcapped_partial。完全取得が必要な実行では明示的に停止し、部分値を解析用の正式出力にしない |
| 全取得後のseed固定層別抽出・同数化・空間間引き | 計算量制御と感度分析用。地域/時期/datasetなどの層と採用確率・重みを残す。最初のページだけを後から無作為化して全体標本と呼ばない |
| 努力量を扱うSDM/点過程、調査データ統合 | 自然分布推定を要する場合に追加。共通biasを仮定するtarget-group backgroundは候補だが、対象群とbias共有の妥当性を検証する。[Phillipsら2009](https://pubmed.ncbi.nlm.nih.gov/19323182/) |
| 2×2 km格子による占有セル集計 | IUCNとの比較が必要な場合の別出力。格子を2 kmに変えるだけで正式なAOO評価にはならない。IUCNは適切な分布・季節などの範囲と2×2 kmの尺度を要求する。[IUCN Guidelines §4.10](https://nc.iucnredlist.org/redlist/content/attachment_files/RedListGuidelines.pdf) |

GBIF公式APIは1ページ300、search全体100,000件を上限とし、超過は非同期downloadを案内している。これは依存側の欠陥ではなくサービス契約である。上限を定数だけ大きくする修正や、細分化queryで無計画に大量取得する案は採用しない。[Occurrence API](https://techdocs.gbif.org/en/openapi/v1/occurrence)。空間偏りが推定結果を変え得ることは元レビューが挙げた [Beckら2014](https://doi.org/10.1016/j.ecoinf.2013.11.002) と整合するが、その研究の効果量をGeneGalleonへ転用しない。

## 段階的実装と所有リポジトリ

### 1. 出力契約・解析への受け渡し

所有元はGeneGalleon。生成器の指標名・説明を揃え、数値の品質列を既定の形質表から分離する。提案するsidecarは種別品質TSVとrun/指標定義JSON。既存artifact provenanceから参照し、独立した未接続metadataにしない。

必須項目：要求名・採用taxon key/rank/accepted name・match根拠、APIとquery、取得開始/終了UTC、download識別子またはsearch snapshot hash、raw件数・取得件数・unique gbifID数・採用件数・除外理由別件数、取得上限・終了理由、期間/季節/地域、日付・不確実性・自然由来情報の欠測数、dataset/basis/eventの内訳、座標基準・格子原点/解像度/面積式・重み・重複規則、ソフトウェア/schema版。未知はunknownとして記録し、0やnativeに変換しない。

下流はmetadataのroleを用い、品質列をtrait=allに入れない。GBIF観測指標の関連解析は明示的な応答選択を要求し、capped/incomplete/未承認の分類照合を既定で不採用にする。除外種と理由は必ず出す。完全取得種だけに限定することで生じる種選択の偏りも感度分析で点検する。metadataなしの古いGBIFファイルを完全なデータと推定しない。

### 2. 取得状態・座標処理の修正

所有元はGeneGalleon。状態を `complete_search`、`complete_download`、`capped_partial`、`interrupted`、`invalid_response`、`taxon_unresolved` 等に分け、ゼロ件・全件品質除外も識別できる理由を残す。「complete」は当該queryの取得完了であり自然分布の完全性ではない。

raw取得数、重複ID、count変動、空/短いページ、終了フラグを照合する。ライブsearchを厳密なsnapshotと主張しない。completeを確認できない結果を成功キャッシュとして再使用しない。

centroidオプションは既知の地理中心点への近接を扱う名前・方向に修正する。種の観測重心からの外れ値削除は真の周縁分布を消すため既定導入しない。欠測時の扱いを明示する。GBIFのフィルタ例も中心点から一定距離以上の点を扱っている。[GBIF filtering guide](https://data-blog.gbif.org/post/gbif-filtering-guide/)。

経度を格子用には半開区間に正規化し、端セルをクリップするか格子の有効範囲を厳密に制約する。有限値・格子幅・対蹠点・極を扱う。凸包は検証した地理的適用範囲を明示し、適用外をNA＋理由にするか、測地学的定義を持つ別方式へ更新する。妥当性を確認せず全世界用へ昇格させない。

### 3. 再現可能な取得と観測過程の監査

所有元はGeneGalleonのGBIFアダプター。download取込、必要フィールドの保持、ID単位の同一レコード重複除去、更新/再利用の明示的選択を実装する。設定だけのcache keyに加えsnapshot・schema・処理版を照合する。データsnapshotの記録は上流プログラムの既定SHA固定とは異なる。

eventDateの期間/精度、dataset/eventに基づく重複候補、座標不確実性、化石/生体標本、establishment/degreeOfEstablishment等の採用可能な公式フィールドを検証して記録する。同一座標だけで別時点の正当な観測を削除しない。自然由来不明を在来とみなさない。栽培/移入メタデータは不完全なので、必要なら地域フロラ等の独立根拠を併用する。

### 4. 感度分析・下流の検証

所有元はGeneGalleon。同じsnapshotから、期間、地域、dataset除外、自然由来のunknown扱い、解像度、record/unique-location/cell重み、標本数を変えた結果を並べる。極値に加える分位点は別名の指標とし、都合のよい定義を結果を見て選ばない。

単純回帰に努力量・ゲノム品質共変量を渡す設計を項目10と調整する。NWKIT側の一般的な回帰・測定誤差・観測モデルの不足や欠陥はNWKITで実装/修正し、GeneGalleonに別の推定器や旧依存向けfallbackを重複実装しない。GBIFの原データ誤りはpublisher/GBIFへ報告し、由来を保存する。

### 5. 任意の自然分布推定

目的が確定し、観測努力・独立調査・環境データの可用性が確認できてから着手する。空間/時期ブロックによる外部検証、観測モデルの仮定、予測対象地域への外挿、面積への変換閾値、不確実性伝播を規定する。推定器の所有リポジトリは採用する依存プログラム、入出力・provenance・検証接続はGeneGalleonとする。今回は新しい依存プログラムの採用までは決めない。

変更時に一括で追従する対象は `generate_species_trait.py`、入力生成entrypoint/core、`gg_entrypoint_config_vars.sh`、形質テンプレート/検証/行列生成、コピー数PGLS・同時選択および形質を読む発現/祖先復元経路、図のラベル、対応テストとvalidation manifest、入力規約・設定・レシピ文書。全producer/consumerを再検索する。coreは既存ファイル内で順序を維持し、stageファイルへ分割しない。既存の二値欠測修正を保持する。

## 試験と受入基準

以下は実装時の試験計画であり、今回の実行結果ではない。

| 層 | 試験 | 受入基準 |
|---|---|---|
| 取得契約 | 模擬APIで0件、上限−1/一致/+1、page境界、途中空ページ、短ページ、count変動、不正payload、重複ID、HTTP失敗、cache再利用 | 全fixtureで件数と終了理由が期待値に一致。incompleteをcompleteとしない。無観測の生物指標はNA。品質除外と取得不足を区別 |
| フィルタ | taxonの高位照合・同名・confidence欠測、日付範囲/欠測、座標不確実性欠測、中心点距離の閾値±ε、栽培/移入/不明、同一eventと別event | 指定policy通り。全除外に理由があり、unknownが既知状態へ化けない |
| 幾何 | ±180同一点、北南極、単一点、同一直線、対蹠点、全球、端数格子、NaN/Inf/非正幅、広域凸包 | 同一地点の重複で占有面積は変化しない。球面格子面積は独立解析解に相対誤差1e-10以内を目安とする。定義不能な平均はNA。凸包の新手法は同一定義の独立実装と照合し、適用範囲を受入条件に含める |
| 下流 | profile→表/sidecar→行列→PGLS/同時選択→図/provenanceの最小統合試験 | trait=allで品質列を検定しない。partialは既定で解析に入らない。明示的に採用した観測指標の意味・除外理由が最終結果まで残る |
| 再現性 | 同一snapshot・設定・seedの反復、schema/期間/格子/フィルタ変更 | 正規化した数値・採用ID集合が一致。変動する実行時刻等は別枠。異なる意味のcacheは再利用されない |

シミュレーションは有限の既知セル集合を真値として、島嶼・日付変更線・極域・分断分布を含める。観測努力を均一/都市・道路集中/地域別ゼロ/年代変化/dataset集中に変え、座標誤り・栽培混入・重複・不確実性欠測・ページの地域順/年代順/無作為順を交差させる。全潜在セル、完全GBIF相当標本、ページ先頭の部分標本を区別して比較する。

評価は面積のbias/RMSE、境界誤差、重心の測地距離誤差、記録数別の飽和曲線、種順位の変化。欠損領域を同数化だけで回復できない負の対照を含める。記述指標の実装合格と、真の分布推定の精度合格を分離する。科学的に許される面積・位置誤差は研究目的で決めるが、目標未達の方法を「不偏」や「補正済み」と表示しない。

遺伝子関連試験では系統相関のある種別分布形質とコピー数をH0下で生成し、研究努力と注釈品質の共有交絡、種の欠測、上限超過を加える。raw観測、完全取得、同数化、努力調整、潜在真値のoracleを比較する。帰無条件ごとに最低2,000反復を初期計画とし、境界では追加する。名目0.05に対してType I errorの片側95%二項信頼上限が0.065以下を暫定基準とする。95%区間被覆率、失敗/NA率、全ファミリー×事前指定応答でのFDRと検出力も別に記録する。成功fitだけに条件付けて合格にしない。観測指標そのものとのH0と、生物学的H0を取り違えない。

既存テストは [GBIF模擬取得](/Users/kf/repos/genegalleon/workflow/tests/test_generate_species_trait.py:797) と [日付変更線](/Users/kf/repos/genegalleon/workflow/tests/test_generate_species_trait.py:922) を保持・拡張する。前者は小標本の値と正の面積を確認するもので、面積精度・上限超過・偏り校正の証拠ではない。

コンテナ検証は [runtime validation](/Users/kf/.codex/worktrees/1098/genegalleon/docs/agent-runtime-validation.md) に従い、変更した現行コードから作ったGeneGalleon DockerでPython/Rの局所試験と統合試験、リポジトリ必須検証を実行する。その後Linux/HPCのGeneGalleon SIFで同じ保存fixtureと最小end-to-endを確認する。Docker結果をSIF互換性と呼ばない。依存プログラムの既定branchを固定tag/SHAに書き換えない。テスト対象snapshotや実行時の依存解決結果は証跡として残す。

## 計算量と予算の見積もり

種iのAPI countをN_i、取得上限をC、page sizeをP≤300とすると、searchの呼出数は概ね Σ_i[2+ceil(min(N_i,C)/P)]（照合1＋count1＋ページ、再試行を除く）。100,000件の種は約336呼出となる。少数種のpilotでpayload bytes、応答時間の中央値/上位分位、再試行率を測り、通信時間・転送量を見積もる。多数種は一括downloadを比較し、サービス負荷・待ち時間を別計上する。[GBIF API利用案内](https://techdocs.gbif.org/en/openapi/)。

N件のフィルタ・格子集計は概ねO(N)、占有セルC_gの保存はO(C_g)。現在の経度区間と凸包はソートを含みO(N log N)、点の保持はO(N)。download保存はO(N×1レコードのbytes)。raw保存・中間・索引・感度分析結果を別々に測る。大規模化時はchunk処理を検討するが、全件の凸包等のメモリ必要量を無視しない。

感度条件K、反復B、ファミリーF、種Sで、回帰fit数は概ねK×B×F×応答数。共分散の密行列分解は1回O(S³)が目安だが、再利用可否や最適化回数は依存実装で測る。小/中規模pilotのwall timeとpeak RSSから予算を外挿し、nested CVはfold数・候補数を追加乗算する。同じ出力・同じ入力のbefore/afterを取るまで高速化率を主張しない。

## 他項目との依存関係

- 項目10：研究努力・ゲノム/トランスクリプトーム品質が双方の観測に影響するため、コピー数誤差・品質共変量の契約を共同で設計する。GBIF単独の修正では終わらない。
- 項目12：少数種、非Gaussian回帰、測定誤差、nested CVの校正と共通のシミュレーション基盤を使う。交絡した観測指標を予測できても自然分布の予測成功とは言えない。
- 項目8：検定集合と応答選択を事前固定し、関連を走査する全ファミリー/応答/手法の補正範囲を明示する。コピー数側の既存global補正は維持する。
- 項目3・9：系統共分散に年代付き樹を使う場合は、年代・較正の不確実性への感度を含める。段階1〜3の取得/metadata修正はこれらを待たずに進められる。
- 既存の形質欠測修正：unresolved・zero records・unknownの意味を二値0へ潰さない契約を共有する。

## ユーザー判断が必要になる点

計画完成のための追加回答は不要。段階1〜3の意味明確化・取得状態・品質列分離・境界修正は証拠から推奨できる。

自然分布推定を実装する段階でのみ、研究上の対象（在来分布、移入後の野生分布、栽培を含む観測分布のどれか）、対象地域/期間/季節、科学的に許容できる位置・面積誤差を確定する必要がある。新規downloadの申請には使用するGBIFアカウントと取得予算が必要になるが、既存downloadファイルの取込は先に用意できる。

今回実施した検証はコード・差分・既存テスト内容・公式資料/原論文の読み取り照合のみ。解析実行、シミュレーション、実GBIFレコードの取得検証、Docker/SIF検証は未実施であり、科学的校正や実行環境互換性を確認済みとはしない。
