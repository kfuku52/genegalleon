# CSUBST 3di20対応の実装計画

2026-09-11。承認された実装計画。実装内容と検証結果は
[3Di利用ガイド](csubst-3di.md)と[検証記録](reviews/upstream-compatibility-2026-09-11/3di/README.md)を参照。
先行修正のsitesへの遺伝暗号転送は実装・Docker検証済み。

## 到達目標と範囲

`csubst_nonsyn_recode=3di20`を選ぶと、GeneGalleonのsearch、scan、通常sites、
scan候補sitesが、同じfull CDS alignmentと対応する推論結果を用いて完了する。
TSVだけでなく、PDF・系統樹パネル・候補一覧のbranch/site対応まで確認する。

- scan、shift、regressの既存の統計設定は変更しない。
- 3Di backend、ASR方式、3Di IQ-TREEモデル、推論バッチサイズは上流既定を使う。
  現在の上流既定は`esm3di-35m`、`direct`、`GTR`。値をGeneGalleonの固定設定へコピーしない。
- 既存の`no`、Dayhoff等のアミノ酸recodingを維持する。
- タンパク質入力だけの解析は、full CDSを要求する本対応の対象外。前提不足は早期に明示する。
- 既存結果の一括無効化や、推論方式変更に伴うキャッシュ移行（旧指摘5）は追加しない。
- coreは`workflow/core/gg_gene_evolution_core.sh`内で編集し、stageファイルに分割しない。

## 確認済みの不足

1. search・scan・sitesはいずれも`--alignment_file`を渡し、3Di必須の
   `--full_cds_alignment_file`を渡していない。
2. 現在の`iqtree_anc` bundleはトリミング後のcodon alignmentから作られる。
   full alignmentに対応するfitとしてそのまま転用できない。
3. コンテナで`csubst[3di]`相当の依存と共有モデルキャッシュを保証していない。
4. `prepare_recoded_site_alignment`は3di20の状態パネルをスキップする。
   既存レポートのsite番号・trimmed/full alignmentの対応も監査が必要。

## 1. 小規模な上流直接実行で契約を確定する

最初に6〜8 tip程度の固定CDS fixtureを用意し、GeneGalleon Docker内で上流の
search → sites、scan → sitesを直接実行する。full/trimmedで列数が違うデータにする。

- 3di20では`--full_cds_alignment_file`のみを入力指定し、`--alignment_file`を付けない。
- full CDSの翻訳、gap除去した各配列への3Di推論、アラインメント位置への復元、
  direct 3Di ASRはCSUBST自身に実行させる。
- full codon fitとdirect 3Di fitの用途・必要ファイルを区別する。
  3Di用MORPH/GTRのIQ-TREE呼び出しをGeneGalleonに再実装しない。
- CSUBSTが出力する3Di alignment、状態・site対応、再利用可能なcache/manifestを確認する。
  一時ファイル名だけに依存するインターフェースを作らない。
- 必要な出力が上流に公開されていなければ、所有するCSUBST側の小さい修正として切り出す。
  非公開ファイルの推測やダミー状態で埋めない。

完了条件：2つの直接実行経路が実backendで成功し、後続sitesに渡すファイルと
site番号の定義が文書化される。ここで得た実出力を以後の比較基準にする。

## 2. full CDSとfitをfamily単位で用意する

主な変更：`workflow/core/gg_gene_evolution_core.sh`。

- `file_og_untrimmed_aln_analysis`をfull CDSの候補入力にする。
  MAFFT後・列トリミング前であり、解析途中で配列端などが失われていないことも確認する。
- 解析に採用したrooted treeのtip集合へ配列を絞る。
  一意なID、tipとの完全対応、同一alignment長、codon frameを検証する。
- 複数塩基gap・曖昧codon・終止codonの扱いは上流仕様に従う。
  欠けた配列や不正frameを無言で補正しない。
- 3di20時だけ、同じrooted topologyと遺伝暗号でfull alignment用codon fitを用意する。
  既存のtrimmed用fitと物理的に分離する。
- 既存`iqtree_anc` archive内に3Di専用サブディレクトリを追加し、full CDSと必要fitを保存する。
  `csubst.input.json`を版管理して拡張し、相対パス・遺伝暗号・入力識別情報を記録する。
  現行schemaの通常解析bundleも読み込めるようにする。
- 新しい3Di入力は既存のartifact provenance契約に登録し、異なるfamilyや
  異なる入力のfitを誤って再利用しない。これは新経路の入力管理であり旧結果の一括移行ではない。

完了条件：bundleだけからfull CDS・対応fit・遺伝暗号を復元でき、
通常解析用のtrimmed入力とは取り違えられない。

## 3. 4つの呼び出し経路を接続する

| 経路 | 主な変更ファイル | 作業 |
| --- | --- | --- |
| search | `workflow/core/gg_gene_evolution_core.sh` | 3Di時にfull CDSとfull codon fitへ入力を切替 |
| scan | 同上 | searchと同じ入力契約を使う |
| 通常sites | `workflow/support/csubst_site_wrapper.py` | bundleの3Di入力を読み共通コマンドビルダーへ渡す |
| scan候補sites | `workflow/support/csubst_scan_candidate_sites.py` | 候補元と同じ3Di bundle・branch IDsを使う |

共通の3Diコマンド形は次のとおり。実際のfitオプションは手順1の確認結果に合わせる。

```bash
csubst <search|scan|sites> \
  --nonsyn_recode 3di20 \
  --full_cds_alignment_file <bundle内のfull CDS> \
  --rooted_tree_file <対応するrooted tree> \
  --genetic_code <familyの遺伝暗号> \
  <対応するfull codon fit・各サブコマンド固有の引数>
```

`--sa_backend`・`--sa_asr_mode`を固定せず、採用された実値は上流の実行情報に残す。
CSUBSTの再利用可能な3Di状態cacheは同一入力の後続処理へ接続する。
通常解析用bundleしかない既存familyで3di20を選んだ場合は、必要なfull入力から
3Di bundleを作る。full入力がない場合は取得・再生成すべき段階を示して停止する。

## 4. コンテナと共有キャッシュを接続する

対象は`container/Dockerfile`、native Apptainerのビルド経路、両者が読む共通依存定義、
`container/scripts/validate_runtime.sh`およびworkflowのruntime環境設定。
実装開始時に共通定義の所在を確定し、片方だけに依存を追加しない。

- 現在の上流extraは`huggingface-hub`、`sentencepiece`、`protobuf`、`torch>=2.6`、
  `transformers`、`peft`。既存のConda/Pip分離を守って不足だけを追加する。
- source wheelの`--no-deps`を外して環境全体を再解決する変更は行わない。
- CSUBSTリソースは`workspace/downloads/csubst`、Hugging Faceリソースは共有downloads内の
  専用ディレクトリへ接続する。既存の明示cache環境設定を尊重する。
- 4経路すべてで同じ環境を使い、モデル識別・検証・排他制御は上流に任せる。
  一時作業ディレクトリへモデル重みを毎回ダウンロードしない。
- ビルド時は依存importと必要APIを検証し、大きなモデル重みの取得は実行時に分離する。
- 初回取得、再利用、事前取得後のオフライン実行、複数familyからの並列取得を検証する。
- 代表的な短い配列でCPU実行の時間・ピークメモリを記録する。
  長配列を黙って切り詰めず、GPU必須とも決め付けない。

## 5. site表示と下流出力を対応させる

主な対象：`csubst_site_wrapper.py`、`csubst_scan_candidate_sites.py`、必要なら
`plot_csubst_aa_change_summary.py`・`generate_orthogroup_database.py`・関連するtreevisパネル。

- full alignment上のsite番号をtrimmed alignmentの同じ番号として使用しない。
  full、trimmed、gap除去後の配列位置を区別する。
- 3Di状態パネルはCSUBSTが実際に使った3Di alignmentから作る。
  20-aa配列から通常のrecoding表で擬似変換しない。
- 元のアミノ酸/CDS表示を併記する場合は対応位置を確認し、3Di状態と明確に区別する。
- scan候補一覧・DB・図の列名や凡例に、3Di状態を通常アミノ酸と誤解させる箇所がないか確認する。
- full fitのbranch IDsとGeneGalleon stat.branchの対応を既存のidentity検証で確認する。
- 既存の出力名・archive構造は維持できる範囲で維持し、変更が必要なら
  producers、consumers、tests、ドキュメントを同時に更新する。

## 検証と完了条件

| 検証 | 合格条件 |
| --- | --- |
| 入力・コマンドテスト | 4経路がfull CDSを使い、trimmed用`--alignment_file`やfitを混用しない |
| full/trimmedで異なる列数 | TSVのsiteとレポートで強調する残基・3Di状態が一致 |
| gap・曖昧codon・除去tip | 上流と同じ扱いになり、不正入力は説明付きで拒否 |
| 遺伝暗号1・2 | familyごとの値が転送される。上流に制限があれば上流側で解決 |
| search → sites | 上流直接実行とbranch/site IDs・採用候補・数値結果を比較し、一致または説明可能な数値許容差内 |
| scan → candidate sites | 同じ候補・支持branch・状態を再現し、TSV/PDF/manifestが揃う |
| cache初回・再利用・並列・offline | 再利用時に重みを再取得せず、並列でも破損しない |
| 通常のno/Dayhoff | 既存search/scan/sitesの回帰テストが通る |
| コンテナ | 現行ソースからのDockerで実backendまで成功。利用可能なLinux環境ではSIFも別途検証 |

固定された小さい推論出力を使う通常テストと、ネットワーク・モデル取得を伴う
実backendテストを分離する。stubだけで対応完了とはしない。
実施できなかったplatformや実backend試験があれば、完了報告で明示する。

## 実装の区切り

1. 上流直接実行とfixture・必要ファイルの確定。
2. family bundle・search/scan/sites入力の接続。
3. コンテナ依存・cacheと3Diレポートの接続。
4. 実backendの一連の実行、通常経路の回帰確認、使用手順の更新。

この計画の承認後に実装・Docker実backend検証・モデル取得を実施。
その後の公開準備と追加検証は[公開検証記録](reviews/upstream-compatibility-2026-09-11/publication/README.md)を参照。
