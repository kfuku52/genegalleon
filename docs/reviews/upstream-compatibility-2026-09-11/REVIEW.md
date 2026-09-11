# 最新上流ツールとの互換性・コンテナ依存レビュー（2026-09-11）

> 後続対応: ユーザー指定に従い、指摘1はCDSKIT側で修正、指摘2はGeneGalleonの旧rate設定を削除した。指摘3・4の[実装プラン](../../csubst-integration-plan.md)を作成し、指摘5は対応しない。以下は修正前のレビュー記録。

**要対応は5件。既定で有効なCDSKIT局在予測に実行エラーがあり、CSUBST scanも有効化すると現行の既定設定で失敗する。** 通常のCSUBST search、NWKITの移行済み機能、CDSKITのAMAS互換統計、AMALGKITのFASTQから発現量集約までの接続は、下記のコンテナ検証で動作した。

製品コード・コンテナ構成は変更せず、レビューと再現証跡のみ追加した。指摘3・4は今回の照合で確認した現存問題であり、最近の更新によって初めて導入されたとは断定しない。

## 対象と検証方法

GeneGalleonの対象は `main` の `10fc03bd9c40c51cc635c5254e6b099c5ad75323`。`origin/main` より先行する23個のローカルコミットを含む。開始時の未追跡 `Rplots.pdf` と `workspace/.gg_cache/` は変更していない。

上流の `origin/master` をfetchし、以下のコミットを `git archive` から取り出してGeneGalleonのLinux ARM64 Docker環境内でビルド・インストールした。ここに記録するSHAは検証用の証跡であり、コンテナのソース既定値を固定するものではない。

| ツール | バージョン | 検証コミット |
| --- | --- | --- |
| CSUBST | 1.16.0 | `cde05f34e5ea51c3f0f1bae721ccd45243e46f38` |
| NWKIT | 0.43.16 | `1996588e561d10aebc47c4daf90be836bbebd315` |
| CDSKIT | 0.31.0 | `20328e15c544b01481e29b998956127416480467` |
| AMALGKIT | 0.16.88 | `21ce1417410ec05aea409370d496ce01896671d0` |

ベースは `local/genegalleon:nwkit-shift-covariance-dev`。NumPy 1.26.4、scikit-learn 1.9.0、PyTorch 2.13.0、ETE 4.4.0。上流wheelは標準のsource-artifactビルドと同じく `--no-deps --no-build-isolation` で生成した。既存イメージに不足していたビルド用コンパイラと、現行R依存宣言に既に含まれる `posterior` を検証環境へ追加した。通常の `pip check` は成功した。

`kftools`、`rkftools`、`kfFractBias` もfetchして照合し、最新コミットがベースイメージの記録と一致することを確認した。正確な環境情報は [manifest.json](evidence/manifest.json) に記録した。

## 1. [P1] 既定の局在モデルを現行scikit-learnでロードできない

場所: [モデル既定値](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/gg_gene_evolution_entrypoint.sh#L180)、[共通実行ヘルパー](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/support/gg_util/04_busco_runtime.sh#L491)、[コンテナのscikit-learn指定](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/container/env/base.required.txt#L14)。

`run_cdskit_localize=1` は既定で有効。GeneGalleonが指定する `targeting5-perox-deeploc21-et-v1` を、最新版CDSKITとコンテナのscikit-learn 1.9.0で読むと失敗する。実際の `gg_run_cdskit_localize` を合成タンパク質2配列で実行し、以下を再現した。

```text
Can't get attribute '__pyx_unpickle_CyHalfBinomialLoss'
on <module 'sklearn._loss._loss' ...>
```

モデルのSHA-256は `d0998df8819d975b4392342ab78dccc0dd95cf301e4d2df8f38c73d0b5aab445` で、現在のCDSKITレジストリと一致する。キャッシュ破損ではない。モデルに含まれるscikit-learn 1.5.2の推定器と、現在のランタイムとの互換性問題である。[実行ログ](evidence/localize.log)

**対応:** 修正の所有者はCDSKITの公開モデル／ローダー。GeneGalleonにpickle互換シムを追加せず、同じモデルの可搬化または上流修正を行う必要がある。上流の新しいESM2既定モデルへ切り替える場合は、予測対象・出力列とML依存も変わるため、現在の `p_noTP/p_SP/p_mTP/p_cTP/p_lTP` やperoxisome列を消費する集約・可視化まで別途検証する。モデル名だけの置換では完了しない。

標準コンテナの検証は現在、`cdskit` のコマンド存在を主に確認するため、このロード失敗を検出できない。使用するモデルそのもののロード・小入力予測を検証対象に加える必要がある。

## 2. [P1] CSUBST scanの既定設定がjoint推論と矛盾する

場所: [既定値](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/gg_gene_evolution_entrypoint.sh#L305)、[exposureの許容値](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/core/gg_gene_evolution_core.sh#L1442)、[実際のCLI](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/core/gg_gene_evolution_core.sh#L6090)。

CSUBST 1.16.0ではjoint endpoint posteriorが既定になった。GeneGalleonは推論方式を指定せず、旧marginal用の `--scan_rate_length n_rescaled --scan_rate_exposure q_weighted` を明示する。その結果、既存の実coreコマンドを実行する統合テストが終了コード2で失敗した。

```text
Joint/bridge observations require --scan_rate_exposure endpoint.
```

さらにGeneGalleon自身の設定検証が `endpoint` を許さないため、exposure値だけをユーザー設定で変更しても先へ進めない。`run_csubst_scan=0` が既定なので、scanを有効にするジョブが対象。[失敗ログ](evidence/runtime.log)

**対応:** 採用する推論方式を明示し、jointなら `endpoint/raw/posterior_sum` に設定・検証・記録を揃える。旧marginalを明示的に選ぶ運用を残すなら、その組み合わせとして扱う。無言のmarginalフォールバックは避ける。

対照実験では、GeneGalleonと同じ `ECMK07+F+R4` のIQ-TREE出力に `endpoint/raw` を指定するとscanが成功し、5候補、2foreground unitsとPDFを生成した。[成功ログ](evidence/scan-joint.log) 上流の契約は [ENDPOINT_POSTERIORS.md](https://github.com/kfuku52/csubst/blob/cde05f34e5ea51c3f0f1bae721ccd45243e46f38/docs/ENDPOINT_POSTERIORS.md) と [SCAN_CTMC.md](https://github.com/kfuku52/csubst/blob/cde05f34e5ea51c3f0f1bae721ccd45243e46f38/docs/SCAN_CTMC.md) を参照。

## 3. [P2] 公開設定の3di20に必要な入力と依存が渡らない

場所: [公開設定](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/gg_gene_evolution_entrypoint.sh#L299)、[search入力](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/core/gg_gene_evolution_core.sh#L5927)、[sitesラッパー](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/support/csubst_site_wrapper.py#L346)、[wheel導入](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/container/scripts/install_source_artifacts.sh#L23)。

GeneGalleonは `3di20` を有効な設定として案内するが、search/scan/sitesへ従来の `--alignment_file` とIQ-TREE bundleを渡し、`--full_cds_alignment_file` を渡していない。最新版で実行すると `--nonsyn_recode 3di20 requires --full_cds_alignment_file.` で失敗する。[再現ログ](evidence/3di.log)

入力を対照実験で修正すると、次に `3Di model downloads require huggingface-hub. Install csubst[3di].` で失敗した。現在のコンテナには `huggingface-hub`、`transformers`、`peft`、`sentencepiece` がなく、source wheelの `--no-deps` 導入ではextrasは増えない。[再現ログ](evidence/3di-full.log)

**対応:** 3Di用のfull CDS alignment、そこからのfit／予測、キャッシュ・provenanceを一貫して接続し、使用する上流backendのextrasをコンテナに明示する。trimmed alignmentで作った既存ASR bundleにfull alignmentの引数だけを付ける修正では、入力の対応が保証されない。通常の `no` やDayhoff等のrecodingまで3Di依存が必要になるわけではない。

## 4. [P2] 非標準遺伝暗号がCSUBST sitesへ伝わらない

場所: [sitesコマンド構築](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/support/csubst_site_wrapper.py#L334)。

gene evolutionのsearch/scanは `--genetic_code` を渡すが、sitesラッパーは引数も転送先も持たない。そのため遺伝暗号2で正常に作ったIQ-TREE bundleを後段で扱うと、sitesだけが遺伝暗号1を用いる。

合成8配列を `GY+F+R4 --seqtype CODON2` で解析し、実際の `build_csubst_sites_command` を使って比較した。

| 実行 | 結果 |
| --- | --- |
| 現在のラッパー | 終了コード2。codon tableとASRの遺伝暗号2の不一致で失敗 |
| 同じコマンドに `--genetic_code 2` を追加 | 終了コード0。sitesのTSV・PDF・output manifestを生成 |

[現行の失敗ログ](evidence/sites-code2-current.log)、[対照の成功ログ](evidence/sites-code2-explicit-code2.log)、[再現スクリプト](evidence/csubst_code2.py)。

**対応:** family単位のASR生成時の遺伝暗号をbundleのprovenanceから受け取り、sitesおよびsite再利用条件まで転送する。異なる遺伝暗号のfamilyをまとめる処理では、単一のグローバル値で代用しない。

## 5. [P2] CSUBSTの推論方式変更が結果の再利用条件に反映されない

場所: [search provenance](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/core/gg_gene_evolution_core.sh#L5830)、[versionの記録方法](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/core/gg_gene_evolution_core.sh#L5980)、[比較対象](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/support/artifact_provenance.py#L668)、[sites ZIPの再利用](https://github.com/kfuku52/genegalleon/blob/10fc03bd9c40c51cc635c5254e6b099c5ad75323/workflow/support/csubst_site_wrapper.py#L1917)。

新しいjoint推論は旧marginal推論と数値的に等価な高速化ではなく、置換数・候補順位・root近傍の適格性が変わる。ところがsearchの出力パラメータに推論方式がなく、`csubst_version` は比較対象外のdiagnosticsにのみ入る。IQ-TREE bundleと明示済みパラメータが同じなら、古いmarginal結果は再利用可能のままで、新規familyはjointで生成される。sitesも旧形式と同じarchive完了マーカーを用い、推論方式を確認しない。

**対応:** `substitution_posterior`／scan observationと、必要な出力意味論の契約を明示・provenance化し、既存出力を移行または明示的に再生成する。単に全依存のversionを再利用キーへ足すという話ではなく、今回の推論意味論の変更を区別する必要がある。searchで生成される `csubst_endpoint_model.json` も現在は管理出力へ保存されないため、fit／推論の監査記録として保持する。

## コンテナ内の不要プログラムの棚卸し

現行のDocker・native Apptainer構成、Conda/APT指定、source artifact、CRAN導入、required-command一覧、ワークフローの実呼び出しを照合した。**今回の対象範囲では、さらに削除できると確認できた直接インストールはなかった。** 移行に伴う主要な削除は対象mainに既に含まれている。

| ツール／依存 | 現行の扱い・残す理由 |
| --- | --- |
| AMAS | 削除済み。`cdskit stats --mode alignment` を使用。DNA/AA両方の生成・読込・既存出力移行テストが成功 |
| NOTUNG | 削除済み。NWKITによるroot候補とLCA loss/reconciliationへ移行済み |
| kfl1ou / l1ou | 明示インストールから除去済み。NWKIT native OUを使用 |
| 単独R RADTE / Rphylopars | 現行インストーラに要求なし。NWKITへ移行済み |
| GeneRax / GRAMPA | tree optimization・polyploidy解析の実呼び出しが残る |
| mapnh / HyPhy / codeml / CAFE5 | 別個の解析として現役。今回の内製置換だけでは不要にならない |
| PAML source / Conda PAML | source側はMCMCtree、Conda側はcodemlを供給。NWKITのgene-tree datingがspecies-tree MCMCtreeを置換したわけではない |
| MAFFT / ClipKIT / trimAl | MSA／選択可能なtrimming経路が現役。CDSKITのstats/backalign等とは役割が異なる |
| fastp / SeqKit / SRA-tools / Kallisto / oarfish / MMseqs2 | GeneGalleonとAMALGKITの取得・filter・quant経路に必要。fastpは今回の実FASTQ検証でも実行された |
| Rのape・phytools・phangorn・igraph・nlme・Rcpp等 | 残存するR解析・可視化、rkftoolsとRパッケージ間の依存がある |
| Conda IQ-TREE | 明示要求は削除済みだがOrthoFinderの間接依存として残る。実行入口は公式ソースからのIQ-TREE 3／worker。強制削除はしない |

過去に作成した個々のイメージに旧パッケージが残っていることと、現在のビルド構成がそれを要求することは区別した。既存の [依存監査](../../container-dependency-audit.md) には過去のkfl1ou検証の記述が残っているが、上部に移行済みとの注記があり、現行インストール要求ではない。

## 実施したチェックと結果

| チェック | 結果 |
| --- | --- |
| 最新4ツールのsource wheel生成・導入、`pip check` | 成功 |
| 抽出可能なshellの直書きCLI | 80箇所・32種類のtool/subcommandで長形式オプション名が一致 |
| support Pythonの上流import | 41箇所でmodule・symbol解決成功 |
| 静的テストlane | 397 passed |
| runtime Python lane | 初回247 passed / 5 failed。4件は既存イメージのR posterior欠落で、現行manifestに合わせて導入後に解消。残る1件は指摘2のscan |
| posterior関連の再検証 | 37 passed。上記4件を含む |
| transcriptome・CSUBST wrapper等の追加チェック | 179 passed / 1 failed。失敗はmacOS共有マウント上のflock挙動で、コンテナ内部の通常filesystemで同一テストを再実行すると成功。上流CLI変更とは別の環境制約 |
| NWKIT | runtime laneのroot/reconciliation、ASR intron、native OU/covariance、RADTE、regression/selection等が成功 |
| CSUBST search | 実coreのコマンドと既定の `ECMK07+F+R4` で成功。b/cb/cb_statsとendpoint modelを生成 |
| CSUBST scanの対照実験 | joint + endpoint/rawで成功しPDFまで生成 |
| CDSKIT localize | 実共通ヘルパー＋現在の既定モデルで指摘1を再現 |
| AMALGKIT | 合成single-end FASTQ 2,000 readsでgetfastq → 実coreのcompletion binding → quant → GeneGalleonのquant検査 → mergeが成功 |
| 非標準遺伝暗号sites | 指摘4の失敗／対照成功を確認 |
| 3Di | 指摘3の入力不足と、入力修正後の依存不足を確認 |

CLI抽出は直書きのshellコマンドを対象とし、動的配列の全経路を自動列挙したものではない。引数名が存在するだけでは値・出力意味論の互換性を保証できないため、上の実行検証と併用した。テスト集合間に重複があるので単純合計はしない。

AMALGKITの公開NCBI/GSAダウンロード、大規模MLモデル推論、全生物データでの数値等価性、native amd64、SIF実行は未検証。ローカルSIFファイルは存在するが、このmacOSホストにApptainer/Singularity実行環境はない。今回は既存GeneGalleon Docker環境への最新上流overlay検証であり、標準イメージ全体のクリーン再ビルド成功やSIF互換性を示すものではない。

初期セットアップで通常のpip build isolationを一度試した際のコンパイラ不足、およびその途中生成CファイルのNumPy不一致は、clean archiveから標準と同じビルド方法で作り直して解消した。製品のビルド欠陥として数えていない。AMALGKIT integrateはtaxonomyダウンロードに進んだため、その確認を中止し、FASTQ→mergeの実行検証には明示的な合成metadataを用いた。

## 証跡と再実行

[evidence/](evidence/) に実行ログ、JUnit XML、CLI/import照合JSON、補助スクリプトを保存した。スクリプトは検証対象repoを `/review`、証跡ディレクトリを `/audit` にmountしたGeneGalleonコンテナ向け。`amalgkit_flow.py`、`csubst_search_flow.py`、`csubst_code2.py` は `/tmp/gg-*` に合成入力を作成する。localize再現には上記ハッシュの公開モデルが必要。

runtime laneの再現コマンド:

```bash
python -m pytest workflow/tests -q -n 3 --gg-suite runtime \
  -p no:cacheprovider --basetemp=/tmp/gg-review-runtime
```

修正時は、指摘1の上流モデル互換性と指摘2のscan設定を先に解消し、推論方式の移行（指摘5）を同時に扱う。その後、3Diと非標準遺伝暗号の経路を接続し、各失敗例を回帰テストへ追加する。
