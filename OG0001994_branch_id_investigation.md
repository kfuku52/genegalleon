# OG0001994 site13: CSUBST / `stat_branch` branch ID 不一致の調査手順

## この文書の対象

`OG0001994` の site 13 について、CSUBSTが出力したbranch IDと
GeneGalleonの`stat_branch`のbranch IDが異なる枝を指している原因を調べる。

具体的には、CSUBST側でforegroundとして選ばれたbranch ID 94が、
`.report.pdf`では`Hexapleomera_sasuke_DN1845-c1-g1`のterminal branchに
誤って対応付けられている問題を対象とする。

この調査が完了するまでは、stat file不足で停止する別のHOGの問題には着手しない。

## AI実行担当者への重要な指示

1. **調査前に元の解析を再実行しないこと。**
2. **元の出力を削除、移動、上書き、`touch`しないこと。**
3. 最初に元のGit revision、ファイルのmtime、サイズ、SHA-256を記録すること。
4. 調査用コピーはGeneGalleon repository内の`tmp/`以下に作ること。
5. 出力ディレクトリを一意に特定できない場合は、推測して先へ進まず候補を報告すること。
6. 不一致をbranch IDの読み替えで補正しないこと。不一致そのものを原因究明の証拠として扱うこと。
7. 調査終了時に「判明した事実」と「まだ推測に留まること」を分けて報告すること。

## 想定する原因

GeneGalleonの旧処理では、既存の`iqtree_anc` archiveが存在すると、
rooted treeまたはtrimmed alignmentが後から更新されてもarchiveを再生成しない場合があった。

その場合、次のようなmixed provenanceが発生し得る。

1. `iqtree_anc`、`csubst_b`、`csubst_scan`は古いrooted tree由来
2. `stat_branch`は新しいrooted tree由来
3. 両者のbranch ID採番自体は正常だが、入力木が違うため同じ番号が別クレードを指す
4. CSUBSTのbranch IDをそのまま`stat_branch`へ渡すと、別の枝が赤く描画される

ただし、元の実行で実際にどのファイルが更新されたかは、以下の調査結果から判定する。

## 1. GeneGalleon revisionの保存

この操作はコードをpullする**前**に行う。

```bash
cd /path/to/genegalleon

mkdir -p "$PWD/tmp/OG0001994_branch-audit"
git rev-parse HEAD \
  > "$PWD/tmp/OG0001994_branch-audit/genegalleon_revision_before_pull.txt"
git status --short \
  > "$PWD/tmp/OG0001994_branch-audit/git_status_before_pull.txt"
```

ローカル変更がある場合は、それも報告対象に含める。勝手にdiscardしない。

## 2. 元出力ディレクトリの特定

`OUT`には、直下に`rooted_tree/`、`iqtree_anc/`、`csubst_b/`、
`csubst_scan/`、`stat_branch/`などが存在するディレクトリを指定する。

場所が分からない場合は、解析workspaceとして妥当な範囲に限定して検索する。

```bash
find /path/to/workspace \
  -type f \
  -name 'OG0001994_iqtree.anc.zip' \
  -print
```

候補が複数ある場合は、各候補について次のファイルの有無を確認し、
勝手に一つを選ばず候補一覧を報告する。

```text
rooted_tree/OG0001994_root.nwk
clipkit/OG0001994_cds.clipkit.fa.gz
iqtree_anc/OG0001994_iqtree.anc.zip
csubst_b/OG0001994_csubst_b.tsv
csubst_cb_stats/OG0001994_csubst_cb_stats.tsv
csubst_scan/OG0001994_csubst_scan.tsv
csubst_scan_units/OG0001994_csubst_scan_units.tsv
stat_branch/OG0001994_stat.branch.tsv
```

一意に特定できた場合のみ、以下を設定する。

```bash
cd /path/to/genegalleon

OG=OG0001994
OUT=/absolute/path/to/directory-containing-rooted_tree-and-stat_branch
AUDIT="$PWD/tmp/${OG}_branch-audit"
mkdir -p "$AUDIT"
```

## 3. 元成果物の証拠保全

以下は元ファイルを変更せず、`cp -a`で調査用コピーを作成する。

```bash
RELATIVE_FILES=(
  "rooted_tree/${OG}_root.nwk"
  "clipkit/${OG}_cds.clipkit.fa.gz"
  "iqtree_anc/${OG}_iqtree.anc.zip"
  "csubst_b/${OG}_csubst_b.tsv"
  "csubst_cb_stats/${OG}_csubst_cb_stats.tsv"
  "csubst_scan/${OG}_csubst_scan.tsv"
  "csubst_scan_units/${OG}_csubst_scan_units.tsv"
  "stat_branch/${OG}_stat.branch.tsv"
)

: > "$AUDIT/timestamps.txt"
: > "$AUDIT/sha256.txt"
: > "$AUDIT/missing_files.txt"

for rel in "${RELATIVE_FILES[@]}"; do
  src="$OUT/$rel"
  if [[ -e "$src" ]]; then
    cp -a "$src" "$AUDIT/"
    stat -c '%y\t%s\t%n' "$src" >> "$AUDIT/timestamps.txt"
    sha256sum "$src" >> "$AUDIT/sha256.txt"
  else
    printf 'MISSING\t%s\n' "$src" >> "$AUDIT/missing_files.txt"
  fi
done

unzip -l "$AUDIT/${OG}_iqtree.anc.zip" \
  > "$AUDIT/iqtree_archive_contents.txt"
```

`missing_files.txt`が空でない場合も調査は中止せず、欠損を結果として報告する。
ただし、`iqtree_anc`または`stat_branch`が欠損している場合はbranch ID照合を実行できないため、
その時点で停止して不足ファイルを報告する。

## 4. 実行設定とログの保存

次の情報を探し、見つかったものを`AUDIT`へコピーする。

- 元解析に使用したGeneGalleon設定ファイル
- 元解析の標準出力・標準エラーログ
- `run_iqtree_anc`、`run_csubst`、`run_csubst_scan`、`run_summary`の値
- 実行日時と実行コマンド
- 使用したcontainer imageまたはSIFの識別情報

設定ファイルを特定できる場合は、少なくとも次を保存する。

```bash
grep -E '^(run_iqtree_anc|run_csubst|run_csubst_scan|run_summary)=' \
  /path/to/original/config \
  > "$AUDIT/relevant_run_flags.txt" || true
```

ファイルが見つからなかった場合は、推測した値を書かず「未発見」と報告する。

## 5. `iqtree_anc` archiveの展開

```bash
EXTRACT="$AUDIT/iqtree-extracted"
if [[ -e "$EXTRACT" ]]; then
  echo "Extraction destination already exists; refusing to overwrite: $EXTRACT" >&2
  exit 1
fi
mkdir -p "$EXTRACT"

unzip -q \
  "$AUDIT/${OG}_iqtree.anc.zip" \
  -d "$EXTRACT"

find "$EXTRACT" \
  -maxdepth 2 \
  -type f \
  -print \
  | sort \
  > "$AUDIT/iqtree_extracted_files.txt"
```

次の2ファイルが存在することを確認する。

```text
tmp/OG0001994_branch-audit/iqtree-extracted/OG0001994.iqtree.anc/csubst.nwk
tmp/OG0001994_branch-audit/iqtree-extracted/OG0001994.iqtree.anc/csubst.treefile
```

## 6. CSUBSTと`stat_branch`の全枝照合

この照合は、表示対象のbranch ID 94だけではなく、両方の木に含まれる**全branch ID**について行う。
各branch IDを、その枝以下に存在する全tip名の集合へ変換して比較する。

修正版GeneGalleonを取得した後、repository rootから次を実行する。

```bash
bash workflow/tests/run_in_sif.sh \
  python - \
  "$AUDIT/${OG}_stat.branch.tsv" \
  "$AUDIT/iqtree-extracted/${OG}.iqtree.anc" \
  > "$AUDIT/branch_identity.stdout.txt" \
  2> "$AUDIT/branch_identity.stderr.txt" <<'PY'
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd() / "workflow" / "support"))
import csubst_site_wrapper as wrapper

wrapper.validate_csubst_stat_branch_identity(
    branch_id_str="8,26,87,94,140,154,198,268,292,303",
    file_stat_branch=sys.argv[1],
    iqtree_anc_dir=sys.argv[2],
)
PY

status=$?
printf '%s\n' "$status" > "$AUDIT/branch_identity.exit_status.txt"
cat "$AUDIT/branch_identity.stdout.txt"
cat "$AUDIT/branch_identity.stderr.txt" >&2
```

シェルに`set -e`が設定されている場合は、上のコマンドを実行する直前に`set +e`を実行し、
不一致によるnon-zero statusでも調査用ファイルを保存できるようにする。

### 正常な場合

exit statusが0で、次のように表示される。

```text
Validated identical CSUBST/stat_branch branch IDs for ... rooted clades.
```

### 不一致の場合

non-zeroのexit statusになり、例えば次の形式で最初の不一致が表示される。

```text
CSUBST and stat_branch branch IDs were generated from inconsistent rooted trees.
Refusing to remap IDs because that would hide stale or mixed-provenance artifacts.
...
94: CSUBST=[...] stat_branch=[...]
```

エラーは期待される調査結果になり得るため、ファイルを削除して再試行しない。

branch ID 94の対応を省略せず保存するため、続けて次を実行する。

```bash
bash workflow/tests/run_in_sif.sh \
  python - \
  "$AUDIT/${OG}_stat.branch.tsv" \
  "$AUDIT/iqtree-extracted/${OG}.iqtree.anc" \
  > "$AUDIT/branch_94_descendant_tips.txt" <<'PY'
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd() / "workflow" / "support"))
import csubst_site_wrapper as wrapper

csubst_signatures = wrapper.get_csubst_branch_clade_signatures(sys.argv[2])
stat_signatures = wrapper.get_stat_branch_clade_signatures(sys.argv[1])

for source, signatures in (
    ("CSUBST", csubst_signatures),
    ("stat_branch", stat_signatures),
):
    descendant_tips = signatures.get(94)
    value = "MISSING" if descendant_tips is None else ";".join(sorted(descendant_tips))
    print(f"{source}\tbranch_id=94\t{value}")
PY

cat "$AUDIT/branch_94_descendant_tips.txt"
```

## 7. mtimeによる更新順序の判定

`timestamps.txt`を確認し、少なくとも次のファイルの順序を整理する。

1. `rooted_tree/OG0001994_root.nwk`
2. `clipkit/OG0001994_cds.clipkit.fa.gz`
3. `iqtree_anc/OG0001994_iqtree.anc.zip`
4. `csubst_b/OG0001994_csubst_b.tsv`
5. `csubst_scan/OG0001994_csubst_scan.tsv`
6. `stat_branch/OG0001994_stat.branch.tsv`

次の条件が揃えば、stale `iqtree_anc` cacheが原因だったと強く判定できる。

- rooted treeまたはtrimmed alignmentが`iqtree_anc` archiveより新しい
- `stat_branch`がその後に生成または更新されている
- 全枝照合が不一致になる
- CSUBST側と`stat_branch`側でtip名そのものは概ね共通している

mtimeが上記を示さなくても、全枝照合が不一致ならmixed provenance自体は確認できる。
コピー、展開、archive作成時にmtimeが保存・変更されることがあるため、
mtimeだけを根拠に「mixed provenanceではない」と結論しない。

## 8. 既存PDF／report cacheの確認

OG0001994に関する既存PDFを検索して一覧を保存する。

```bash
find "$(dirname "$OUT")" \
  -type f \
  -path '*OG0001994*' \
  \( -name '*OG0001994*report.pdf' \
     -o -name '*OG0001994*focused_tree_site.pdf' \
     -o -name 'csubst.tree_site.pdf' \) \
  -print \
  > "$AUDIT/existing_OG0001994_pdfs.txt"
```

全枝照合が正常なのに既存PDFだけが誤っている場合は、
PDFまたはcandidate report directoryが古い入力から作られたまま再利用された可能性を調べる。
この場合も、既存PDFを先に削除せずmtime、SHA-256、親ディレクトリ名を保存する。

## 9. 判定基準

### 判定A: stale `iqtree_anc` / mixed provenance

- 全枝照合が不一致
- rooted tree/alignmentと`iqtree_anc`の生成時期が異なる
- CSUBSTと`stat_branch`で同じIDが異なるdescendant-tip集合を指す

この場合、branch ID採番ロジックそのものの故障ではなく、
異なるrooted treeから作られた成果物を組み合わせたことが原因である。

### 判定B: report/PDF cacheのみが古い

- 全枝照合は正常
- 現在のCSUBSTと`stat_branch`ではbranch ID 94が同じ枝を指す
- 既存PDFだけが異なる枝を赤く表示している

この場合、candidate reportまたはPDFの再利用経路を追加調査する。

### 判定C: 現時点では特定不能

- 必要ファイルが欠損
- 元revisionまたは元設定が不明
- 複数の出力ディレクトリがあり、実際にreportへ使われたものを特定できない

この場合は推測でAまたはBに分類せず、追加で必要なファイルを列挙する。

## 10. 調査後の修復

修復は、上記の証拠保全と原因判定が終わってから行う。

修正版では次のように処理される。

1. rooted treeまたはtrimmed alignmentが新しければ`iqtree_anc`を再生成
2. `iqtree_anc`またはtrait fileが新しければ`csubst_b`と`csubst_scan`を再生成
3. 更新されたCSUBST出力から`stat_branch`を再生成
4. report生成前にCSUBSTと`stat_branch`の全branch IDを照合
5. 不一致ならbranch IDを読み替えず、report生成を停止

既存の誤ったcandidate report/PDFは、証拠保全後に別名へ退避してから再生成する。
元成果物を直接上書きしてはならない。

## 11. AI実行担当者が返す報告形式

最終報告には、最低限以下を含める。

```text
1. GeneGalleon revision:
2. 調査したOUTの絶対パス:
3. 見つからなかった必要ファイル:
4. rooted tree / alignment / iqtree_anc / csubst / stat_branchのmtime順序:
5. 全枝照合のexit status:
6. 全枝照合の最初の不一致（最大5件）:
7. branch ID 94がCSUBST側で指すdescendant tips:
8. branch ID 94がstat_branch側で指すdescendant tips:
9. 判定（A / B / C）:
10. 判明した事実:
11. まだ推測に留まる点:
12. 追加で必要なデータ:
```

`AUDIT`ディレクトリは調査完了後も削除せず、必要に応じてarchive化して共有する。
