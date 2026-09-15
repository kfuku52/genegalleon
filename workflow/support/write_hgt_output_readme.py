#!/usr/bin/env python3

"""Write a column guide for the tabular HGT candidate outputs."""

import argparse
import csv
import os
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from scaffold_taxonomy import CONTEXT_COLUMNS
from score_hgt_candidates import (
    BRANCH_OUTPUT_COLUMNS,
    GENE_OUTPUT_COLUMNS,
    ORTHOGROUP_OUTPUT_COLUMNS,
    TAXONOMIC_RANKS,
    taxonomy_lineage_column,
    taxonomy_rank_column,
    taxonomy_rank_list_column,
)

ColumnSpec = Tuple[str, str]


COMMON_SPECS: Dict[str, ColumnSpec] = {
    "orthogroup": ("Orthogroup（遺伝子ファミリー）のID。", "文字列"),
}

BRANCH_SPECS: Dict[str, ColumnSpec] = {
    **COMMON_SPECS,
    "branch_id": (
        "orthogroupの`stat_branch`表における候補枝のID。",
        "整数",
    ),
    "node_name": (
        "遺伝子系統樹で候補枝に対応するノード名。",
        "文字列",
    ),
    "generax_event": (
        "GeneRaxイベントを正規化したラベル。候補抽出では`H`がtransferイベントを表し、周辺の系統樹では`S`=speciation、`D`=duplication、`L`=leafです。",
        "`H`, `S`, `D`, `L` など",
    ),
    "generax_transfer": (
        "`stat_branch`から保持した元のGeneRax transfer注釈。`Y`で始まる値はtransfer注釈を示し、`Y@src@dest`の形式なら有向transfer plotのsource/destinationとして利用できます。",
        "文字列（通常 `Y...` / 空欄）",
    ),
    "generax_event_parent": (
        "この候補枝の親枝におけるGeneRaxイベントの正規化ラベル。",
        "`H`, `S`, `D`, `L` など / 空欄",
    ),
    "taxon": (
        "`stat_branch`で候補枝に対応付けられたtaxonラベル。",
        "文字列 / 空欄",
    ),
    "spnode_coverage": (
        "`stat_branch`から引き継いだspecies-tree node / coverageラベル。",
        "文字列 / 空欄",
    ),
    "candidate_gene_count": (
        "候補GeneRax枝の下流にあるgene labelの数。",
        "0以上の整数",
    ),
    "matched_leaf_count": (
        "候補遺伝子のうち、`stat_branch`のleaf行に対応付けられた数。",
        "0以上の整数",
    ),
    "candidate_genes": (
        "この候補枝に割り当てられたgene IDのセミコロン区切り一覧。",
        "遺伝子IDの `; ` 区切り",
    ),
    "besthit_gene_count": (
        "best-hit accession、organism、またはtaxidのいずれかがある対応候補遺伝子の数。",
        "0以上の整数",
    ),
    "besthit_taxid_count": (
        "best-hitのtaxonomy IDを利用できたbest-hit比較の数。",
        "0以上の整数",
    ),
    "besthit_taxonomy_method": (
        "focal taxonとbest hitの比較に使われた方法の最頻値。taxonomy databaseまたはname heuristicです。",
        "`taxonomy_db`, `name_heuristic` など / 空欄",
    ),
    "besthit_same_superkingdom_fraction": (
        "比較可能なbest-hitのうち、hitがfocal lineageと同じsuperkingdomに属する割合。名前だけの比較では算出できません。",
        "0--1 / 空欄（比較不能）",
    ),
    "besthit_lca_rank_mode": (
        "best-hit比較で得られたLCA rank関係の最頻値。",
        "`species`, `genus`, `genus_mismatch` など / 空欄",
    ),
    "intron_measured_gene_count": (
        "イントロン状態を観測できた候補遺伝子の数。祖先状態のimputationは測定数に含めません。",
        "0以上の整数",
    ),
    "intron_supported_gene_count": (
        "測定された候補遺伝子のうち、観測イントロン数が1以上の遺伝子の数。",
        "0以上の整数",
    ),
    "intron_support_fraction": (
        "`intron_supported_gene_count / intron_measured_gene_count`。イントロン状態が未測定なら空欄。",
        "0--1 / 空欄（未測定）",
    ),
    "expression_measured_gene_count": (
        "少なくとも1つの数値expression測定値を持つ対応候補遺伝子の数。",
        "0以上の整数",
    ),
    "expression_measured_fraction": (
        "`expression_measured_gene_count / matched_leaf_count`。",
        "0--1 / 空欄（未測定）",
    ),
    "clade_min_expression_pearsoncor": (
        "expression-awareな枝統計から引き継いだ、cladeレベルPearson相関の最小値。",
        "通常 -1--1 / 空欄",
    ),
    "synteny_measured_gene_count": (
        "数値synteny support scoreを持つ候補遺伝子の数。",
        "0以上の整数",
    ),
    "synteny_supported_gene_count": (
        "測定された候補遺伝子のうち、synteny support scoreが正の遺伝子の数。",
        "0以上の整数",
    ),
    "synteny_support_fraction": (
        "`synteny_supported_gene_count / synteny_measured_gene_count`。",
        "0--1 / 空欄（未測定）",
    ),
    "synteny_mean_support_score": (
        "測定された候補遺伝子のsynteny support scoreの平均。scoreは、non-zero offsetのsynteny edgeのうち別の候補遺伝子と共有されるedgeの割合です。",
        "0--1 / 空欄（未測定）",
    ),
    "contamination_measured_gene_count": (
        "利用可能なcontamination lineage compatibility判定を持つ候補遺伝子の数。",
        "0以上の整数",
    ),
    "contamination_incompatible_gene_count": (
        "contamination QCでfocal lineageと非互換と判定された、測定済み候補遺伝子の数。",
        "0以上の整数",
    ),
    "contamination_incompatible_fraction": (
        "`contamination_incompatible_gene_count / contamination_measured_gene_count`。",
        "0--1 / 空欄（未測定）",
    ),
    "contamination_top_lca_taxid": (
        "contamination非互換遺伝子で最も多かったLCAのtaxonomy ID。",
        "Taxonomy ID / 空欄",
    ),
    "contamination_top_lca_sciname": (
        "contamination非互換遺伝子で最も多かったLCAの学名。",
        "学名 / 空欄",
    ),
    "contamination_top_lca_fraction": (
        "`contamination_top_lca_taxid` / `contamination_top_lca_sciname`に該当する、contamination非互換遺伝子の割合。",
        "0--1 / 空欄",
    ),
    "representative_gene_id": (
        "この候補枝の代表best-hit注釈を持つgene ID。best-hit accession / organism / taxidの組み合わせが最も多いgeneを選び、同数なら注釈の充足度とgene IDで決定します。",
        "gene ID / 空欄（候補geneなし）",
    ),
    "representative_gene_taxon": (
        "`representative_gene_id`のfocal taxonラベル。",
        "文字列 / 空欄",
    ),
    "representative_besthit_accession": (
        "`representative_gene_id`に対応する代表best-hit accession。",
        "accession / 空欄",
    ),
    "representative_besthit_organism": (
        "`representative_gene_id`に対応する代表best-hit organism名。",
        "文字列 / 空欄",
    ),
    "representative_besthit_taxid": (
        "`representative_gene_id`に対応する代表best-hit taxonomy ID。",
        "Taxonomy ID / 空欄",
    ),
    "representative_besthit_taxonomy_method": (
        "代表best-hitのtaxonomy比較に使った方法。",
        "`taxonomy_db`, `name_heuristic` など / 空欄",
    ),
    "representative_besthit_lca_rank": (
        "代表best-hitとfocal lineageのLCA rank関係。",
        "`species`, `genus`, `genus_mismatch` など / 空欄",
    ),
    "representative_besthit_same_superkingdom": (
        "代表best-hitがfocal lineageと同じsuperkingdomかどうか。",
        "0 / 1 / 空欄（比較不能）",
    ),
}

GENE_SPECS: Dict[str, ColumnSpec] = {
    **COMMON_SPECS,
    "gene_id": ("候補遺伝子のID（leaf / node name）。", "文字列"),
    "gene_taxon": ("候補遺伝子のtaxon / speciesラベル。", "文字列 / 空欄"),
    "candidate_branch_count": (
        "この遺伝子を含む、orthogroup内の異なる候補枝の数。",
        "0以上の整数",
    ),
    "candidate_branch_ids": (
        "この遺伝子を含む候補枝のID一覧（セミコロン区切り）。",
        "整数IDの `; ` 区切り",
    ),
    "besthit_accession": ("annotation tableから保持したbest-hit accession。", "文字列 / 空欄"),
    "besthit_organism": ("best-hitのorganism名。", "文字列 / 空欄"),
    "besthit_taxid": ("利用可能な場合のbest-hit taxonomy ID。", "Taxonomy ID / 空欄"),
    "besthit_taxonomy_method": (
        "best-hit関係の分類に使った方法。taxonomy databaseまたはname heuristicです。",
        "`taxonomy_db`, `name_heuristic` など / 空欄",
    ),
    "besthit_lca_rank": (
        "この遺伝子のbest hitに対するLCA rank関係。",
        "`species`, `genus`, `genus_mismatch` など / 空欄",
    ),
    "besthit_same_superkingdom": (
        "best hitがfocal lineageと同じsuperkingdomかどうか。1=同じ、0=異なる。名前だけの比較ではsuperkingdomを確定できないため空欄です。",
        "0 / 1 / 空欄（比較不能）",
    ),
    "intron_supported": (
        "観測イントロン状態がイントロンありを支持するか。観測イントロン数が0より大きければ`True`、それ以外は`False`。imputed状態は使いません。イントロン入力がない場合も`False`になるため、枝単位の測定数と併読してください。",
        "`True` / `False`",
    ),
    "expression_measured": (
        "この遺伝子について少なくとも1つの数値expression値が測定されたかどうか。expression入力自体がない場合も`False`です。",
        "`True` / `False`",
    ),
    "synteny_support_score": (
        "この遺伝子のsynteny support。non-zero offsetのsynteny edge数に対する、共有synteny edge数の割合。",
        "0--1 / 空欄（未測定）",
    ),
    "contamination_lca_taxid": (
        "contamination QCがこの遺伝子について報告したLCA taxonomy ID。",
        "Taxonomy ID / 空欄",
    ),
    "contamination_lca_sciname": (
        "contamination QCがこの遺伝子について報告したLCA学名。",
        "学名 / 空欄",
    ),
    "contamination_is_compatible_lineage": (
        "contamination QCにおけるfocal lineageとの互換性。`True`=互換、`False`=非互換。空欄は利用可能なQC結果がないことを示します。",
        "`True` / `False` / 空欄（未測定）",
    ),
}

ORTHOGROUP_SPECS: Dict[str, ColumnSpec] = {
    **COMMON_SPECS,
    "hgt_branch_count": (
        "orthogroup内の異なるGeneRax transfer候補枝の数。",
        "0以上の整数",
    ),
    "hgt_gene_count": (
        "それらの枝から集約された異なる候補遺伝子の数。",
        "0以上の整数",
    ),
    "candidate_branch_ids": (
        "orthogroup内の全transfer候補枝のID一覧（セミコロン区切り）。",
        "整数IDの `; ` 区切り",
    ),
}


for _rank in TAXONOMIC_RANKS:
    _rank_label = _rank.replace("_", " ")
    _recipient_gene_column = taxonomy_rank_column("recipient", _rank)
    _donor_gene_column = taxonomy_rank_column("donor", _rank)
    _recipient_list_column = taxonomy_rank_list_column("recipient", _rank)
    _donor_list_column = taxonomy_rank_list_column("donor", _rank)
    BRANCH_SPECS[_recipient_list_column] = (
        f"この候補枝に含まれるcandidate geneのfocal taxonをNCBI taxonomyの{_rank_label} rankへ解決した、重複なしの一覧。",
        f"{_rank_label}名の `; ` 区切り / 空欄（未解決）",
    )
    BRANCH_SPECS[_donor_list_column] = (
        f"この候補枝に含まれるbest-hit organism / taxidをNCBI taxonomyの{_rank_label} rankへ解決した、重複なしの一覧。推定donor側の分類であり、HGT donorの確定判定ではありません。",
        f"{_rank_label}名の `; ` 区切り / 空欄（未解決）",
    )
    GENE_SPECS[_recipient_gene_column] = (
        f"候補遺伝子の`gene_taxon`をNCBI taxonomyの{_rank_label} rankへ解決した名前。recipientはfocal側です。",
        f"{_rank_label}名 / 空欄（未解決）",
    )
    GENE_SPECS[_donor_gene_column] = (
        f"候補遺伝子のbest-hit organism / taxidをNCBI taxonomyの{_rank_label} rankへ解決した名前。donorはbest-hit側のproxyであり、確定donorを意味しません。",
        f"{_rank_label}名 / 空欄（未解決）",
    )
    ORTHOGROUP_SPECS[_recipient_list_column] = (
        f"orthogroup内の候補遺伝子のfocal taxonをNCBI taxonomyの{_rank_label} rankへ解決した、重複なしの一覧。",
        f"{_rank_label}名の `; ` 区切り / 空欄（未解決）",
    )
    ORTHOGROUP_SPECS[_donor_list_column] = (
        f"orthogroup内のbest-hit organism / taxidをNCBI taxonomyの{_rank_label} rankへ解決した、重複なしの一覧。推定donor側の分類であり、HGT donorの確定判定ではありません。",
        f"{_rank_label}名の `; ` 区切り / 空欄（未解決）",
    )
    BRANCH_SPECS[f"representative_{_recipient_gene_column}"] = (
        f"`representative_gene_id`のfocal taxonを{_rank_label} rankへ解決した名前。",
        f"{_rank_label}名 / 空欄（未解決）",
    )
    BRANCH_SPECS[f"representative_{_donor_gene_column}"] = (
        f"`representative_gene_id`のbest-hit lineageを{_rank_label} rankへ解決した名前。推定donor側の分類です。",
        f"{_rank_label}名 / 空欄（未解決）",
    )

_host_common = {
    "host_scaffold_status": ("scaffold分類組成の取得状態。measuredは取得済みであり高確度HGTを意味しません。partial=recipientの一部のみ対応、missing_scaffold_taxonomy=入力なし、gene_not_mapped=座標対応なし、recipient_unresolved=transfer/tree対応不明、no_recipient_scaffold=recipient側座標なし。", "状態ラベル"),
}
GENE_SPECS.update(_host_common)
BRANCH_SPECS.update(_host_common)
GENE_SPECS.update({
    "host_scaffold_id": ("遺伝子が載るscaffold。gene_taxonとの組で識別します。", "文字列 / 空欄"),
    "host_scaffold_locus_id": ("GFF上のlocus ID。対応不明ならCDS IDを使います。", "文字列 / 空欄"),
    "host_scaffold_count_unit": ("gff_locus=明示的GFF親子関係でisoformを統合、cds_id=関係不明でCDS ID単位。後者はisoform過大計数の可能性があります。", "gff_locus / cds_id / 空欄"),
})
BRANCH_SPECS.update({
    "host_scaffold_recipient_gene_count": ("候補遺伝子のうち、GeneRax recipient枝配下の現生種に属する数。", "整数 / 空欄"),
    "host_scaffold_mapped_gene_count": ("recipient側候補のうちscaffold分類表に対応した遺伝子数。", "整数 / 空欄"),
    "host_scaffold_count": ("対応したspecies+scaffoldの重複なしの数。複数候補が同じscaffoldにあっても1件。", "整数 / 空欄"),
})
_metric_descriptions = {
    "total_count": "集計対象のlocus数（locus不明はCDS ID単位）。",
    "compatible_count": "当該rankでhost taxidと一致するlocus数。",
    "incompatible_count": "当該rankでhost taxidと異なるlocus数。",
    "unresolved_count": "遺伝子またはhostが当該rankまで分類できないlocus数。isoform間分類不一致も含む。",
    "cds_id_count": "total_countのうち明示的GFF locusに対応できずCDS ID単位で数えた数。isoform過大計数の可能性を点検するための値。",
    "classified_fraction": "(compatible_count + incompatible_count) / total_count。分類可能率。",
    "compatible_fraction": "compatible_count / (compatible_count + incompatible_count)。分類可能なものに対する宿主分類群一致率。",
    "compatible_all_fraction": "compatible_count / total_count。未分類も分母に含む宿主分類群一致率。",
}
for _column in CONTEXT_COLUMNS:
    _suffix = _column.removeprefix("host_scaffold_")
    _background = _suffix.startswith("background_")
    _rank, _metric = _suffix.removeprefix("background_").split("_", 1)
    _scope = "全orthogroupのHGT候補locusの和集合を除外した背景。" if _background else "候補も含むscaffold全体。"
    _description = f"rank={_rank}。{_scope}{_metric_descriptions[_metric]}分母0なら空欄。HGT確率ではありません。"
    _range = "0--1 / 空欄" if "fraction" in _metric else "0以上の整数 / 空欄"
    GENE_SPECS[_column] = (_description, _range)
    BRANCH_SPECS[_column] = ("recipient側の重複なしspecies+scaffoldについてlocus数を合算。" + _description, _range)

BRANCH_SPECS[taxonomy_lineage_column("recipient", plural=True)] = (
    "候補枝に含まれるfocal taxonの、taxonomy DBが返した全named rankと名前を重複なしで保持します。各値は`rank:name`形式です。",
    "`rank:name; ...` / 空欄（未解決）",
)
BRANCH_SPECS[taxonomy_lineage_column("donor", plural=True)] = (
    "候補枝に含まれるbest-hit lineageの、taxonomy DBが返した全named rankと名前を重複なしで保持します。推定donor側のproxyです。",
    "`rank:name; ...` / 空欄（未解決）",
)
GENE_SPECS[taxonomy_lineage_column("recipient")] = (
    "候補遺伝子のfocal taxonについて、taxonomy DBが返した全named rankと名前を保持します。各値は`rank:name`形式です。",
    "`rank:name; ...` / 空欄（未解決）",
)
GENE_SPECS[taxonomy_lineage_column("donor")] = (
    "候補遺伝子のbest-hit lineageについて、taxonomy DBが返した全named rankと名前を保持します。推定donor側のproxyです。",
    "`rank:name; ...` / 空欄（未解決）",
)
BRANCH_SPECS["representative_recipient_taxonomy"] = (
    "代表遺伝子のfocal taxonについて、taxonomy DBが返した全named rankと名前を保持します。",
    "`rank:name; ...` / 空欄（未解決）",
)
BRANCH_SPECS["representative_donor_taxonomy"] = (
    "代表遺伝子のbest-hit lineageについて、taxonomy DBが返した全named rankと名前を保持します。推定donor側のproxyです。",
    "`rank:name; ...` / 空欄（未解決）",
)


for _side in ("recipient", "donor"):
    _lineage_column = taxonomy_lineage_column(_side, plural=True)
    ORTHOGROUP_SPECS[_lineage_column] = BRANCH_SPECS[_lineage_column]


TABLES: Sequence[Tuple[str, str, Sequence[str], Dict[str, ColumnSpec]]] = (
    (
        "hgt_branch_candidates.tsv",
        "候補枝単位（1行 = 1 orthogroup × 1 transfer-candidate branch）",
        BRANCH_OUTPUT_COLUMNS,
        BRANCH_SPECS,
    ),
    (
        "hgt_gene_candidates.tsv",
        "遺伝子単位（1行 = 1 orthogroup × 1 gene_id。複数枝に現れる場合は集約）",
        GENE_OUTPUT_COLUMNS,
        GENE_SPECS,
    ),
    (
        "hgt_orthogroup_summary.tsv",
        "orthogroup単位（1行 = 1 orthogroup）",
        ORTHOGROUP_OUTPUT_COLUMNS,
        ORTHOGROUP_SPECS,
    ),
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Write the README for GeneGalleon HGT output tables.")
    parser.add_argument("--output", required=True, type=str, help="README.md path to write")
    parser.add_argument("--branch_tsv", required=True, type=str)
    parser.add_argument("--gene_tsv", required=True, type=str)
    parser.add_argument("--orthogroup_tsv", required=True, type=str)
    return parser


def read_header(path: str) -> List[str]:
    if not os.path.isfile(path):
        return []
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            row = next(csv.reader(handle, delimiter="\t"), [])
    except (OSError, UnicodeError, csv.Error):
        return []
    return [str(value).strip() for value in row if str(value).strip()]


def format_column_rows(
    columns: Sequence[str], expected_columns: Sequence[str], specs: Dict[str, ColumnSpec]
) -> List[str]:
    ordered = list(columns) if columns else list(expected_columns)
    for column in expected_columns:
        if column not in ordered:
            ordered.append(column)
    rows = []
    expected_set = set(expected_columns)
    for column in ordered:
        if column in specs:
            meaning, value_range = specs[column]
        else:
            meaning = "このバージョンで定義されていない追加列。生成元の変更内容を確認してください。"
            value_range = "未定義"
        if columns and column not in expected_set:
            meaning = "追加列。" + meaning
        if columns and column not in columns:
            meaning = "期待される列ですが、このファイルのヘッダーには存在しません。" + meaning
        rows.append(f"| `{column}` | {meaning} | {value_range} |")
    return rows


def build_readme(paths: Dict[str, str]) -> str:
    lines = [
        "# GeneGalleon HGT output tables",
        "",
        "このREADMEは、HGT評価で生成される3つのTSVの列定義です。",
        "",
        "## まず押さえる点",
        "",
        "- HGT候補の一次条件はGeneRaxのtransferイベント（`generax_event == H`、またはtransfer注釈が`Y`で始まる枝）です。",
        "- これらの表は候補と補助的な証拠を整理したもので、単一のHGT確率や校正済み総合スコアではありません。",
        "- `*_fraction` と `*_score` は基本的に0--1ですが、分母となる測定値がない場合は空欄です。空欄は「陰性」ではなく「未測定・比較不能」を意味します。",
        "- `hgt_branch_candidates.tsv` → `hgt_gene_candidates.tsv` → `hgt_orthogroup_summary.tsv` の順に、枝・遺伝子・orthogroupへ集約されています。",
        "- `candidate_branch_count` / `hgt_branch_count` は枝数、`candidate_gene_count` / `hgt_gene_count` は遺伝子数です。これらはカウントであり、信頼度スコアではありません。",
        "- taxonomy rankで表を絞り込む場合は、gene表の`recipient_<rank>` / `donor_<rank>`を使ってください。branch・orthogroup表には同じrankの重複なしリスト列があります。",
        "- rank列はtaxonomy DBで該当rankを解決できた場合だけ埋まります。空欄は「そのrankではない」ではなく、taxonomy DB・taxid・lineageのいずれかが不足していることを示します。`recipient_taxonomy` / `donor_taxonomy` は標準列にないnamed rankも含む全lineageです。",
        "- branch表の`representative_*`列は枝内の全gene注釈を置き換えるものではありません。best-hit注釈の最頻組み合わせから1件を抜き出した代表値なので、全遺伝子の詳細はgene表で確認してください。",
        "- plot出力の`hgt_transfer_edges.tsv`は`generax_transfer`の`Y@src@dest`を方向別イベント数へ集約した表です。`hgt_transfer_tree.pdf`は双方向を1本の曲線で示し、各先端側半分の太さがその先端へ向かうイベント数を表します（最小0.35 pt）。遠距離ほど濃い青です。表示は件数順位と距離順位を交互に採用し、既存の逆方向を追加します。`phylogenetic_distance`は端点ノード間の経路長、`distance_metric`はbranch_lengthまたはtopology_edges、`selection_reason`はcount/distance/all/reciprocal/not_displayedです。未対応ペアの距離は欠損です。詳しくはplots/README.mdを参照してください。",
        "",
        "## 表同士の対応",
        "",
        "| 表 | 行の単位 | 主な対応関係 |",
        "|---|---|---|",
        "| `hgt_branch_candidates.tsv` | orthogroup × 候補枝 | `candidate_genes` に枝の下流遺伝子を列挙 |",
        "| `hgt_gene_candidates.tsv` | orthogroup × gene_id | `candidate_branch_ids` にその遺伝子を含む候補枝を列挙 |",
        "| `hgt_orthogroup_summary.tsv` | orthogroup | 候補枝数・候補遺伝子数を集約 |",
        "",
    ]
    for filename, grain, expected_columns, specs in TABLES:
        path = paths[filename]
        actual_columns = read_header(path)
        status = "ヘッダー確認済み" if actual_columns else "ファイル未生成またはヘッダーを読めません"
        lines.extend(
            [
                f"## `{filename}`",
                "",
                f"粒度: {grain}。{status}。",
                "",
                "| 列 | 意味 | 値域・欠測時の扱い |",
                "|---|---|---|",
            ]
        )
        lines.extend(format_column_rows(actual_columns, expected_columns, specs))
        lines.append("")
    lines.extend(
        [
            "## 解釈上の注意",
            "",
            "- best-hitの分類が名前ヒューリスティックだけで行われた場合、superkingdomの比較はできないため、`besthit_same_superkingdom` 系の値が空欄になることがあります。",
            "- イントロンは観測値を優先し、祖先状態のimputationだけではsupportに数えません。入力`stat_branch`に認識可能なイントロン列がない場合は未測定扱いです。",
            "- gene表の`intron_supported=False`は、観測イントロンがない場合だけでなく入力が未測定の場合にも現れます。枝表の`intron_measured_gene_count`で測定の有無を確認してください。",
            "- expressionの`True`は少なくとも1つの数値測定があることだけを示し、発現量の大きさを表しません。`False`は入力がない場合も含みます。",
            "- `hgt_` evidence heatmapのPDF表示では、各列をleaf内の最大値で割ってcolorbarを0--1に正規化しています。TSVのカウント値自体は変更されません。",
            "- 欠測値を0として扱って候補を否定しないでください。まず対応する`*_measured_*`列または入力データの有無を確認してください。",
            "",
            "このREADMEはGeneGalleonのHGT出力ステージが、現在のTSVヘッダーと出力スキーマから自動生成します。",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    args = build_arg_parser().parse_args()
    paths = {
        "hgt_branch_candidates.tsv": args.branch_tsv,
        "hgt_gene_candidates.tsv": args.gene_tsv,
        "hgt_orthogroup_summary.tsv": args.orthogroup_tsv,
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(build_readme(paths), encoding="utf-8")


if __name__ == "__main__":
    main()
