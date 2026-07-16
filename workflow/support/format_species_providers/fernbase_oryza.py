"""Provider resolver implementation: fernbase oryza."""

import os
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote, urlparse

from format_species_annotations import (
    choose_first_gff_attribute,
    collapse_transcript_suffix,
    extract_header_tag_value,
    extract_provider_id,
    iter_fasta_records,
    merge_coordinate_intervals,
    parse_gff_attributes,
    resolve_feature_gene_token,
    strip_suffix_case_insensitive,
)
from format_species_common import (
    FERNBASE_GFF_EXCLUDE_PATTERN,
    provider_candidate_sort_key,
)
from format_species_constants import (
    FASTA_EXTENSIONS,
    FERNBASE_COMBINED_FILENAME_MARKER,
    FERNBASE_CONFIDENCE_MODE_FIELD,
    FERNBASE_CONFIDENCE_MODE_HIGH_LOW_COMBINED,
    FERNBASE_CONFIDENCE_MODE_HIGH_ONLY,
    GENBANK_EXTENSIONS,
    GFF_EXTENSIONS,
    ORYZA_MINUTA_CANONICAL_SOURCE_ID,
    ORYZA_MINUTA_DEFAULT_SPECIES_KEY,
    ORYZA_MINUTA_SOURCE_ID_ALIASES,
    ORYZA_MINUTA_SUBGENOME_BUNDLES,
)
from format_species_provider_config import (
    ORYZA_MINUTA_PROVIDER,
)
from format_species_provider_urls import (
    fetch_text_with_headers,
    is_url_like,
    parse_links_from_document,
    render_id_url_template,
    resolve_oryza_minuta_gramene_base_url,
    strip_provider_prefix,
)
from format_species_taxonomy import (
    normalize_species_key_for_runtime,
)
from format_species_writers import open_text, write_fasta_records, write_text_lines

from .common import (
    is_probable_cds_url,
    is_probable_gff_url,
    resolve_urls_from_index_url,
    select_best_url_for_label,
    source_id_candidates,
)


def parse_fernbase_confidence_mode(raw_value):
    text = str(raw_value or "").strip()
    if text == "":
        return ""
    normalized = re.sub(r"[_-]+", " ", text.lower())
    normalized = re.sub(r"\s+", " ", normalized).strip()
    if normalized in ("high confidence only", "high only", "high"):
        return FERNBASE_CONFIDENCE_MODE_HIGH_ONLY
    if normalized in ("high low combined", "combined", "merge", "merged"):
        return FERNBASE_CONFIDENCE_MODE_HIGH_LOW_COMBINED
    raise ValueError(
        "invalid {} '{}'; expected one of: {}, {}".format(
            FERNBASE_CONFIDENCE_MODE_FIELD,
            raw_value,
            FERNBASE_CONFIDENCE_MODE_HIGH_ONLY,
            FERNBASE_CONFIDENCE_MODE_HIGH_LOW_COMBINED,
        )
    )


def effective_fernbase_confidence_mode(raw_value):
    parsed = parse_fernbase_confidence_mode(raw_value)
    if parsed != "":
        return parsed
    return FERNBASE_CONFIDENCE_MODE_HIGH_ONLY


def build_fernbase_combined_filename(name):
    text = str(name or "").strip()
    if text == "":
        return text
    if FERNBASE_COMBINED_FILENAME_MARKER in text.lower():
        return text
    replaced = re.sub(r"highconfidence", FERNBASE_COMBINED_FILENAME_MARKER, text, count=1, flags=re.IGNORECASE)
    if replaced != text:
        return replaced
    replaced = re.sub(r"lowconfidence", FERNBASE_COMBINED_FILENAME_MARKER, text, count=1, flags=re.IGNORECASE)
    if replaced != text:
        return replaced
    suffixes = FASTA_EXTENSIONS + GFF_EXTENSIONS + GENBANK_EXTENSIONS
    stem = strip_suffix_case_insensitive(text, suffixes)
    if stem != text:
        return stem + "." + FERNBASE_COMBINED_FILENAME_MARKER + text[len(stem) :]
    return text + "." + FERNBASE_COMBINED_FILENAME_MARKER


def fernbase_feature_gene_token_from_header(header):
    gene_symbol = extract_header_tag_value(header, "gene")
    extracted = extract_provider_id("fernbase", header)
    collapsed = collapse_transcript_suffix("fernbase", extracted)
    if gene_symbol != "":
        if (
            collapsed != ""
            and collapsed != gene_symbol
            and re.search(r"(?:^|[._-]){}$".format(re.escape(gene_symbol)), collapsed, re.IGNORECASE)
        ):
            return collapsed
        return gene_symbol
    if collapsed != "":
        return collapsed
    return extracted


def list_directory_file_urls(reference, timeout, headers):
    text = str(reference or "").strip()
    if text == "":
        return []
    parsed = urlparse(text)
    scheme = parsed.scheme.lower()
    if scheme == "file":
        path = Path(unquote(parsed.path)).expanduser().resolve()
        directory = path if path.is_dir() else path.parent
        return [child.resolve().as_uri() for child in sorted(directory.iterdir()) if child.is_file()]

    index_url = text
    if not index_url.endswith("/"):
        index_url = index_url.rsplit("/", 1)[0] + "/"
    index_text = fetch_text_with_headers(index_url, timeout, headers)
    links = parse_links_from_document(index_url, index_text)
    out = []
    seen = set()
    for link in links:
        path = urlparse(link).path
        if path.endswith("/"):
            continue
        if link in seen:
            continue
        seen.add(link)
        out.append(link)
    return out


def select_fernbase_confidence_candidate(label, candidates, confidence):
    desired = str(confidence or "").strip().lower()
    filtered = []
    for candidate in candidates:
        path_lower = urlparse(candidate).path.lower()
        if label == "CDS":
            if not is_probable_cds_url("fernbase", path_lower):
                continue
        elif label == "GFF":
            if not is_probable_gff_url(path_lower):
                continue
            if FERNBASE_GFF_EXCLUDE_PATTERN.search(path_lower):
                continue
        else:
            continue
        if desired == "high" and "highconfidence" not in path_lower:
            continue
        if desired == "low" and "lowconfidence" not in path_lower:
            continue
        filtered.append(candidate)
    if len(filtered) == 0:
        return ""
    return sorted(
        filtered,
        key=lambda x: provider_candidate_sort_key("fernbase", label, urlparse(x).path.split("/")[-1]),
    )[0]


def resolve_fernbase_split_bundle(cds_url, gff_url, genome_url, timeout, headers):
    references = [gff_url, cds_url, genome_url]
    candidates = []
    last_error = None
    for reference in references:
        text = str(reference or "").strip()
        if text == "":
            continue
        try:
            candidates = list_directory_file_urls(text, timeout, headers)
        except Exception as exc:
            last_error = exc
            continue
        if len(candidates) > 0:
            break
    if len(candidates) == 0:
        if last_error is not None:
            raise ValueError(last_error)
        return {"is_split": False}

    high_cds_url = select_fernbase_confidence_candidate("CDS", candidates, "high")
    low_cds_url = select_fernbase_confidence_candidate("CDS", candidates, "low")
    high_gff_url = select_fernbase_confidence_candidate("GFF", candidates, "high")
    low_gff_url = select_fernbase_confidence_candidate("GFF", candidates, "low")
    resolved_genome_url = str(genome_url or "").strip()
    if resolved_genome_url == "":
        resolved_genome_url = select_best_url_for_label("fernbase", "GENOME", candidates)
    return {
        "is_split": all(
            value != ""
            for value in (
                high_cds_url,
                low_cds_url,
                high_gff_url,
                low_gff_url,
            )
        ),
        "high_cds_url": high_cds_url,
        "low_cds_url": low_cds_url,
        "high_gff_url": high_gff_url,
        "low_gff_url": low_gff_url,
        "genome_url": resolved_genome_url,
        "high_cds_filename": Path(urlparse(high_cds_url).path).name if high_cds_url != "" else "",
        "low_cds_filename": Path(urlparse(low_cds_url).path).name if low_cds_url != "" else "",
        "high_gff_filename": Path(urlparse(high_gff_url).path).name if high_gff_url != "" else "",
        "low_gff_filename": Path(urlparse(low_gff_url).path).name if low_gff_url != "" else "",
    }


def interval_overlaps_merged_intervals(intervals, start, end):
    lo = 0
    hi = len(intervals)
    while lo < hi:
        mid = (lo + hi) // 2
        if intervals[mid][1] < start:
            lo = mid + 1
        else:
            hi = mid
    if lo >= len(intervals):
        return False
    return intervals[lo][0] <= end


def load_fernbase_high_gene_intervals(gff_path):
    intervals_by_seqid = defaultdict(list)
    with open_text(gff_path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line == "":
                continue
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or str(parts[2] or "").strip().lower() != "gene":
                continue
            try:
                start = int(parts[3])
                end = int(parts[4])
            except Exception:
                continue
            if start > end:
                start, end = end, start
            intervals_by_seqid[str(parts[0] or "").strip()].append((start, end))
    return {seqid: merge_coordinate_intervals(intervals) for seqid, intervals in intervals_by_seqid.items()}


def select_low_confidence_gene_ids_without_high_overlap(low_gff_path, high_gene_intervals):
    kept_gene_ids = set()
    total_low_genes = 0
    with open_text(low_gff_path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line == "":
                continue
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or str(parts[2] or "").strip().lower() != "gene":
                continue
            try:
                start = int(parts[3])
                end = int(parts[4])
            except Exception:
                continue
            if start > end:
                start, end = end, start
            attrs = parse_gff_attributes(parts[8])
            gene_id = choose_first_gff_attribute(attrs, ("ID", "gene_id", "locus_tag", "gene", "Name"))
            if gene_id == "":
                continue
            total_low_genes += 1
            intervals = high_gene_intervals.get(str(parts[0] or "").strip(), ())
            if interval_overlaps_merged_intervals(intervals, start, end):
                continue
            kept_gene_ids.add(gene_id)
    return kept_gene_ids, total_low_genes


def record_matches_fernbase_gene_ids(record, kept_gene_ids, feature_records, gene_cache):
    attrs = record.get("attrs", {})
    tokens = set()
    explicit = choose_first_gff_attribute(attrs, ("gene", "gene_id", "locus_tag", "geneName"))
    if explicit != "":
        tokens.add(explicit)
    feature_type = str(record.get("feature_type", "") or "").strip().lower()
    if feature_type == "gene":
        direct = choose_first_gff_attribute(attrs, ("ID", "Name"))
        if direct != "":
            tokens.add(direct)
    feature_id = str(record.get("feature_id", "") or "").strip()
    if feature_id != "":
        tokens.add(resolve_feature_gene_token(feature_id, feature_records, "fernbase", gene_cache, set()))
    for parent_id in record.get("parents", ()):
        parent_token = resolve_feature_gene_token(parent_id, feature_records, "fernbase", gene_cache, set())
        if parent_token != "":
            tokens.add(parent_token)
    return any(str(token or "").strip() in kept_gene_ids for token in tokens)


def iter_fernbase_filtered_low_gff_lines(low_gff_path, kept_gene_ids):
    records = []
    feature_records = {}
    with open_text(low_gff_path, "rt") as handle:
        for raw_line in handle:
            if raw_line.startswith("##FASTA"):
                break
            line = raw_line.rstrip("\n\r")
            if line == "" or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            attrs = parse_gff_attributes(parts[8])
            feature_id = choose_first_gff_attribute(attrs, ("ID", "transcript_id", "protein_id", "Name"))
            parents = attrs.get("Parent", ())
            record = {
                "line": raw_line if raw_line.endswith("\n") else raw_line + "\n",
                "feature_type": str(parts[2] or "").strip().lower(),
                "attrs": attrs,
                "feature_id": feature_id,
                "parents": parents,
            }
            records.append(record)
            if feature_id != "":
                feature_records[feature_id] = {
                    "feature_type": record["feature_type"],
                    "parents": parents,
                    "attrs": attrs,
                }
    gene_cache = {}
    for record in records:
        if record_matches_fernbase_gene_ids(record, kept_gene_ids, feature_records, gene_cache):
            yield record["line"]


def iter_fernbase_combined_gff_lines(high_gff_path, low_gff_path, kept_gene_ids):
    with open_text(high_gff_path, "rt") as handle:
        for raw_line in handle:
            if raw_line.startswith("##FASTA"):
                break
            yield raw_line if raw_line.endswith("\n") else raw_line + "\n"
    yield from iter_fernbase_filtered_low_gff_lines(low_gff_path, kept_gene_ids)


def iter_fernbase_combined_cds_records(high_cds_path, low_cds_path, kept_gene_ids):
    yield from iter_fasta_records(high_cds_path)
    for header, sequence in iter_fasta_records(low_cds_path):
        gene_token = fernbase_feature_gene_token_from_header(header)
        if gene_token in kept_gene_ids:
            yield header, sequence


def merge_fernbase_confidence_bundle(
    species_key,
    high_cds_path,
    low_cds_path,
    high_gff_path,
    low_gff_path,
    combined_cds_path,
    combined_gff_path,
    overwrite,
):
    if (
        combined_cds_path.exists()
        and combined_cds_path.stat().st_size > 0
        and combined_gff_path.exists()
        and combined_gff_path.stat().st_size > 0
        and not overwrite
    ):
        return {
            "warnings": [
                "[download:fernbase] {} combined confidence bundle already exists. Skipping merge.".format(species_key)
            ],
            "errors": [],
        }

    high_gene_intervals = load_fernbase_high_gene_intervals(high_gff_path)
    kept_low_gene_ids, total_low_genes = select_low_confidence_gene_ids_without_high_overlap(
        low_gff_path,
        high_gene_intervals,
    )
    write_text_lines(
        combined_gff_path,
        iter_fernbase_combined_gff_lines(high_gff_path, low_gff_path, kept_low_gene_ids),
    )
    write_fasta_records(
        combined_cds_path,
        iter_fernbase_combined_cds_records(high_cds_path, low_cds_path, kept_low_gene_ids),
    )
    return {
        "warnings": [
            "[download:fernbase] {} merged FernBase confidence bundle: kept {} of {} low-confidence genes".format(
                species_key,
                len(kept_low_gene_ids),
                total_low_genes,
            )
        ],
        "errors": [],
    }


def normalize_oryza_minuta_source_id(source_id):
    text = strip_provider_prefix(source_id, ORYZA_MINUTA_PROVIDER)
    normalized = re.sub(r"[^A-Za-z0-9]+", "_", str(text or "").strip().lower()).strip("_")
    if normalized in ORYZA_MINUTA_SOURCE_ID_ALIASES:
        return ORYZA_MINUTA_CANONICAL_SOURCE_ID
    raise ValueError(
        "unsupported id '{}' for provider={}; expected one of: {}".format(
            source_id,
            ORYZA_MINUTA_PROVIDER,
            ", ".join(sorted(value for value in ORYZA_MINUTA_SOURCE_ID_ALIASES if value != "")),
        )
    )


def build_oryza_minuta_merged_filename(label):
    normalized = str(label or "").strip().lower()
    if normalized == "cds":
        return "Oryza_minuta.gramene_tetraploids.merged.cds.all.fa.gz"
    if normalized == "gff":
        return "Oryza_minuta.gramene_tetraploids.merged.gff3.gz"
    if normalized == "genome":
        return "Oryza_minuta.gramene_tetraploids.merged.dna.toplevel.fa.gz"
    raise ValueError("unknown Oryza minuta merged bundle label: {}".format(label))


def resolve_oryza_minuta_download_bundle_from_id(source_id, species_key):
    canonical_id = normalize_oryza_minuta_source_id(source_id)
    resolved_species_key = normalize_species_key_for_runtime(species_key or ORYZA_MINUTA_DEFAULT_SPECIES_KEY)
    base_url = resolve_oryza_minuta_gramene_base_url()
    download_targets = []
    representative = {"cds_url": "", "gff_url": "", "genome_url": ""}

    for bundle in ORYZA_MINUTA_SUBGENOME_BUNDLES:
        token = bundle["token"]
        cds_url = "{}/fasta/{}/cds/{}".format(base_url, token, bundle["cds_filename"])
        gff_url = "{}/gff3/{}/{}".format(base_url, token, bundle["gff_filename"])
        genome_url = "{}/fasta/{}/dna/{}".format(base_url, token, bundle["genome_filename"])
        if representative["cds_url"] == "":
            representative["cds_url"] = cds_url
            representative["gff_url"] = gff_url
            representative["genome_url"] = genome_url
        download_targets.extend(
            (
                {
                    "label": "CDS_{}".format(bundle["label"]),
                    "url": cds_url,
                    "filename": bundle["cds_filename"],
                },
                {
                    "label": "GFF_{}".format(bundle["label"]),
                    "url": gff_url,
                    "filename": bundle["gff_filename"],
                },
                {
                    "label": "GENOME_{}".format(bundle["label"]),
                    "url": genome_url,
                    "filename": bundle["genome_filename"],
                },
            )
        )

    return {
        "canonical_id": canonical_id,
        "species_key": resolved_species_key,
        "cds_url": representative["cds_url"],
        "gff_url": representative["gff_url"],
        "genome_url": representative["genome_url"],
        "cds_filename": build_oryza_minuta_merged_filename("cds"),
        "gff_filename": build_oryza_minuta_merged_filename("gff"),
        "genome_filename": build_oryza_minuta_merged_filename("genome"),
        "download_targets": tuple(download_targets),
    }


def iter_oryza_minuta_combined_gff_lines(gff_paths):
    emitted_version = False
    for path in gff_paths:
        with open_text(path, "rt") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n\r")
                if line == "":
                    continue
                if line.startswith("##FASTA"):
                    break
                if line.startswith("##gff-version"):
                    if emitted_version:
                        continue
                    emitted_version = True
                    yield "##gff-version 3\n"
                    continue
                yield raw_line if raw_line.endswith("\n") else raw_line + "\n"
    if not emitted_version:
        yield "##gff-version 3\n"


def iter_oryza_minuta_combined_fasta_records(paths):
    for path in paths:
        yield from iter_fasta_records(path)


def merge_oryza_minuta_multisource_bundle(
    species_key,
    cds_paths,
    gff_paths,
    genome_paths,
    combined_cds_path,
    combined_gff_path,
    combined_genome_path,
    overwrite,
):
    combined_outputs = (combined_cds_path, combined_gff_path, combined_genome_path)
    if all(path.exists() and path.stat().st_size > 0 for path in combined_outputs) and not overwrite:
        return {
            "warnings": [
                "[download:{}] {} merged Gramene tetraploid bundle already exists. Skipping merge.".format(
                    ORYZA_MINUTA_PROVIDER,
                    species_key,
                )
            ],
            "errors": [],
        }

    write_fasta_records(combined_cds_path, iter_oryza_minuta_combined_fasta_records(cds_paths))
    write_text_lines(combined_gff_path, iter_oryza_minuta_combined_gff_lines(gff_paths))
    write_fasta_records(combined_genome_path, iter_oryza_minuta_combined_fasta_records(genome_paths))
    return {
        "warnings": [
            "[download:{}] {} merged {} public Gramene tetraploid subgenome bundles".format(
                ORYZA_MINUTA_PROVIDER,
                species_key,
                len(tuple(ORYZA_MINUTA_SUBGENOME_BUNDLES)),
            )
        ],
        "errors": [],
    }


def resolve_oryza_minuta_download_urls_from_id(source_id, species_key, timeout, headers):
    del timeout
    del headers
    return resolve_oryza_minuta_download_bundle_from_id(source_id, species_key)


def fernbase_release_sort_key(name):
    lower = str(name or "").lower()
    version_tokens = tuple(int(token) for token in re.findall(r"[0-9]+", lower))
    has_version = 1 if re.search(r"(?:^|[_-])(?:asm[_-])?v?[0-9]+", lower) else 0
    return (has_version, version_tokens, lower)


def infer_fernbase_species_key(source_id, species_key):
    explicit = str(species_key or "").strip()
    if explicit != "":
        return explicit
    raw = strip_provider_prefix(source_id, "fernbase")
    if raw == "":
        raw = str(source_id or "").strip()
    if is_url_like(raw):
        parts = [unquote(part) for part in urlparse(raw).path.split("/") if part]
        if "ftp" in parts:
            ftp_index = parts.index("ftp")
            parts = parts[ftp_index + 1 :]
        if len(parts) > 0:
            return parts[0]
        return ""
    return raw.split("/", 1)[0]


def resolve_fernbase_download_urls_from_id(source_id, species_key, timeout, headers):
    source_clean = str(source_id or "").strip()
    root_template = os.environ.get("GG_FERNBASE_ID_URL_TEMPLATE", "").strip()
    if root_template == "":
        root_template = "https://fernbase.org/ftp/{id}/"

    index_urls = []
    if is_url_like(source_clean):
        index_urls.append(source_clean)
    id_candidates = source_id_candidates("fernbase", source_id, species_key)
    if len(id_candidates) == 0:
        id_candidates = [source_clean]
    for candidate in id_candidates:
        rendered = render_id_url_template(root_template, "fernbase", candidate, species_key)
        if rendered != "":
            index_urls.append(rendered)

    deduped_index_urls = []
    seen_index_urls = set()
    for index_url in index_urls:
        normalized = str(index_url or "").strip()
        if normalized == "":
            continue
        if normalized in seen_index_urls:
            continue
        seen_index_urls.add(normalized)
        deduped_index_urls.append(normalized)

    inferred_species_key = infer_fernbase_species_key(source_id, species_key)
    last_error = None
    for index_url in deduped_index_urls:
        try:
            resolved = resolve_urls_from_index_url("fernbase", index_url, timeout, headers)
        except Exception as exc:
            last_error = exc
            continue
        if resolved.get("gff_url", "") != "" and resolved.get("genome_url", "") != "":
            if inferred_species_key != "":
                resolved["species_key"] = inferred_species_key
            return resolved

        try:
            index_text = fetch_text_with_headers(index_url, timeout, headers)
        except Exception as exc:
            last_error = exc
            continue

        subdir_candidates = []
        for link in parse_links_from_document(index_url, index_text):
            path = urlparse(link).path
            if not path.endswith("/"):
                continue
            subdir_name = unquote(path.rstrip("/").split("/")[-1])
            if subdir_name in ("", "..", ".", inferred_species_key, "chloroplast_genome"):
                continue
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", subdir_name):
                continue
            subdir_candidates.append((fernbase_release_sort_key(subdir_name), link))

        for _sort_key, subdir_url in sorted(subdir_candidates, reverse=True):
            try:
                resolved = resolve_urls_from_index_url("fernbase", subdir_url, timeout, headers)
            except Exception as exc:
                last_error = exc
                continue
            if resolved.get("gff_url", "") != "" and resolved.get("genome_url", "") != "":
                if inferred_species_key != "":
                    resolved["species_key"] = inferred_species_key
                return resolved

    if last_error is not None:
        raise ValueError("FernBase id '{}' did not resolve: {}".format(source_id, last_error))
    raise ValueError("FernBase id '{}' did not resolve to downloadable GFF/genome URLs".format(source_id))
