#!/usr/bin/env python

import argparse
import csv
import gzip
import os
import re


OUTPUT_COLUMNS = [
    "accession",
    "signal_start",
    "signal_end",
    "transmem_aa",
    "transmem_count",
    "transmem_regions",
    "kegg_gene",
    "kegg_orthology",
    "go_ids",
    "go_aspects",
    "go_terms",
    "go_evidence",
    "gene_name_primary",
    "gene_name_synonyms",
    "ec_numbers",
    "subcellular_location",
    "keywords",
    "interpro_ids",
    "pfam_ids",
    "reactome_ids",
    "organism",
    "taxid",
]


def open_text(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def dedup_keep_order(values):
    seen = set()
    out = []
    for value in values:
        if value == "":
            continue
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def parse_region_token(token):
    if token == "":
        return None
    token = token.replace("<", "").replace(">", "").replace("?", "")
    match = re.match(r"^(\d+)\.\.(\d+)$", token)
    if not match:
        return None
    start = int(match.group(1))
    end = int(match.group(2))
    if end < start:
        start, end = end, start
    return start, end


def parse_ft_region(line):
    if len(line) < 22:
        return None
    payload = line[21:].strip()
    if payload == "":
        return None
    token = payload.split()[0]
    return parse_region_token(token)


def clean_subcell_text(text):
    cleaned = re.sub(r"\{[^{}]*\}", "", text)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    cleaned = cleaned.rstrip(".").strip()
    return cleaned


def normalize_organism_text(parts):
    text = " ".join(parts).strip()
    text = re.sub(r"\s+", " ", text)
    text = text.rstrip(".")
    return text


def build_record_for_accession(
    accession,
    signal_regions,
    tm_regions,
    kegg_gene_ids,
    kegg_orthology_ids,
    go_ids,
    go_aspects,
    go_terms,
    go_evidence,
    gene_name_primary,
    gene_name_synonyms,
    ec_numbers,
    subcellular_location,
    keywords,
    interpro_ids,
    pfam_ids,
    reactome_ids,
    organism,
    taxid,
):
    signal_starts = [str(start) for start, _ in signal_regions]
    signal_ends = [str(end) for _, end in signal_regions]
    tm_segments = [f"{start}-{end}" for start, end in tm_regions]
    tm_total_len = sum((end - start + 1) for start, end in tm_regions)

    return {
        "accession": accession,
        "signal_start": ";".join(signal_starts),
        "signal_end": ";".join(signal_ends),
        "transmem_aa": str(tm_total_len) if tm_regions else "",
        "transmem_count": str(len(tm_regions)) if tm_regions else "",
        "transmem_regions": ";".join(tm_segments),
        "kegg_gene": ";".join(dedup_keep_order(kegg_gene_ids)),
        "kegg_orthology": ";".join(dedup_keep_order(kegg_orthology_ids)),
        "go_ids": ";".join(dedup_keep_order(go_ids)),
        "go_aspects": ";".join(dedup_keep_order(go_aspects)),
        "go_terms": "; ".join(dedup_keep_order(go_terms)),
        "go_evidence": ";".join(dedup_keep_order(go_evidence)),
        "gene_name_primary": gene_name_primary,
        "gene_name_synonyms": ";".join(dedup_keep_order(gene_name_synonyms)),
        "ec_numbers": ";".join(dedup_keep_order(ec_numbers)),
        "subcellular_location": "; ".join(dedup_keep_order(subcellular_location)),
        "keywords": ";".join(dedup_keep_order(keywords)),
        "interpro_ids": ";".join(dedup_keep_order(interpro_ids)),
        "pfam_ids": ";".join(dedup_keep_order(pfam_ids)),
        "reactome_ids": ";".join(dedup_keep_order(reactome_ids)),
        "organism": organism,
        "taxid": taxid,
    }


def emit_records(
    writer,
    accessions,
    signal_regions,
    tm_regions,
    kegg_gene_ids,
    kegg_orthology_ids,
    go_ids,
    go_aspects,
    go_terms,
    go_evidence,
    gene_name_primary,
    gene_name_synonyms,
    ec_numbers,
    subcellular_location,
    keywords,
    interpro_ids,
    pfam_ids,
    reactome_ids,
    organism,
    taxid,
):
    accessions = dedup_keep_order(accessions)
    if len(accessions) == 0:
        return
    signal_regions = dedup_keep_order(signal_regions)
    tm_regions = dedup_keep_order(tm_regions)
    for accession in accessions:
        record = build_record_for_accession(
            accession=accession,
            signal_regions=signal_regions,
            tm_regions=tm_regions,
            kegg_gene_ids=kegg_gene_ids,
            kegg_orthology_ids=kegg_orthology_ids,
            go_ids=go_ids,
            go_aspects=go_aspects,
            go_terms=go_terms,
            go_evidence=go_evidence,
            gene_name_primary=gene_name_primary,
            gene_name_synonyms=gene_name_synonyms,
            ec_numbers=ec_numbers,
            subcellular_location=subcellular_location,
            keywords=keywords,
            interpro_ids=interpro_ids,
            pfam_ids=pfam_ids,
            reactome_ids=reactome_ids,
            organism=organism,
            taxid=taxid,
        )
        writer.writerow(record)


def parse_uniprot_dat(dat_path, out_path):
    with open_text(dat_path, "rt") as fin, open_text(out_path, "wt") as fout:
        writer = csv.DictWriter(fout, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        accessions = []
        signal_regions = []
        tm_regions = []
        kegg_gene_ids = []
        kegg_orthology_ids = []
        go_ids = []
        go_aspects = []
        go_terms = []
        go_evidence = []
        gene_name_primary = ""
        gene_name_synonyms = []
        ec_numbers = []
        subcellular_location = []
        keywords = []
        interpro_ids = []
        pfam_ids = []
        reactome_ids = []
        organism_parts = []
        taxid = ""
        in_subcellular_location = False

        for raw_line in fin:
            line = raw_line.rstrip("\n")
            if line == "//":
                emit_records(
                    writer=writer,
                    accessions=accessions,
                    signal_regions=signal_regions,
                    tm_regions=tm_regions,
                    kegg_gene_ids=kegg_gene_ids,
                    kegg_orthology_ids=kegg_orthology_ids,
                    go_ids=go_ids,
                    go_aspects=go_aspects,
                    go_terms=go_terms,
                    go_evidence=go_evidence,
                    gene_name_primary=gene_name_primary,
                    gene_name_synonyms=gene_name_synonyms,
                    ec_numbers=ec_numbers,
                    subcellular_location=subcellular_location,
                    keywords=keywords,
                    interpro_ids=interpro_ids,
                    pfam_ids=pfam_ids,
                    reactome_ids=reactome_ids,
                    organism=normalize_organism_text(organism_parts),
                    taxid=taxid,
                )
                accessions = []
                signal_regions = []
                tm_regions = []
                kegg_gene_ids = []
                kegg_orthology_ids = []
                go_ids = []
                go_aspects = []
                go_terms = []
                go_evidence = []
                gene_name_primary = ""
                gene_name_synonyms = []
                ec_numbers = []
                subcellular_location = []
                keywords = []
                interpro_ids = []
                pfam_ids = []
                reactome_ids = []
                organism_parts = []
                taxid = ""
                in_subcellular_location = False
                continue

            if in_subcellular_location:
                if line.startswith("CC       "):
                    cleaned = clean_subcell_text(line[5:].strip())
                    if cleaned != "":
                        subcellular_location.append(cleaned)
                    continue
                in_subcellular_location = False

            if line.startswith("AC   "):
                for token in line[5:].split(";"):
                    token = token.strip()
                    if token != "":
                        accessions.append(token)
                continue

            if line.startswith("GN   "):
                payload = line[5:].strip()
                for match in re.finditer(r"Name=([^;]+);", payload):
                    for token in match.group(1).split(","):
                        token = token.strip()
                        if token == "":
                            continue
                        if gene_name_primary == "":
                            gene_name_primary = token
                        elif token != gene_name_primary:
                            gene_name_synonyms.append(token)
                for match in re.finditer(r"Synonyms=([^;]+);", payload):
                    for token in match.group(1).split(","):
                        token = token.strip()
                        if token != "" and token != gene_name_primary:
                            gene_name_synonyms.append(token)
                continue

            if line.startswith("DE   "):
                for ec in re.findall(r"EC=([^;]+);", line):
                    ec = ec.strip()
                    if ec != "":
                        ec_numbers.append(ec)
                continue

            if line.startswith("OS   "):
                text = line[5:].strip()
                if text != "":
                    organism_parts.append(text)
                continue

            if line.startswith("OX   "):
                match = re.search(r"NCBI_TaxID=([0-9]+)", line)
                if match and taxid == "":
                    taxid = match.group(1)
                continue

            if line.startswith("KW   "):
                for token in line[5:].split(";"):
                    token = token.strip().rstrip(".")
                    if token != "":
                        keywords.append(token)
                continue

            if line.startswith("CC   -!- SUBCELLULAR LOCATION:"):
                in_subcellular_location = True
                cleaned = clean_subcell_text(line.split("SUBCELLULAR LOCATION:", 1)[1].strip())
                if cleaned != "":
                    subcellular_location.append(cleaned)
                continue

            if line.startswith("DR   GO;"):
                parts = [p.strip() for p in line[5:].split(";")]
                if len(parts) >= 2 and parts[1] != "":
                    go_ids.append(parts[1])
                if len(parts) >= 3 and parts[2] != "":
                    aspect_term = parts[2]
                    if ":" in aspect_term:
                        aspect, term = aspect_term.split(":", 1)
                        if aspect != "":
                            go_aspects.append(aspect)
                        if term != "":
                            go_terms.append(term)
                if len(parts) >= 4 and parts[3] != "":
                    evidence = parts[3].split(":", 1)[0].strip().rstrip(".")
                    if evidence != "":
                        go_evidence.append(evidence)
                continue

            if line.startswith("DR   KEGG;"):
                parts = [p.strip() for p in line[5:].split(";")]
                if len(parts) >= 2 and parts[1] != "":
                    kegg_gene_ids.append(parts[1])
                continue

            if line.startswith("DR   KO;"):
                parts = [p.strip() for p in line[5:].split(";")]
                if len(parts) >= 2 and parts[1] != "":
                    kegg_orthology_ids.append(parts[1])
                continue

            if line.startswith("DR   InterPro;"):
                parts = [p.strip() for p in line[5:].split(";")]
                if len(parts) >= 2 and parts[1] != "":
                    interpro_ids.append(parts[1])
                continue

            if line.startswith("DR   Pfam;"):
                parts = [p.strip() for p in line[5:].split(";")]
                if len(parts) >= 2 and parts[1] != "":
                    pfam_ids.append(parts[1])
                continue

            if line.startswith("DR   Reactome;"):
                parts = [p.strip() for p in line[5:].split(";")]
                if len(parts) >= 2 and parts[1] != "":
                    reactome_ids.append(parts[1])
                continue

            if line.startswith("FT   SIGNAL"):
                region = parse_ft_region(line)
                if region is not None:
                    signal_regions.append(region)
                continue

            if line.startswith("FT   TRANSMEM"):
                region = parse_ft_region(line)
                if region is not None:
                    tm_regions.append(region)
                continue

        # Last record if no trailing '//'
        emit_records(
            writer=writer,
            accessions=accessions,
            signal_regions=signal_regions,
            tm_regions=tm_regions,
            kegg_gene_ids=kegg_gene_ids,
            kegg_orthology_ids=kegg_orthology_ids,
            go_ids=go_ids,
            go_aspects=go_aspects,
            go_terms=go_terms,
            go_evidence=go_evidence,
            gene_name_primary=gene_name_primary,
            gene_name_synonyms=gene_name_synonyms,
            ec_numbers=ec_numbers,
            subcellular_location=subcellular_location,
            keywords=keywords,
            interpro_ids=interpro_ids,
            pfam_ids=pfam_ids,
            reactome_ids=reactome_ids,
            organism=normalize_organism_text(organism_parts),
            taxid=taxid,
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--uniprot_dat", metavar="PATH", required=True, type=str, help="UniProt Swiss-Prot DAT(.gz) path.")
    parser.add_argument("--outfile", metavar="PATH", required=True, type=str, help="Output metadata TSV(.gz).")
    args = parser.parse_args()

    if (not os.path.exists(args.uniprot_dat)) or (os.path.getsize(args.uniprot_dat) == 0):
        raise SystemExit(f"Input DAT file was not found or empty: {args.uniprot_dat}")

    parse_uniprot_dat(dat_path=args.uniprot_dat, out_path=args.outfile)
    if (not os.path.exists(args.outfile)) or (os.path.getsize(args.outfile) == 0):
        raise SystemExit(f"Output metadata file was not created: {args.outfile}")


if __name__ == "__main__":
    main()
