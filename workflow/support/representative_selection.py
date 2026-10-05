"""Read immutable locus-to-transcript choices for downstream annotation readers.

Exact identifiers win over documented GFF/species aliases. Alias ambiguity and
missing choices are errors: a supplied manifest never reselects longest CDS.
"""

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

try:
    from format_species_annotation.grouping_identity import strip_gff_feature_prefix
    from species_labeling import extract_species_label
except ImportError:  # package imports in tests
    from .format_species_annotation.grouping_identity import strip_gff_feature_prefix
    from .species_labeling import extract_species_label

REQUIRED_COLUMNS = (
    "species", "gene_id", "candidate_id", "source_transcript_id", "status",
    "score", "margin", "reason",
)


def identifier_aliases(identifier, species=""):
    """Only structural prefixes and the exact declared species may be removed."""
    identifier = str(identifier)
    aliases = {identifier, strip_gff_feature_prefix(identifier)}
    if species:
        aliases.update(value.removeprefix(species + "_") for value in tuple(aliases))
    aliases.update(strip_gff_feature_prefix(value) for value in tuple(aliases))
    return frozenset(aliases - {""})


def matching_identifier(identifier, choices, species=""):
    """Return a unique actual identifier, prioritising an exact source ID."""
    if identifier in choices:
        return identifier
    aliases = identifier_aliases(identifier, species)
    matches = [value for value in choices if aliases & identifier_aliases(value, species)]
    if len(matches) > 1:
        raise ValueError(f"Ambiguous representative transcript alias {identifier}: {sorted(matches)}")
    return matches[0] if matches else None


class RepresentativeMap:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.rows = {}
        self.by_species = {}
        self.alias_rows = {}
        candidate_owners, transcript_owners = {}, {}
        with self.path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            if len(set(fields)) != len(fields) or set(REQUIRED_COLUMNS) - set(fields):
                raise ValueError("Representative map requires unique columns: " + ",".join(REQUIRED_COLUMNS))
            for line, row in enumerate(reader, 2):
                if None in row or any(row[key] is None for key in REQUIRED_COLUMNS):
                    raise ValueError(f"Malformed representative map row {line}")
                for key in ("species", "gene_id", "candidate_id", "source_transcript_id", "status"):
                    if not row[key] or any(char.isspace() for char in row[key]):
                        raise ValueError(f"Invalid representative map {key} at row {line}")
                for key in ("score", "margin"):
                    if row[key]:
                        try:
                            value = float(row[key])
                        except ValueError as exc:
                            raise ValueError(f"Invalid representative map {key} at row {line}") from exc
                        if not math.isfinite(value):
                            raise ValueError(f"Nonfinite representative map {key} at row {line}")
                key = row["species"], row["gene_id"]
                if key in self.rows:
                    raise ValueError(f"Duplicate representative choice for {key[0]}:{key[1]}")
                for column, owners in (("candidate_id", candidate_owners),
                                       ("source_transcript_id", transcript_owners)):
                    identity = row["species"], row[column]
                    if identity in owners:
                        raise ValueError(f"Conflicting representative {column} ownership: {row[column]}")
                    owners[identity] = row["gene_id"]
                self.rows[key] = row
                self.by_species.setdefault(key[0], {})[key[1]] = row
                for alias in identifier_aliases(key[1], key[0]):
                    self.alias_rows.setdefault((key[0], alias), []).append(row)
        if not self.rows:
            raise ValueError("Empty representative map")

    def choice(self, gene_id, species=""):
        if not species:
            inferred = extract_species_label(gene_id)
            species = inferred if inferred in self.by_species else ""
        groups = ([self.by_species.get(species, {})] if species else self.by_species.values())
        matches = []
        for rows in groups:
            if not rows:
                continue
            if gene_id in rows:
                matches.append(rows[gene_id])
                continue
            aliases = identifier_aliases(gene_id, species)
            group_species = species or next(iter(rows.values()))["species"]
            group_matches = {}
            for alias in aliases:
                for row in self.alias_rows.get((group_species, alias), ()):
                    group_matches[row["gene_id"]] = row
            matches.extend(group_matches.values())
        if len(matches) != 1:
            detail = "ambiguous" if matches else "missing"
            raise ValueError(f"Representative map {detail} choice for {species}:{gene_id}")
        return matches[0]


def load_representative_map(value):
    if isinstance(value, RepresentativeMap):
        return value
    # Public package imports and executable support-script imports can coexist
    # in one interpreter. Preserve an already validated map from this exact
    # module file rather than interpreting the equivalent class as a path.
    module = sys.modules.get(type(value).__module__)
    other_class = getattr(module, 'RepresentativeMap', None)
    if (other_class is not None and isinstance(value, other_class)
            and Path(module.__file__).resolve() == Path(__file__).resolve()):
        return value
    return RepresentativeMap(value) if value else None


EFFECTIVE_ROLES = ("cds", "protein", "gff", "genome", "representative_map")
ANALYSIS_ROLES = ('analysis_cds', 'analysis_gff')


def file_digest(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def verify_published_inputs(path, rows, source_hashes):
    """Bind copied manifests to a published view while allowing external TSVs."""
    roots = set()
    for row in rows.values():
        for role in EFFECTIVE_ROLES:
            source = Path(row[role])
            root = (source.parent if role == 'representative_map' else source.parent.parent)
            canonical = (source.name == 'representative_map.tsv' if role == 'representative_map'
                         else source.parent.name == 'species_' + role)
            published = (root / 'receipt.json').exists()
            if canonical and not published:
                try:
                    frozen = json.loads((root.parent / 'plan.json').read_text())
                    published = bool(frozen['request']['implementation'].get('gene_model_refinement.py'))
                except (OSError, ValueError, TypeError, KeyError):
                    pass
            if canonical and published:
                roots.add(root)
    if not roots:
        return
    if len(roots) != 1:
        raise ValueError('Effective inputs reference multiple published views')
    root = roots.pop()
    try:
        receipt = json.loads((root / 'receipt.json').read_text())
        files = receipt['files']
        if not isinstance(files, dict) or not files or not receipt.get('key'):
            raise ValueError('Missing publication contract')
        for relative, expected in files.items():
            source = (root / relative).resolve()
            if not source.is_relative_to(root) or not source.is_file():
                raise ValueError('Invalid publication file')
            if source not in source_hashes:
                source_hashes[source] = file_digest(source)
            if source_hashes[source] != expected:
                raise ValueError('Changed publication file')
        if file_digest(path) != files.get('inputs.tsv'):
            raise ValueError('Effective input manifest differs from the published view')
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError('Effective view is incomplete or corrupt') from exc


def load_effective_inputs(path):
    """Verify the selected input bundle as one manifest-bound collection."""
    path = Path(path).resolve()
    required = {"species", "genetic_code", *EFFECTIVE_ROLES,
                *(role + "_sha256" for role in EFFECTIVE_ROLES)}
    result = {}
    maps = {}
    source_hashes = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or required - set(fields):
            raise ValueError("Effective inputs manifest missing required paths, hashes, or species columns")
        analysis_fields = {*ANALYSIS_ROLES, *(role + '_sha256' for role in ANALYSIS_ROLES)}
        has_analysis = bool(analysis_fields & set(fields))
        if has_analysis and analysis_fields - set(fields):
            raise ValueError('Effective analysis view requires CDS/GFF paths and both hashes')
        if has_analysis:
            required.update(analysis_fields)
        for row in reader:
            if None in row or any(row[key] is None or not row[key] for key in required):
                raise ValueError("Malformed effective inputs manifest row")
            species = row["species"]
            if species in result or any(char.isspace() for char in species):
                raise ValueError(f"Invalid or duplicate effective inputs species: {species}")
            try:
                int(row["genetic_code"])
            except ValueError as exc:
                raise ValueError(f"Invalid effective inputs genetic code for {species}") from exc
            for role in (*EFFECTIVE_ROLES, *(ANALYSIS_ROLES if has_analysis else ())):
                source = Path(row[role])
                source = (source if source.is_absolute() else path.parent / source).resolve()
                if source not in source_hashes:
                    source_hashes[source] = file_digest(source)
                if source_hashes[source] != row[role + "_sha256"]:
                    raise ValueError(f"Changed effective inputs {role} for {species}: {source}")
                row[role] = str(source)
            selection_path = row["representative_map"]
            if selection_path not in maps:
                maps[selection_path] = RepresentativeMap(selection_path)
            if species not in maps[selection_path].by_species:
                raise ValueError(f"Effective representative map lacks species: {species}")
            result[species] = row
    if not result:
        raise ValueError("Empty effective inputs manifest")
    verify_published_inputs(path, result, source_hashes)
    return result
