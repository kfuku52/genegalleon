"""Conservative, copy-aware selection of coding isoforms on a frozen synteny graph.

The scores and default cutoffs are experimental ranking rules, not calibrated
probabilities.  This module selects one record per input locus and never changes
the correspondence graph or annotation.  Ambiguous or weak decisions retain the
source representative and expose the reason in the result.
"""

from __future__ import annotations

import itertools
import math
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from typing import Any

from Bio.Align import PairwiseAligner

HARD_FAILURES = ("internal_stop", "phase_conflict", "sequence_mismatch", "frameshift", "ownership_conflict",
                 "translation_uncertain")
INDEPENDENT_SUPPORT = ("rna_supported", "full_length_supported", "principal_supported", "domain_supported")


def _protein(candidate: dict) -> str:
    return str(candidate.get("protein", "")).upper().removesuffix("*")


def _length(candidate: dict) -> int:
    """Use the formatter's corrected, unpadded length when it is supplied."""
    return int(candidate.get("corrected_cds_length", len(candidate.get("cds", "")) or 3 * len(_protein(candidate))))


def _independent_support(candidate: dict) -> bool:
    quality = candidate.get("quality", {})
    return any(quality.get(key, False) for key in INDEPENDENT_SUPPORT) or float(quality.get("independent_support", 0)) > 0


def _eligible(candidate: dict) -> bool:
    quality = candidate.get("quality", {})
    protein = _protein(candidate)
    if not protein or "*" in protein or any(quality.get(key, False) for key in HARD_FAILURES):
        return False
    if not quality.get("usable", True):
        return False
    # Accepting a homology-only prediction into a catalog is not permission to
    # replace the representative.  The caller owns this independent gate.
    if candidate.get("origin", "original") != "original" and not quality.get("representative_eligible", False):
        return False
    return True


def _junction_positions(candidate: dict) -> tuple[tuple[int, int], ...]:
    """Project coding junctions to protein coordinates, independent of chromosome.

    GFF blocks are 0-based half-open.  Their first phase is trimmed only at the
    partial transcript start; later phases describe continuation of a codon and
    must not be subtracted from each exon.  Positions include the split-codon
    remainder, so phase differences cannot masquerade as conserved junctions.
    """
    if "protein_junctions" in candidate:
        return tuple(sorted((int(position), int(phase)) for position, phase in candidate["protein_junctions"]))
    blocks = candidate.get("blocks", [])
    if not blocks:
        return ()
    ordered = sorted(blocks, key=lambda block: (int(block[0]), int(block[1])),
                     reverse=candidate.get("strand", "+") == "-")
    offset = -int(ordered[0][2] or 0)
    positions = []
    for block in ordered[:-1]:
        offset += int(block[1]) - int(block[0])
        if offset >= 0:
            positions.append(divmod(offset, 3))
    return tuple(positions)


def _quality(candidate: dict, longest: int) -> float:
    quality = candidate.get("quality", {})
    completeness = min(1.0, _length(candidate) / max(1, longest))
    result = 0.30 * completeness
    if quality.get("valid_orf", True) and not quality.get("partial", False):
        result += 0.20
    else:
        result -= 0.20
    ambiguous = quality.get("ambiguous", False)
    if ambiguous:
        result -= 0.20
    protein = _protein(candidate)
    result -= 0.50 * (protein.count("X") / max(1, len(protein)))
    if _independent_support(candidate):
        result += 0.25
    return result


@dataclass(frozen=True)
class PairScore:
    score: float
    identity: float
    coverage_a: float
    coverage_b: float
    junction_similarity: float
    aligned_residues: int
    internal_loss_a: bool
    internal_loss_b: bool
    bounded: bool = False

    def as_dict(self) -> dict:
        return dict(self.__dict__)


def _aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    return aligner


def pair_score(candidate_a: dict, candidate_b: dict, *, max_alignment_cells: int = 25_000_000) -> PairScore:
    """Global protein alignment with reciprocal coverage and projected junctions."""
    protein_a, protein_b = _protein(candidate_a), _protein(candidate_b)
    if not protein_a or not protein_b or len(protein_a) * len(protein_b) > max_alignment_cells:
        return PairScore(0.0, 0.0, 0.0, 0.0, 0.0, 0, False, False, True)
    return _score_signatures((protein_a, _junction_positions(candidate_a)),
                             (protein_b, _junction_positions(candidate_b)),
                             max_alignment_cells=max_alignment_cells, aligner_factory=_aligner)


def _project_junctions(junctions: tuple, segments: list, protein_length: int, columns: int) -> set:
    """Project only requested boundaries through monotone alignment segments."""
    positions, index = set(), 0
    for residue, phase in sorted(junctions):
        if residue == protein_length:
            positions.add((columns, phase))
        elif 0 <= residue < protein_length:
            while index < len(segments) and residue >= segments[index][1]:
                index += 1
            if index < len(segments):
                start, end, column = segments[index]
                if start <= residue < end:
                    positions.add((column + residue - start, phase))
    return positions


def _score_signatures(signature_a: tuple, signature_b: tuple, *, max_alignment_cells: int,
                      aligner_factory) -> PairScore:
    protein_a, boundaries_a = signature_a
    protein_b, boundaries_b = signature_b
    if not protein_a or not protein_b or len(protein_a) * len(protein_b) > max_alignment_cells:
        return PairScore(0.0, 0.0, 0.0, 0.0, 0.0, 0, False, False, True)
    aligned = matches = column = 0
    internal_loss_a = internal_loss_b = False
    segments_a, segments_b = [], []
    if protein_a == protein_b and protein_a.isascii():
        # A full identical alignment has no gaps under these scoring rules.
        # Keep the native input contract for non-ASCII strings and apply the
        # cell cap before this shortcut, including perfectly matching pairs.
        aligned = column = len(protein_a)
        matches = sum(residue not in "XBZJ" for residue in protein_a)
        if boundaries_a:
            segments_a.append((0, aligned, 0))
        if boundaries_b:
            segments_b.append((0, aligned, 0))
    else:
        coordinates = aligner_factory().align(protein_a, protein_b)[0].coordinates
        for index in range(coordinates.shape[1] - 1):
            start_a, end_a = int(coordinates[0, index]), int(coordinates[0, index + 1])
            start_b, end_b = int(coordinates[1, index]), int(coordinates[1, index + 1])
            width_a, width_b = end_a - start_a, end_b - start_b
            width = max(width_a, width_b)
            if width_a and width_b:
                aligned += width_a
                matches += sum(a == b and a not in "XBZJ" for a, b in
                               zip(protein_a[start_a:end_a], protein_b[start_b:end_b], strict=True))
            elif width_a >= 9 and 0 < start_b < len(protein_b):
                internal_loss_a = True
            elif width_b >= 9 and 0 < start_a < len(protein_a):
                internal_loss_b = True
            if boundaries_a and width_a:
                segments_a.append((start_a, end_a, column))
            if boundaries_b and width_b:
                segments_b.append((start_b, end_b, column))
            column += width
    junctions_a = _project_junctions(boundaries_a, segments_a, len(protein_a), column)
    junctions_b = _project_junctions(boundaries_b, segments_b, len(protein_b), column)
    union = junctions_a | junctions_b
    junction_similarity = len(junctions_a & junctions_b) / len(union) if union else 1.0
    identity = matches / max(1, aligned)
    coverage_a, coverage_b = aligned / len(protein_a), aligned / len(protein_b)
    # The identity contribution is also multiplied by reciprocal coverage:
    # a perfect short overlap cannot receive the same score as a complete match.
    reciprocal = min(coverage_a, coverage_b)
    score = reciprocal * (0.85 * identity + 0.15 * junction_similarity)
    return PairScore(score, identity, coverage_a, coverage_b, junction_similarity, aligned,
                     internal_loss_a, internal_loss_b)


class _PairScores:
    def __init__(self, max_alignment_cells: int, cached: bool, cache_limit: int):
        self.max_alignment_cells = max_alignment_cells
        self.cached = cached
        self.cache_limit = cache_limit
        self.memo: OrderedDict[tuple, PairScore] = OrderedDict()
        self.signatures: dict[int, tuple] = {}
        self.aligner = None
        self.requests = self.evaluations = self.hits = self.bounded = 0
        self.peak_entries = self.evictions = 0

    def signature(self, candidate: dict) -> tuple:
        identifier = id(candidate)
        if identifier not in self.signatures:
            self.signatures[identifier] = (_protein(candidate), _junction_positions(candidate))
        return self.signatures[identifier]

    def get_aligner(self) -> PairwiseAligner:
        if self.aligner is None:
            self.aligner = _aligner()
        return self.aligner

    def get(self, a: dict, b: dict) -> PairScore:
        self.requests += 1
        signature_a, signature_b = self.signature(a), self.signature(b)
        reverse = signature_a > signature_b
        key = (signature_b, signature_a) if reverse else (signature_a, signature_b)
        if self.cached and key in self.memo:
            self.hits += 1
            result = self.memo[key]
            self.memo.move_to_end(key)
        else:
            result = _score_signatures(*key, max_alignment_cells=self.max_alignment_cells,
                                       aligner_factory=self.get_aligner)
            self.evaluations += 1
            self.bounded += int(result.bounded)
            if self.cached:
                self.memo[key] = result
                if len(self.memo) > self.cache_limit:
                    self.memo.popitem(last=False)
                    self.evictions += 1
                self.peak_entries = max(self.peak_entries, len(self.memo))
        if reverse:
            return PairScore(result.score, result.identity, result.coverage_b, result.coverage_a,
                             result.junction_similarity, result.aligned_residues,
                             result.internal_loss_b, result.internal_loss_a, result.bounded)
        return result


def _candidate_id(candidate: dict) -> str:
    return str(candidate["candidate_id"])


def _node_key(species: Any, gene: Any) -> tuple[str, str]:
    return str(species), str(gene)


def _baseline(locus: dict, candidates: list[dict]) -> dict:
    requested = locus.get("baseline_candidate_id") or locus.get("source_baseline_candidate_id")
    by_id = {_candidate_id(candidate): candidate for candidate in candidates}
    if requested:
        if requested not in by_id:
            raise ValueError(f"Baseline candidate {requested!r} is absent from locus {locus['gene_id']!r}")
        return by_id[requested]
    originals = [candidate for candidate in candidates if candidate.get("origin", "original") == "original"]
    return min(originals or candidates, key=lambda candidate: (-_length(candidate), _candidate_id(candidate)))


def _selection(candidate: dict, node: tuple[str, str], status: str, reason: str,
               score: float | None = None, margin: float | None = None) -> dict:
    return {"species": node[0], "gene_id": node[1], "candidate_id": _candidate_id(candidate),
            "source_transcript_id": candidate.get("source_transcript_id", _candidate_id(candidate)),
            "status": status, "score": score, "margin": margin, "reason": reason}


def locus_can_vote(locus: dict, candidate_limit: int = 32) -> bool:
    """Admission predicate shared by in-memory and streaming graph selection."""
    records = locus.get("candidates", [])
    return len(records) <= candidate_limit and any(_eligible(candidate) for candidate in records)


def prepare_correspondence(node_keys, edges, *, coordinate_ambiguity=(), nonvoting=(), copy_ambiguity=()) -> dict:
    """Resolve ambiguity on the complete graph before removing unusable votes.

    Only identifiers, weights and boolean admission metadata are required.  In
    particular, a streaming caller must perform the one-to-many test before
    splitting components, so an omitted edge cannot hide a second gene copy.
    """
    nodes = set(node_keys)
    coordinates = set(coordinate_ambiguity) & nodes
    ambiguity = coordinates | (set(copy_ambiguity) & nodes)
    nonvoting = set(nonvoting) & nodes
    raw_pairs, omitted = {}, {}
    for edge in edges:
        a = _node_key(edge["species_a"], edge["gene_a"])
        b = _node_key(edge["species_b"], edge["gene_b"])
        if a not in nodes or b not in nodes:
            raise ValueError(f"Correspondence edge refers to an absent locus: {a}, {b}")
        if a[0] == b[0]:
            raise ValueError(f"Crossspecies correspondence requires different species: {a}, {b}")
        weight = float(edge.get("weight", 1.0))
        if not math.isfinite(weight) or weight <= 0:
            raise ValueError("Correspondence weights must be finite and positive")
        pair = tuple(sorted((a, b)))
        if edge.get("ambiguous", False):
            ambiguity.update((a, b))
            omitted[pair] = "ambiguous_copy_correspondence"
        else:
            raw_pairs[pair] = max(weight, raw_pairs.get(pair, 0.0))
    copies = defaultdict(lambda: defaultdict(set))
    for a, b in raw_pairs:
        copies[a][b[0]].add(b[1])
        copies[b][a[0]].add(a[1])
    ambiguity.update(node for node, species_copies in copies.items()
                     if any(len(genes) > 1 for genes in species_copies.values()))
    retained = {}
    for (a, b), weight in raw_pairs.items():
        if a in coordinates or b in coordinates:
            omitted[(a, b)] = "ambiguous_locus_coordinates"
        elif a in ambiguity or b in ambiguity:
            omitted[(a, b)] = "ambiguous_copy_correspondence"
        elif a in nonvoting or b in nonvoting:
            omitted[(a, b)] = "ineligible_donor_candidate"
        else:
            retained[(a, b)] = weight
    return {"edge_by_pair": retained, "ambiguity": ambiguity, "coordinate_ambiguity": coordinates,
            "omitted_edges": [{"a": list(a), "b": list(b), "reason": reason}
                              for (a, b), reason in sorted(omitted.items())]}


def select_representatives(catalogs: list[dict], edges: list[dict], policy: str = "conserved", *,
                           min_margin: float = 0.10, min_support: int = 2, min_identity: float = 0.35,
                           min_coverage: float = 0.75, candidate_limit: int = 32,
                           max_alignment_cells: int = 25_000_000, cache_pair_scores: bool = True,
                           pair_cache_limit: int = 50_000, exact_limit: int = 0, _copy_ambiguity=()) -> dict:
    """Select one representative per locus, with deterministic abstention.

    ``exact_limit`` optionally enumerates small component state spaces for
    validation.  The production default is a bounded multi-seed coordinate
    heuristic and does not claim a global optimum.  ``cache_pair_scores=False``
    is an equivalent reference implementation for performance comparisons.
    """
    if policy not in {"longest", "conserved"}:
        raise ValueError(f"Unknown representative policy: {policy}")
    if not math.isfinite(min_margin) or not 0 <= min_margin or not 0 <= min_identity <= 1 or not 0 <= min_coverage <= 1:
        raise ValueError("Selection margin, identity or coverage is out of range")
    if min_support < 1 or candidate_limit < 1 or max_alignment_cells < 1 or pair_cache_limit < 1 or exact_limit < 0:
        raise ValueError("Selection limits must be positive (exact_limit may be zero)")
    loci, candidates, baselines, quality, options, owners = {}, {}, {}, {}, {}, {}
    eligibility = {}

    def is_eligible(candidate: dict) -> bool:
        return eligibility[id(candidate)]

    for catalog in catalogs:
        if catalog.get("schema", 1) != 1:
            raise ValueError("Unsupported gene model catalog schema")
        for locus in catalog.get("loci", []):
            node = _node_key(locus.get("species", catalog.get("species", "")), locus["gene_id"])
            if node in loci:
                raise ValueError(f"Duplicate locus: {node}")
            records = sorted(locus.get("candidates", []), key=_candidate_id)
            if not node[0] or not records:
                raise ValueError(f"Locus requires species and at least one candidate: {node}")
            if len({_candidate_id(candidate) for candidate in records}) != len(records):
                raise ValueError(f"Duplicate candidate ID at locus: {node}")
            for candidate in records:
                identifier = node[0], _candidate_id(candidate)
                if identifier in owners:
                    raise ValueError(f"Candidate has multiple owning loci: {identifier}")
                owners[identifier] = node
                if id(candidate) not in eligibility:
                    eligibility[id(candidate)] = _eligible(candidate)
            loci[node], candidates[node], baselines[node] = locus, records, _baseline(locus, records)
            eligible = [candidate for candidate in records if is_eligible(candidate)]
            # Ineligible coding paths are retained for source provenance and
            # baseline comparisons, but cannot change the quality scale of
            # usable alternatives.  Compare an ineligible source on that same
            # scale, retaining its actual length and other quality penalties.
            longest = max((_length(candidate) for candidate in eligible), default=max(1, _length(baselines[node])))
            quality[node] = {_candidate_id(candidate): _quality(candidate, longest) for candidate in records}
            options[node] = eligible
            if len(records) > candidate_limit:
                options[node] = []
    nodes = sorted(loci)
    if policy == "longest":
        # Preserve the legacy longest contract for source annotations.  New
        # predictions additionally require explicit representative permission;
        # catalog acceptance alone only allows retaining an alternate path.
        longest_options = {node: [candidate for candidate in candidates[node]
                                 if candidate.get("origin", "original") == "original" or is_eligible(candidate)]
                           for node in nodes}
        return {"schema": 1, "policy": policy, "selections": [
            _selection(min(longest_options[node] or [baselines[node]],
                           key=lambda candidate: (-_length(candidate), _candidate_id(candidate))),
                       node, "longest", "longest_corrected_cds") for node in nodes],
                "scores": [], "audit": [], "omitted_edges": [],
                "metrics": {"loci": len(nodes), "pair_evaluations": 0}}

    prepared = prepare_correspondence(
        nodes, edges, coordinate_ambiguity={node for node in nodes if loci[node].get("ambiguous_coordinates", False)},
        nonvoting={node for node in nodes if not options[node]},
        copy_ambiguity=_copy_ambiguity)
    coordinate_ambiguity, ambiguity = prepared["coordinate_ambiguity"], prepared["ambiguity"]
    omitted_edges, trusted_edges = prepared["omitted_edges"], prepared["edge_by_pair"]
    adjacency = defaultdict(list)
    for (a, b), weight in sorted(trusted_edges.items()):
        adjacency[a].append((b, weight))
        adjacency[b].append((a, weight))
    pairs = _PairScores(max_alignment_cells, cache_pair_scores, pair_cache_limit)
    graph_state = None
    voting, totals, normalized, incident = {}, {}, {}, {}

    def voting_graph(state: dict) -> None:
        nonlocal graph_state, voting, totals, normalized, incident
        if state is graph_state:
            return
        # Coordinate-optimization states contain only eligible options; their
        # eligibility cannot change during an update.  Adoption states are
        # immutable snapshots.  Cache one graph view per snapshot so scoring
        # each candidate remains proportional to its degree.
        graph_state = state
        # Trial states contain one complete connected component; adoption
        # snapshots contain every locus.  Neither scoring path needs to scan
        # loci or edges outside the supplied state's scope.
        scoped_nodes = sorted(state)
        voting = {node: is_eligible(state[node]) for node in scoped_nodes}
        totals = {node: sum(weight for neighbor, weight in adjacency[node] if voting[neighbor]) for node in scoped_nodes}
        normalized = {(a, b): weight * 0.5 * (1 / totals[a] + 1 / totals[b])
                      for a in scoped_nodes for b, weight in adjacency[a] if a < b and voting[a] and voting[b]}
        incident = dict.fromkeys(scoped_nodes, 0.0)
        for (a, b), weight in sorted(normalized.items()):
            incident[a] += weight
            incident[b] += weight

    def local_score(node: tuple, candidate: dict, state: dict) -> float:
        score = quality[node][_candidate_id(candidate)]
        if not is_eligible(candidate):
            return score
        voting_graph(state)
        weights = []
        for neighbor, weight in adjacency[node]:
            if not voting[neighbor]:
                continue
            if voting[node]:
                coefficient = normalized[tuple(sorted((node, neighbor)))]
            else:
                # Score an eligible alternative to an ineligible baseline on
                # the graph that would result from adopting that alternative.
                coefficient = weight * 0.5 * (1 / totals[node] + 1 / (totals[neighbor] + weight))
            weights.append((neighbor, coefficient))
        denominator = sum(weight for _, weight in weights)
        for neighbor, weight in weights:
            score += weight / denominator * pairs.get(candidate, state[neighbor]).score
        return score

    def objective(component: list, state: dict) -> float:
        # F_g = Q_g + sum_h(e_gh / I_g) C_gh is comparable across graph
        # degrees, where I_g = sum_h e_gh.  The symmetric global objective
        # J = sum_g I_g Q_g + sum_edges e_gh C_gh satisfies delta J = I_g
        # delta F_g for every coordinate update.  Isolated loci use Q alone.
        voting_graph(state)
        value = sum((incident[node] or 1.0) * quality[node][_candidate_id(state[node])]
                    for node in component if voting[node])
        for a in component:
            for b, _ in adjacency[a]:
                if a < b and voting[a] and voting[b]:
                    value += normalized[(a, b)] * pairs.get(state[a], state[b]).score
        return value

    state = dict(baselines)
    components, visited = [], set()
    for node in nodes:
        if node in visited:
            continue
        stack, component = [node], []
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            component.append(current)
            stack.extend(neighbor for neighbor, _ in adjacency[current] if neighbor not in visited)
        components.append(sorted(component))
    optimization_audit = []
    optimized_nodes = set()
    for component in components:
        needs_comparison = any(len(options[node]) > 1 or
                               (options[node] and baselines[node] not in options[node]) for node in component)
        if not needs_comparison or not any(adjacency[node] for node in component):
            continue
        optimized_nodes.update(component)
        active = {node: options[node] or [baselines[node]] for node in component}
        independent = {node: min(active[node], key=lambda candidate: (-quality[node][_candidate_id(candidate)],
                                                                    _candidate_id(candidate))) for node in component}
        seeds = [{node: baselines[node] if baselines[node] in active[node] else independent[node] for node in component},
                 independent]
        for rank in range(min(3, max(map(len, active.values())))):
            seeds.append({node: active[node][min(rank, len(active[node]) - 1)] for node in component})
        best, best_objective, best_signature = None, -math.inf, ()
        unique_seeds = set()
        rounds = 0
        for seed in seeds:
            signature = tuple(_candidate_id(seed[node]) for node in component)
            if signature in unique_seeds:
                continue
            unique_seeds.add(signature)
            trial = dict(seed)
            for iteration in range(20):
                changed = False
                for node in component:
                    current = trial[node]
                    evaluated = [(local_score(node, candidate, trial), candidate) for candidate in active[node]]
                    selected_score, selected = min(evaluated, key=lambda item: (-item[0], item[1] is not current,
                                                                               _candidate_id(item[1])))
                    current_score = next(score for score, candidate in evaluated if candidate is current)
                    if selected_score > current_score + 1e-12:
                        trial[node] = selected
                        changed = True
                rounds = max(rounds, iteration + 1)
                if not changed:
                    break
            value = objective(component, trial)
            signature = tuple(_candidate_id(trial[node]) for node in component)
            # Exact objective ties prefer retaining more baseline records, then
            # stable IDs, independently of catalog and edge order.
            retention = sum(trial[node] == baselines[node] for node in component)
            tie_signature = (-retention, signature)
            if value > best_objective + 1e-12 or (abs(value - best_objective) <= 1e-12 and tie_signature < best_signature):
                best, best_objective, best_signature = trial, value, tie_signature
        count = math.prod(len(active[node]) for node in component)
        exact_value = None
        if exact_limit and count <= exact_limit:
            exact_value = -math.inf
            for values in itertools.product(*(active[node] for node in component)):
                exact_value = max(exact_value, objective(component, dict(zip(component, values, strict=True))))
        state.update({node: best[node] for node in component})
        component_audit = {"loci": [list(node) for node in component], "state_count": count,
                           "seeds": len(unique_seeds), "rounds": rounds, "objective": best_objective,
                           "exact_objective": exact_value,
                           "heuristic_matches_exact": None if exact_value is None else
                           abs(exact_value - best_objective) <= 1e-10}
        if count.bit_length() > 12_000:
            # A malformed giant family can exceed Python's safe JSON integer
            # conversion limit even though the heuristic itself stays bounded.
            component_audit["state_count"] = None
            component_audit["log10_state_count"] = sum(math.log10(len(active[node])) for node in component)
        optimization_audit.append(component_audit)

    # Adoption is synchronous and is then rechecked after abstaining neighbors
    # return to their originals.  A proposal cannot depend on unadopted donors.
    decisions, score_rows = {}, []

    def comparison_scores(node: tuple, candidate: dict, frozen: dict) -> list[float]:
        alternatives = [other for other in options[node] if other != candidate]
        baseline = baselines[node]
        if baseline != candidate and baseline not in options[node]:
            # A sole usable repair is still compared with its actual source
            # baseline, rather than inventing a margin for a one-element rank.
            alternatives.append(baseline)
        return [local_score(node, other, frozen) for other in alternatives]

    def decision(node: tuple, frozen: dict) -> dict:
        baseline, proposed = baselines[node], frozen[node]
        if node in ambiguity:
            reason = "ambiguous_locus_coordinates" if node in coordinate_ambiguity else "unresolved_paralog_copy"
            return _selection(baseline, node, "ambiguous_correspondence", reason)
        if len(candidates[node]) > candidate_limit:
            return _selection(baseline, node, "insufficient_evidence", "candidate_limit_exceeded")
        if not options[node]:
            return _selection(baseline, node, "insufficient_evidence", "no_eligible_candidate")
        if not adjacency[node]:
            return _selection(baseline, node, "insufficient_evidence", "no_trusted_crossspecies_correspondence")
        if node not in optimized_nodes:
            return _selection(baseline, node, "unchanged", "single_eligible_candidate")
        ranking = sorted(((local_score(node, candidate, frozen), candidate) for candidate in options[node]),
                         key=lambda item: (-item[0], item[1] != baseline, _candidate_id(item[1])))
        best_score, best_candidate = ranking[0]
        proposed_score = local_score(node, proposed, frozen)
        alternatives = comparison_scores(node, proposed, frozen)
        margin = proposed_score - max(alternatives) if alternatives else None

        def abstain(reason: str) -> dict:
            baseline_score = local_score(node, baseline, frozen)
            baseline_alternatives = comparison_scores(node, baseline, frozen)
            baseline_margin = baseline_score - max(baseline_alternatives) if baseline_alternatives else None
            result = _selection(baseline, node, "insufficient_evidence", reason, baseline_score, baseline_margin)
            result.update(proposed_candidate_id=_candidate_id(proposed), proposed_score=proposed_score,
                          proposed_margin=margin)
            return result

        if proposed == baseline:
            return _selection(baseline, node, "unchanged", "source_representative_retained", proposed_score, margin)
        if proposed != best_candidate or margin is None or margin < min_margin:
            return abstain("candidate_score_margin")
        support = set()
        for neighbor, _ in adjacency[node]:
            donor = frozen[neighbor]
            if not is_eligible(donor):
                continue
            pair = pairs.get(proposed, donor)
            if not pair.bounded and pair.identity >= min_identity and min(
                    pair.coverage_a, pair.coverage_b) >= min_coverage:
                support.add(neighbor[0])
        if len(support) < min_support:
            return abstain("insufficient_independent_donor_species")
        # Avoid making every species' annotation match a shared short fragment.
        if _length(proposed) < 0.70 * _length(baseline) and not _independent_support(proposed):
            return abstain("shorter_than_source_major_coding_span")
        if is_eligible(baseline) and not _independent_support(proposed):
            baseline_pair = pairs.get(baseline, proposed)
            if baseline_pair.internal_loss_a:
                return abstain("possible_species_specific_internal_region")
        return _selection(proposed, node, "conserved", "supported_crossspecies_coding_isoform", proposed_score, margin)

    adopted_state, rejected = dict(state), set()
    for _ in range(len(nodes) + 1):
        for node in nodes:
            if node not in rejected:
                decisions[node] = decision(node, adopted_state)
                if decisions[node]["candidate_id"] != _candidate_id(adopted_state[node]):
                    rejected.add(node)
        next_state = {node: next(candidate for candidate in candidates[node]
                                 if _candidate_id(candidate) == decisions[node]["candidate_id"]) for node in nodes}
        if next_state == adopted_state:
            break
        adopted_state = next_state
    else:
        raise RuntimeError("Conservative adoption failed to converge")
    for node in sorted(optimized_nodes):
        selected_score = local_score(node, adopted_state[node], adopted_state)
        alternatives = comparison_scores(node, adopted_state[node], adopted_state)
        decisions[node]["score"] = selected_score
        decisions[node]["margin"] = selected_score - max(alternatives) if alternatives else None
        if "proposed_candidate_id" in decisions[node]:
            proposal = next(candidate for candidate in candidates[node]
                            if _candidate_id(candidate) == decisions[node]["proposed_candidate_id"])
            proposal_score = local_score(node, proposal, adopted_state)
            alternatives = comparison_scores(node, proposal, adopted_state)
            decisions[node].update(proposed_score=proposal_score,
                                   proposed_margin=proposal_score - max(alternatives) if alternatives else None)
        for candidate in candidates[node]:
            score_rows.append({"species": node[0], "gene_id": node[1], "candidate_id": _candidate_id(candidate),
                               "eligible": is_eligible(candidate), "quality_score": quality[node][_candidate_id(candidate)],
                               "score": local_score(node, candidate, adopted_state) if is_eligible(candidate) else None,
                               "selected": candidate == adopted_state[node]})
    return {"schema": 1, "policy": policy, "selections": [decisions[node] for node in nodes],
            "scores": score_rows, "audit": optimization_audit, "omitted_edges": omitted_edges,
            "parameters": {"min_margin": min_margin, "min_support": min_support, "min_identity": min_identity,
                           "min_coverage": min_coverage, "candidate_limit": candidate_limit,
                           "max_alignment_cells": max_alignment_cells, "short_span_fraction": 0.70,
                           "cutoffs_calibrated": False},
            "metrics": {"loci": len(nodes), "trusted_edges": len(trusted_edges), "ambiguous_loci": len(ambiguity),
                        "components": len(components), "optimized_loci": len(optimized_nodes),
                        "pair_requests": pairs.requests, "pair_evaluations": pairs.evaluations,
                        "pair_cache_hits": pairs.hits, "bounded_alignments": pairs.bounded,
                        "pair_cache_limit": pair_cache_limit, "pair_cache_peak_entries": pairs.peak_entries,
                        "pair_cache_evictions": pairs.evictions,
                        "changed_representatives": sum(adopted_state[node] != baselines[node] for node in nodes)}}
