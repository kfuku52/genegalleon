"""Filter modeled species branches using existing host-species taxonomy."""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

from hgt_species_tree import event_species_branch

DIRECTION_CHOICES = ('off', 'any', 'custom', 'non_arthropoda_to_insecta', 'non_arthropoda_to_arthropoda')
DISABLED = ('off', 'any')
EVENT_FIELDS = ['direction_filter_status', 'direction_filter_reason',
                'direction_donor_arthropoda_status', 'direction_recipient_insecta_status',
                'direction_recipient_arthropoda_status', 'direction_donor_status',
                'direction_recipient_status', 'direction_rules_sha256']
MISSING = {'', '.', 'na', 'nan', 'none', 'null', 'unknown', 'unavailable'}


def key(value):
    return str(value).strip().replace(' ', '_')


def known(value):
    return str(value).strip().lower() not in MISSING


def read_taxonomy(path, ranks=('domain', 'phylum', 'class')):
    raw = Path(path).read_bytes()
    import io
    reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig')), delimiter='\t')
    fields = reader.fieldnames or []
    required = {'species', 'resolution_status', *ranks}
    if len(fields) != len(set(fields)) or not required <= set(fields):
        raise ValueError('Missing or duplicate species-taxonomy columns')
    aliases, species = {}, set()
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise ValueError('Malformed species-taxonomy row')
        name = key(row['species'])
        if not name or name in species:
            raise ValueError('Duplicate or empty taxonomy species identity')
        species.add(name)
        if row['resolution_status'] == 'resolved' and row.get('class') == 'Insecta' and row.get('phylum', 'Arthropoda') != 'Arthropoda':
            raise ValueError('Insecta classification requires Arthropoda phylum: ' + name)
        names = {name}
        if row.get('tree_status') == 'mapped' and known(row.get('tree_tip', '')):
            names.add(key(row['tree_tip']))
        for alias in names:
            if alias in aliases:
                raise ValueError('Ambiguous taxonomy species/tree-tip alias: ' + alias)
            aliases[alias] = row
    return aliases, hashlib.sha256(raw).hexdigest()


def tip_state(row, group):
    if not row or row['resolution_status'] != 'resolved':
        return 'unknown'
    phylum, domain, klass = (row[name].strip() for name in ('phylum', 'domain', 'class'))
    if group == 'arthropoda':
        return 'within' if phylum == 'Arthropoda' else 'outside' if known(phylum) or domain in {'Bacteria', 'Archaea'} else 'unknown'
    if klass == 'Insecta':
        return 'within'
    return 'outside' if known(klass) or domain in {'Bacteria', 'Archaea'} or known(phylum) and phylum != 'Arthropoda' else 'unknown'


def classify_branches(nodes, taxonomy):
    rows, states = [], {}
    for name, tips in sorted(nodes.items()):
        result = dict(species_branch=name, clade_tip_labels='; '.join(tips), descendant_tip_count=len(tips))
        for group in ('arthropoda', 'insecta'):
            counts = Counter(tip_state(taxonomy.get(tip), group) for tip in tips)
            status = ('mixed' if counts['within'] and counts['outside'] else 'unknown' if counts['unknown'] or not tips
                      else 'within' if counts['within'] else 'outside')
            result[group + '_status'] = status
            for state in ('within', 'outside', 'unknown'):
                result[f'{group}_{state}_tip_count'] = counts[state]
            result[group + '_unknown_tip_labels'] = '; '.join(tip for tip in tips if tip_state(taxonomy.get(tip), group) == 'unknown')
        states[name] = result
        rows.append(result)
    return states, rows


def _filter_legacy_events(events, nodes, taxonomy_path, direction_filter):
    """Evaluate each modeled branch's entire descendant clade, preserving unknowns."""
    if direction_filter not in DIRECTION_CHOICES:
        raise ValueError('Invalid focused HGT direction filter')
    recipient_group = 'arthropoda' if direction_filter == 'non_arthropoda_to_arthropoda' else 'insecta'
    taxonomy, source_sha = read_taxonomy(taxonomy_path)
    states, branches = classify_branches(nodes, taxonomy)
    selected, audit = [], []
    for event in events:
        donor, recipient = (event_species_branch(event, side) for side in ('donor', 'recipient'))
        if event['generax_transfer'] != f"Y@{event['generax_donor_node']}@{event['generax_recipient_node']}":
            raise ValueError('Event transfer token disagrees with donor/recipient branch IDs')
        ds = states.get(donor, {}).get('arthropoda_status', 'unmapped')
        rs = states.get(recipient, {}).get(recipient_group + '_status', 'unmapped')
        reason = ('' if direction_filter == 'any' else
                  'event_mapping_unresolved' if event.get('mapping_status', 'matched') != 'matched'
                  else 'donor_' + ds + '_arthropoda' if ds != 'outside'
                  else 'recipient_' + rs + '_' + recipient_group if rs != 'within' else '')
        status = ('passed' if not reason else 'withheld' if event.get('mapping_status', 'matched') != 'matched'
                  or ds in {'unknown', 'mixed', 'unmapped'} or rs in {'unknown', 'mixed', 'unmapped'} else 'excluded_direction')
        row = dict(event, direction_filter_status=status,
                   direction_filter_reason=reason or ('direction_filter_disabled' if direction_filter == 'any' else
                       'all_donor_tips_outside_arthropoda_all_recipient_tips_in_' + recipient_group),
                   direction_donor_arthropoda_status=ds,
                   direction_recipient_insecta_status=states.get(recipient, {}).get('insecta_status', 'unmapped'),
                   direction_recipient_arthropoda_status=states.get(recipient, {}).get('arthropoda_status', 'unmapped'))
        audit.append(row)
        if status == 'passed':
            selected.append(row)
    return selected, audit, branches, source_sha


def _unique_object(pairs):
    result = {}
    for name, value in pairs:
        if name in result:
            raise ValueError('Duplicate direction rule key: ' + name)
        result[name] = value
    return result


def read_rules(path):
    """Read an explicit data-only donor/recipient selection contract."""
    raw = Path(path).read_bytes()
    rules = json.loads(raw.decode('utf-8-sig'), object_pairs_hook=_unique_object)
    if not isinstance(rules, dict) or set(rules) != {'schema_version', 'donor', 'recipient'}:
        raise ValueError('Direction rules require schema_version, donor and recipient only')
    if type(rules['schema_version']) is not int or rules['schema_version'] != 1:
        raise ValueError('Unsupported direction rules schema_version')
    for side in ('donor', 'recipient'):
        rule = rules[side]
        if not isinstance(rule, dict) or rule.get('kind') not in {'any', 'taxonomy', 'tips', 'branches'}:
            raise ValueError('Invalid direction selector kind: ' + side)
        if rule['kind'] == 'any':
            if set(rule) != {'kind'}:
                raise ValueError('An any selector cannot contain constraints')
            continue
        allowed = {'kind', 'values', 'match'} | ({'rank'} if rule['kind'] == 'taxonomy' else set())
        if set(rule) - allowed or not {'kind', 'values'} <= set(rule):
            raise ValueError('Invalid direction selector fields: ' + side)
        if rule.get('match', 'within') not in {'within', 'outside'}:
            raise ValueError('Direction match must be within or outside')
        values = rule['values']
        if not isinstance(values, list) or not values or any(not isinstance(v, str) or not known(v) for v in values):
            raise ValueError('Direction values must be nonempty strings')
        values = [v.strip() if rule['kind'] == 'taxonomy' else key(v) for v in values]
        if len(values) != len(set(values)):
            raise ValueError('Duplicate direction selector values')
        rule['values'] = values
        if rule['kind'] == 'taxonomy':
            rank = rule.get('rank')
            if (not isinstance(rank, str) or not rank.strip() or rank != rank.strip()
                    or rank in {'species', 'resolution_status', 'tree_status', 'tree_tip', 'taxid'}
                    or rank.endswith('_taxid')):
                raise ValueError('Direction taxonomy rank must name a taxonomy rank column (species_rank for species)')
    return rules, hashlib.sha256(raw).hexdigest()


def validate_options(direction_filter, direction_rules='', taxonomy_path=''):
    if direction_filter not in DIRECTION_CHOICES:
        raise ValueError('Invalid focused HGT direction filter')
    if direction_filter == 'custom':
        if not direction_rules:
            raise ValueError('Custom direction filtering requires direction_rules')
        rules, _ = read_rules(direction_rules)
        if any(rules[side]['kind'] == 'taxonomy' for side in ('donor', 'recipient')) and not taxonomy_path:
            raise ValueError('Focused direction filtering requires existing species taxonomy')
    elif direction_rules:
        raise ValueError('direction_rules requires the explicit custom direction filter')
    elif direction_filter not in DISABLED and not taxonomy_path:
        raise ValueError('Focused direction filtering requires existing species taxonomy')


def filter_label(direction_filter):
    if direction_filter == 'custom':
        return 'Configured donor → recipient'
    if direction_filter in DISABLED:
        return 'Direction filter off'
    group = 'Arthropoda' if direction_filter == 'non_arthropoda_to_arthropoda' else 'Insecta'
    return f'Non-Arthropoda donor → {group} recipient'


def _rank_state(row, rule):
    if not row or row['resolution_status'] != 'resolved':
        return 'unknown'
    value = row[rule['rank']].strip()
    if not known(value):
        return 'unknown'
    values = json.loads(value) if value.startswith('[') else [value]
    if not isinstance(values, list) or not values or any(not isinstance(v, str) or not known(v) for v in values):
        raise ValueError('Malformed multi-valued taxonomy rank')
    return 'within' if set(values) & set(rule['values']) else 'outside'


def _custom_branch(name, tips, rule, taxonomy):
    kind = rule['kind']
    if kind == 'branches':
        status = 'within' if name in rule['values'] else 'outside'
        return status, {}
    counts = Counter('within' if kind == 'any' or kind == 'tips' and tip in rule['values'] else
                     'outside' if kind == 'tips' else _rank_state(taxonomy.get(tip), rule) for tip in tips)
    status = ('mixed' if counts['within'] and counts['outside'] else 'unknown' if counts['unknown'] or not tips
              else 'within' if counts['within'] else 'outside')
    return status, counts


def filter_events(events, nodes, taxonomy_path='', direction_filter='off', direction_rules=''):
    """Opt-in endpoint rules; legacy names remain explicit compatibility presets."""
    validate_options(direction_filter, direction_rules, taxonomy_path)
    if direction_filter not in DISABLED and direction_filter != 'custom':
        return _filter_legacy_events(events, nodes, taxonomy_path, direction_filter)
    if direction_filter in DISABLED:
        audit = []
        for event in events:
            if event['generax_transfer'] != f"Y@{event['generax_donor_node']}@{event['generax_recipient_node']}":
                raise ValueError('Event transfer token disagrees with donor/recipient branch IDs')
            audit.append(dict(event, direction_filter_status='passed', direction_filter_reason='direction_filter_disabled'))
        return audit, audit, [], ''
    rules, rules_sha = read_rules(direction_rules)
    tips = {tip for descendants in nodes.values() for tip in descendants}
    for side in ('donor', 'recipient'):
        rule = rules[side]
        allowed = tips if rule['kind'] == 'tips' else set(nodes)
        if rule['kind'] in {'tips', 'branches'} and set(rule['values']) - allowed:
            raise ValueError('Direction selector identities absent from analysis species tree: ' + side)
    ranks = {rules[s]['rank'] for s in ('donor', 'recipient') if rules[s]['kind'] == 'taxonomy'}
    taxonomy, taxonomy_sha = read_taxonomy(taxonomy_path, ranks) if ranks else ({}, '')
    states, branches = {}, []
    for name, descendants in sorted(nodes.items()):
        row = dict(species_branch=name, clade_tip_labels='; '.join(descendants), descendant_tip_count=len(descendants))
        for side in ('donor', 'recipient'):
            status, counts = _custom_branch(name, descendants, rules[side], taxonomy)
            if not descendants:
                status = 'unknown'
            row[side + '_status'] = status
            for state in ('within', 'outside', 'unknown'):
                row[f'{side}_{state}_tip_count'] = counts.get(state, 0) if rules[side]['kind'] != 'branches' else ''
        branches.append(row)
        states[name] = row
    selected, audit = [], []
    for event in events:
        if event['generax_transfer'] != f"Y@{event['generax_donor_node']}@{event['generax_recipient_node']}":
            raise ValueError('Event transfer token disagrees with donor/recipient branch IDs')
        actual = {s: states.get(event_species_branch(event, s), {}).get(s + '_status', 'unmapped')
                  for s in ('donor', 'recipient')}
        reason = 'event_mapping_unresolved' if event.get('mapping_status', 'matched') != 'matched' else ''
        if not reason:
            reason = next((s + '_' + actual[s] for s in ('donor', 'recipient')
                           if actual[s] != rules[s].get('match', 'within')), '')
        status = ('passed' if not reason else 'withheld' if event.get('mapping_status', 'matched') != 'matched'
                  or any(v in {'mixed', 'unknown', 'unmapped'} for v in actual.values()) else 'excluded_direction')
        row = dict(event, direction_filter_status=status, direction_filter_reason=reason or 'configured_endpoint_rules_matched',
                   direction_donor_status=actual['donor'], direction_recipient_status=actual['recipient'],
                   direction_rules_sha256=rules_sha)
        audit.append(row)
        if status == 'passed':
            selected.append(row)
    return selected, audit, branches, taxonomy_sha
