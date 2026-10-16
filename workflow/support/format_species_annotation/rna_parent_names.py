"""Recover an RNA ID only from its explicit Name and exact child references."""
import hashlib
import os
import tempfile
from collections import Counter
from pathlib import Path

from format_species_writers import open_text


def repair_rna_parent_names(source, scratch, candidates, anchors, ids, protected):
    """A self-parented gene/RNA export can name its actual RNA in Name."""
    from .locus_identity import contains, geometry, rewrite_attributes, rows
    selected = {}
    names = Counter(attrs["Name"][0] for _number, _fields, attrs in candidates)
    for number, fields, attrs in candidates:
        identifier, name = attrs['ID'][0], attrs['Name'][0]
        models = anchors[identifier]
        genes = [scope for scope in models if scope[0] == 'gene']
        rnas = [scope for scope in models if scope[0] == 'mrna']
        if (identifier in protected or name in protected or len(genes) != 1 or len(rnas) != 1
                or len(models) != 2 or not contains(genes[0][1], rnas[0][1])
                or name in ids or names[name] != 1):
            continue
        selected[name] = (number, fields, attrs)
    if not selected:
        return Path(source), []
    children, blocked = set(), set()
    for _number, fields, attrs in rows(source):
        names = set(attrs.get('Parent', ())) & selected.keys()
        for name in names:
            parent = selected[name][1]
            if (len(attrs.get('Parent', ())) != 1 or not contains(geometry(parent), geometry(fields))
                    or attrs.get('exception') or attrs.get('part') or attrs.get('is_ordered') == ('true',)):
                blocked.add(name)
            elif fields[2].lower() == 'cds':
                children.add(name)
    selected = {name: item for name, item in selected.items() if name in children and name not in blocked}
    if not selected:
        return Path(source), []
    renamed = {number: (name, fields, attrs) for name, (number, fields, attrs) in selected.items()}
    with Path(source).open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    output = Path(scratch) / (digest + '.rna-parent-names.gff')
    mappings = []
    with open_text(source, 'rt') as src, tempfile.NamedTemporaryFile(mode='w', dir=scratch, delete=False) as dest:
        dest.write('##genegalleon-rna-parent-names version=1 source_sha256=' + digest + '\n')
        for number, line in enumerate(src, 1):
            if number in renamed:
                name, fields, attrs = renamed[number]
                original = attrs['ID'][0]
                # Replace only ID. The gene Parent and every coding coordinate
                # already refer to the author-declared model and stay intact.
                changed = fields[:]
                fields_id_only = ';'.join(item for item in fields[8].split(';') if item.partition('=')[0].strip() == 'ID')
                replacement = rewrite_attributes(fields_id_only, {original: name})
                changed[8] = ';'.join(replacement if item.partition('=')[0].strip() == 'ID' else item
                                      for item in fields[8].split(';'))
                line = '\t'.join(changed) + '\n'
                mappings.append(dict(source_id=original, target_id=name, feature_type='mrna',
                    seqid=fields[0], strand=fields[6], start=int(fields[3]), end=int(fields[4]),
                    line=number, reason='explicit_rna_name_parent'))
            dest.write(line)
    os.replace(dest.name, output)
    return output, mappings
