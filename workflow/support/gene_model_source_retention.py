"""Keep known source coding structure when homology repairs an incomplete locus.

A damaged whole ORF does not invalidate every retained exon or its complete
termini. Changes to those coordinates need target whole-path RNA support;
adding missing exons or extending an incomplete terminus remains possible.
"""
from Bio.Data import CodonTable


def source_structure_retention(candidate, originals, genetic_code):
    if not candidate.get('blocks'):
        return {'preserved': False, 'sources': []}
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    checks = []
    for original in originals:
        if (original.get('seqid', candidate['seqid']) != candidate['seqid']
                or original.get('strand', candidate['strand']) != candidate['strand']
                or not original.get('blocks')):
            continue
        blocks = sorted(original['blocks'], reverse=candidate['strand'] == '-')
        predicted = sorted(candidate['blocks'], reverse=candidate['strand'] == '-')
        dna = original.get('cds', '').upper()
        explicit_partial = original.get('quality', {}).get('annotated_partial', False)
        complete_start = not explicit_partial and blocks[0][2] == 0 and dna[:3] in table.start_codons
        complete_stop = not explicit_partial and dna[-3:] in table.stop_codons
        changed = []
        for i, (a, b, _phase) in enumerate(blocks):
            matches = []
            for j, (x, y, _predicted_phase) in enumerate(predicted):
                left_ok = x == a
                right_ok = y == b
                if candidate['strand'] == '+':
                    left_ok |= i == j == 0 and not complete_start and x <= a
                    right_ok |= i == len(blocks) - 1 and j == len(predicted) - 1 and not complete_stop and y >= b
                else:
                    right_ok |= i == j == 0 and not complete_start and y >= b
                    left_ok |= i == len(blocks) - 1 and j == len(predicted) - 1 and not complete_stop and x <= a
                if left_ok and right_ok:
                    matches.append(j)
            if not matches:
                changed.append([a, b])
        five_prime = (lambda block: block[0]) if candidate['strand'] == '+' else (lambda block: block[1])
        three_prime = (lambda block: block[1]) if candidate['strand'] == '+' else (lambda block: block[0])
        start_preserved = not complete_start or five_prime(blocks[0]) == five_prime(predicted[0])
        stop_preserved = not complete_stop or three_prime(blocks[-1]) == three_prime(predicted[-1])
        checks.append({'source_candidate_id': original['candidate_id'], 'changed_source_blocks': changed,
                       'complete_start_preserved': start_preserved, 'complete_stop_preserved': stop_preserved,
                       'preserved': not changed and start_preserved and stop_preserved})
    return {'preserved': any(row['preserved'] for row in checks), 'sources': checks}
