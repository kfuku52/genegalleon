# Donor-supported terminal completion

Missing start/stop codons can trigger a bounded search of the contiguous
genomic flank of the first/last coding exon. The search preserves the existing
coding sequence, reading frame, exon junctions and phases. It uses the species'
genetic code and never manufactures a codon, masks a stop or crosses an
in-frame stop or an ambiguous assembly base. The existing validator separately
handles an exact next genomic stop codon.

Defaults are 300 nucleotides per terminus, at most 32 start candidates and
25 million protein alignment cells per candidate. The donor query's previously
unaligned terminal residues must support the extension in a new global protein
alignment. Added amino acids must map to that terminal segment, both donor
termini must align, and a missing internal donor segment cannot be called
terminal recovery. Identity and reciprocal coverage must meet the existing
species-specific rescue thresholds. Permitted alternative initiation codons
use methionine for the protein alignment, without changing genomic DNA.

When the original alignment already reaches the donor's last residue, a target
species can have a short C-terminal overhang before its genomic stop. This is
limited to **two amino acids**, requires that the unchanged CDS's last residue
still align to the donor's last residue, and retains the same reciprocal
coverage and identity gates. `terminal_max_unaligned_c_overhang=0` disables
that exception; values above two are refused. Long unsupported tails stay
partial evidence.

Competing eligible starts with a score difference at or below 0.02 stay
ambiguous. The ranking score is identity multiplied by the lower reciprocal
coverage; it is a ranking rule, not a calibrated probability of correctness.
Up to eight optimal protein alignments are examined. Different terminal
projections or unexamined additional optimal paths remain ambiguous evidence;
the first traceback is never treated as uniquely supporting a coding start.
No extension can remove an existing frameshift, internal stop, invalid phase,
splice defect, assembly ambiguity or sequence mismatch. After independent
donor alignment, the original strict genomic validator checks the completed
chain again. Synteny, genomic ownership and copy correspondence remain the
caller's additional admission gates.

The `terminal_completion` audit preserves the original coding chain, source
DNA hash and miniprot alignment, each attempted chain, the exact added genomic
DNA, recalculated protein evidence and rejection reasons. A successful model
identifies its metrics as `alignment_source=terminal_completion_global`;
the original PAF is retained as provenance. Failed or ambiguous models retain
their original sequence and problems, with explicit `partial_evidence` that
does not grant intact representative eligibility. Partial evidence is not
proof of a pseudogene or true gene loss.
