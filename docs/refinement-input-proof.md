# Refinement input verification within a command

Refinement retains full-byte-first input proofs only inside its explicit CLI
invocation. The first use of each input scope reads and checks its frozen SHA-256.
Later stage guards retain the ordinary names selection: global request files and
the requested species' FASTA, GFF and genome. Library calls outside that command
scope continue to perform fresh full-byte verification.

Small plan JSON is read every time. Its whole value, original generation and
pathname must match the first frozen copy. Implementation, runtime dependency
and predictor identities are still checked ordinarily on each load. Input and
plan changes are sticky failures, including same-byte atomic replacements;
restoring an observed generation does not refresh a proof within that command.

The existing preparation proof fences exact resolved targets, file identities,
permissions, symlinks and ancestor identities. Ordinary ancestors can gain
unrelated output siblings. Shared physical aliases reuse a full-byte hash only
after its original proof and each new alias pathname have been fenced.

Input proofs are checked again after output hashing and before publication, and
all used scopes are checked before successful command exit. No proof persists
to a later CLI, and dependency receipt/content verification remains independent.
This changes neither scientific admission nor output keys except the normal
implementation fingerprint. Existing frozen plans require their original
implementation; validation of a changed implementation uses a fresh plan.
