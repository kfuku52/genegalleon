# GIFT acquisition integrity audit

The September 2026 audit follows the 40-species pilot and covers scientific-name
resolution, trait identity, pagination, cache refresh/resume, source parsing,
aggregation and publication through the input-generation workflow.

## Reproduced defects and corrections

The initial adversarial suite reproduced failures in 24 of 25 cases. These were
multiple cases of the defect classes below, not 24 independent root causes.

| Defect | Correction | Regression evidence |
| --- | --- | --- |
| Interrupted refresh could combine a fresh first page with an old tail; concurrent readers could change generations mid-run | Query-scoped generations, pinned per reader; refresh starts a new namespace before requesting data | Interrupted refresh/resume and interleaved reader/writer tests |
| A damaged cache envelope could raise `AttributeError` instead of recovering | Validate envelope type, schema, exact URL, timestamp and payload checksum before reuse | List/null/string/envelope corruption cases |
| JSON error rows with a `work_ID`, invalid IDs or invalid agreement scores could be accepted | Validate required value fields and scalar types, positive identifiers and finite 0–1 agreement values | Malformed API-row cases; no cached data written |
| HTTP-date `Retry-After` was ignored; truncated HTTP bodies were not retried | Honor both header forms and retry transient incomplete reads within bounds | Future HTTP-date and truncated-body tests |
| An ambiguous trait group could select its most common member | Require an exact name resolving to one ID; ambiguity is an error | Min/max height group test |
| ID/name aliases could duplicate one physical trait's observations | Fetch each ID once and retain an alias-to-ID map for plan filtering | Helper and full CLI binary-sum tests |
| Explicitly unresolved taxonomic hits and qualified targets could be treated as exact species matches | Exclude unresolved flags and report qualified taxonomic scope separately | Unresolved-name and infraspecific-target tests |
| Configuration typos, duplicate headers or extra fields could silently change interpretation | Validate declarations, reserved names and duplicate database definitions | Malformed plan and mapping-header cases |
| Malformed/infinite numeric values could disappear or enter the output; six-digit formatting lost precision | Reject invalid numeric observations and overflowing aggregates; use round-trip float formatting | Numeric invalid/overflow/precision cases |
| Free text could be coerced to numeric or reduced to one category | Explicit `text` type with a sorted JSON array of distinct values; retain literal missing words in text input | Text aggregation and CSV-reader cases |
| Output/stats paths could alias an input or each other; publication could leave a partial result on ordinary errors | Preflight path aliases and file types, prepare all payloads, install atomically per file and roll back ordinary failures | Input/output alias and injected second-installation failure tests |
| Strict dry-run expected data that it deliberately did not retrieve | Validate configuration but skip observation-dependent checks during dry-run | Full CLI dry-run preserves the existing result |

The existing GBIF plan producer now declares its binary aggregation explicitly
instead of depending on an undocumented fallback from `median` to `any`.
Non-strict output tables retain requested all-missing columns so downstream
observation counts distinguish unavailable data from absent schema.

## Dataset validation

All 109 public GIFT 3.2 catalog items were replayed through the production CLI
using checksum-verified saved API responses. The resulting 40-species table
matched the reviewed reference: 81 observed traits and 1,096 observed cells.
Numeric values were compared numerically and text arrays were compared with the
previous complete sets of strings. A following offline run produced a
byte-identical TSV and made no fixture or remote requests. This replay validates
processing of the full saved dataset; it is not a fresh 109-trait network download.

The live check separately queries woodiness for all 40 species using the public
API and compares a subsequent offline run with that live result and the reviewed
reference. The reviewed synonym mappings recover four previously missed species;
the broader Citrus aggregate remains excluded.

## Boundaries

The API still requires global pages on a cold acquisition. GIFT releases and
reviewed taxonomic decisions are explicit in retrieval reports; no atomic remote
snapshot is implied for a changing server, particularly beta. Missing public
observations are not imputed. File rollback covers ordinary publication errors,
not a hard process kill between separate file replacements. The workflow records
a successful output receipt only after completion.

Automated regression suites and workflow integration run in the GeneGalleon
Docker runtime. These checks do not establish SIF runtime compatibility.
See [retrieval controls and schemas](gift-trait-retrieval.md) and the executable
[test suite](../workflow/tests/test_gift_retrieval_audit.py).
