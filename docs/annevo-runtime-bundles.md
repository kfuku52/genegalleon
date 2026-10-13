# ANNEVO runtime bundle production

The main GeneGalleon image carries the loader. A separately published runtime
contains an isolated Conda environment, the official native ANNEVO inference
code, and its original license notices. Model weights are separate first-use
downloads; their official URLs, sizes and SHA-256 values are bound to the same
resolved source commit. Training code, example data and weights do not enter the
runtime archive.

## Build contract

Resolve the official moving branch once at build start, with an explicit commit
override only for a one-off reproduction. Do not save resolved upstream commits
or dependency versions as repository defaults. Build the dedicated environment
inside the target GeneGalleon Linux architecture/base-ABI container. Packaging
tools belong to the builder, never to the sealed input environment.
Temporary archives and CUDA-library spooling use the bound build storage;
contained `/tmp` may be a small RAM filesystem and is not a bundle staging area.

The base ABI JSON records `system`, `machine`, `libc: {name, version}` and the
qualified GeneGalleon `image_digest`. The producer requires the environment's
CPU-hidden Python/Torch metadata to match this ABI and contain CUDA-built
Torch; a CPU-only environment cannot masquerade as a GPU profile. It also checks
all three native CLI helps. These checks do not constitute GPU qualification.

```bash
python container/scripts/build_annevo_runtime_bundle.py \
  --environment /build/annevo-environment \
  --source /build/official-ANNEVO \
  --source-ref main \
  --base-abi-json /build/genegalleon-base-abi.json \
  --profile linux-amd64-cuda-profile \
  --output /build/annevo-distribution
```

`conda-pack` must be available in the builder. Its missing-file and editable
checks remain enabled; there is no `--dest-prefix`, ignored missing file or
CPU fallback. If an old environment was modified by pip and its Conda metadata
is inconsistent, recreate the dedicated environment from the resolved build
profile. Do not weaken packaging checks or modify the original sealed runtime.

The current official profile makes Intel and LLVM OpenMP packages both own
`lib/libiomp5.so` and `lib/libomptarget.so`. Conda-pack reads their distinct
package-cache copies, which would create conflicting archive entries. For only
these two demonstrated collisions, the producer retains the effective physical
prefix files once, after their SHA matches a recorded provider and their type,
mode and absence of relocation placeholders pass verification. Conda-pack's
ordinary metadata checks still run before excluding those duplicate names;
provider/package metadata, licenses and the complete conflict proof are retained.
Unknown overlaps and conflicting entries remain errors. Remove this handling
when the official profile or upstream packaging no longer has these collisions,
after verifying ordinary strict packaging preserves the effective prefix.

For a private qualification copy, use a fresh target directory and regular
copies/reflinks, never hardlinks that share mutable inodes. Keep the original
inventory and all source bytes read-only. A copied old Conda environment still
needs strict packaging validation; copying does not repair its package metadata.

## Archive and delivery

The compressed tar contains `environment/`, native `source/`, and `provenance/`.
It includes `environment/bin/conda-unpack`, a package/dependency manifest, source
commit metadata, original license notices, and a complete pre-relocation payload
manifest. Unsafe paths, escaping links and special files are rejected. Each
public part defaults to 1 GiB; the manifest binds both every part and the complete
joined archive. OCI layer digests and sizes preserve these exact part identities.
Input PAX path/link overrides are validated before relocation; emitted PAX
headers are regenerated from the final names and normalized metadata. Before
publishing parts, the producer reads back every serialized header and file hash
and checks the exact inventory, links, modes and declared size. A private full
consumer extraction and CPU relocation check must pass before a new candidate's
GPU qualification; construction alone does not verify relocated imports.

The builder emits no speculative download URLs. A publisher adds real locators
only after durable upload and anonymous download/hash verification. The
intended tag-free delivery is a public GHCR OCI artifact, with an immutable
content-addressed locator and a moving latest index. CI publishes immutable
candidates only: the default alias is updated with `--promote-existing` after
that exact candidate's reserved GPU qualification passes. Publication defaults
to immutable-only; promotion requires both the actual GPU exit proof and cold
prepare receipt. The publisher checks the complete resolved manifest, source,
selected model, archive identity, CUDA profile, successful cold installation,
download-free warm reuse and byte-identical prediction outputs before changing
the alias:

```bash
python container/scripts/publish_annevo_runtime_bundle.py \
  --manifest /build/annevo-distribution/<exact-manifest>.json \
  --immutable-only
# After the exact immutable candidate's reserved GPU qualification succeeds:
python container/scripts/publish_annevo_runtime_bundle.py \
  --manifest /build/annevo-distribution/<exact-manifest>.json \
  --promote-existing \
  --gpu-qualification-proof /qualification/qualification_exit.json \
  --prepared-receipt /qualification/cold_prepare.json
```

Existing release-SIF
publication requires an existing release tag and is not an implicit release
authorization. Until a public runtime is actually published and verified, a
prebuilt default is not ready. Do not silently fall back to an older runtime,
another model, CPU inference or source installation after a prebuilt failure.

Official references: [conda-pack relocation caveats](https://conda.github.io/conda-pack/)
and [anonymous public GHCR pulls](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry).

## Install at its final address

Verify the manifest, part hashes, whole archive and raw payload before executing
anything. The transport verifies each opaque part against its declared SHA-256
and byte count, including reused and resumed files. A gzip prefix in a split
part does not make that part a complete gzip document. This explicit byte
contract never decompresses or normalizes part bytes; ordinary scientific
downloads retain their structural gzip/ZIP and HTML validation. Extract directly
into the final private cache directory. Run that
directory's Python with `environment/bin/conda-unpack`, then verify native
imports and metadata. Write readiness atomically only after verification; update
a small active pointer to publish it. **Never rename the runtime directory after
`conda-unpack`**: the repaired binary prefixes refer to that final path.

Fetch only the requested model from its source-bound metadata, verify its size
and SHA-256, and preserve the ANNEVO license/provenance. Keep outputs, Torch and
Numba caches, and temporary socket paths outside the runtime. Reserved GPU
qualification of the final unpacked prefix/image is separate from construction;
it must establish actual CUDA model execution rather than CPU fallback.

ANNEVO's current repository uses noncommercial terms. Preserve the original
license and copyright notices in bundles and model provenance; a combined
distribution must retain each component's terms. Do not relabel ANNEVO as MIT,
GPL, or unrestricted commercial software.
