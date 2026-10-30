# Readonly immutable launch code

An external launch producer may explicitly seal
`launch_storage=readonly-hardlinks-v1` into an immutable runtime manifest. Every
inventoried source file must then be owned by the executing user and have mode
0400, 0444 or 0555. Exact paths, inode, owner/group, mode, size, mtime, source
hashes, required closure and snapshot contract remain verified. Manifest and
directory identities retain their original full fences.

Concurrent creation of another package can add a hardlink to these read-only
code inodes. The bootstrap permits only the resulting link-count/ctime changes;
it hashes the body again after such a change during a read or after an earlier
inventory read. A changed body, restored mtime, changed permission, replaced
inode, redirected directory or leaf still stops startup. This is conditional on
the explicit storage contract; legacy manifests retain all original fences.

The bootstrap exposes `GG_LAUNCH_STORAGE_CAPABILITY=readonly-hardlinks-v1` for
producers to require before enabling the format. GeneGalleon creates no shared
cache, modifies no existing launch package and performs no object garbage
collection. Producers retain their ownership, publication, recovery and saved
request contracts. Execution settings and scientific outputs are unchanged.
