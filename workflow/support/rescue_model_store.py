"""Lossless, bounded rescue-model publications with selective verified readers.

Candidate order and every JSON field are retained. Shared bodies contain only
exactly equal predictor/DNA values; evidence and decisions remain per candidate.
The manifest is an index, not an independent trust anchor: verified readers bind
it and their consumed members to the enclosing producer receipt. A streaming
reader must be exhausted to complete its checksum and generation fence.
"""
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sqlite3
import tempfile
from collections import OrderedDict
from pathlib import Path, PurePosixPath

FORMAT = "gg_rescue_model_store_v1"
KINDS = {"models", "accepted", "partial", "revision"}
BODY_FIELDS = frozenset({
    "seqid", "strand", "cds", "sequence", "assembly_ambiguous_bases",
    "coverage", "identity", "query_span_coverage", "query_start", "query_end",
    "query_length", "frameshift", "paf_seqid", "paf_strand", "search",
})
MAX_RECORD_BYTES = 64 * 1024 * 1024
_SHA = re.compile(r"[0-9a-f]{64}\Z")
_PAF = re.compile(r"\A(##?PAF\t)?([^\t\r\n]+)(\t[^\r\n]*)(\r?\n)?\Z")


def _json(value, *, sort=False):
    return json.dumps(value, ensure_ascii=False, allow_nan=False,
                      separators=(",", ":"), sort_keys=sort).encode("utf-8")


def _codec(value):
    if value not in {"auto", "gzip", "zstd"}:
        raise ValueError("Unsupported model-store compression")
    if value != "gzip":
        try:
            import zstandard
        except ImportError:
            if value == "zstd":
                raise ValueError("Model-store zstd codec requires zstandard") from None
        else:
            return "zstd", zstandard.__version__
    return "gzip", "stdlib"


def _compress(value, codec, level):
    if codec == "gzip":
        return gzip.compress(value, compresslevel=level, mtime=0)
    import zstandard
    return zstandard.ZstdCompressor(level=level, write_checksum=True).compress(value)


def _decompress(value, codec):
    if codec == "gzip":
        reader = gzip.GzipFile(fileobj=io.BytesIO(value), mode="rb")
    else:
        import zstandard
        reader = zstandard.ZstdDecompressor().stream_reader(io.BytesIO(value))
    with reader:
        result = reader.read(MAX_RECORD_BYTES + 1)
    if len(result) > MAX_RECORD_BYTES:
        raise ValueError("Model-store record exceeds bounded buffer")
    return result


def _writer(path, codec, level):
    if codec == "gzip":
        return gzip.GzipFile(filename=path, mode="wb", compresslevel=level, mtime=0)
    import zstandard
    return zstandard.ZstdCompressor(level=level, write_checksum=True).stream_writer(Path(path).open("wb"))


def _digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _sync_directory(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _safe_member(value):
    if (not isinstance(value, str) or not value or "\\" in value
            or any(ord(c) < 32 for c in value)):
        raise ValueError("Unsafe model-store member")
    p = PurePosixPath(value)
    if p.is_absolute() or str(p) != value or any(c in {"", ".", ".."} for c in p.parts):
        raise ValueError("Unsafe model-store member")
    return value


def _member(root, value):
    value = _safe_member(value)
    current = root
    for part in PurePosixPath(value).parts:
        current = current / part
        if current.is_symlink():
            raise ValueError("Symlink model-store member")
    if not current.is_file():
        raise ValueError(f"Missing model-store member: {value}")
    return current


def _snapshot(path):
    raw = Path(path).read_bytes()
    try:
        value = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("Malformed model-store metadata") from error
    return value, hashlib.sha256(raw).hexdigest()


def _layout(directory):
    directory = Path(directory).resolve()
    if (directory / "manifest.json").is_file():
        return directory.parent, directory
    return directory, directory / "model_store"


def _split_fields(record):
    shared, context = {}, {}
    for key, value in record.items():
        if key in BODY_FIELDS:
            shared[key] = value
        elif key == "paf" and isinstance(value, str) and (match := _PAF.fullmatch(value)):
            shared[key] = [match[1], match[3], match[4]]
            context[key] = {"paf_query_token": match[2]}
        else:
            context[key] = value
    return shared, context, list(record)


def _pack(record):
    if not isinstance(record, dict) or any(not isinstance(k, str) for k in record):
        raise ValueError("Model-store records must be JSON objects")
    shared, context, order = _split_fields(record)
    raw_order = None
    if isinstance(record.get("raw_prediction"), dict):
        raw_body, raw_context, raw_order = _split_fields(record["raw_prediction"])
        shared["raw_prediction"] = raw_body
        context["raw_prediction"] = raw_context
    body = _json(shared, sort=True)
    if len(body) > MAX_RECORD_BYTES:
        raise ValueError("Model-store body exceeds bounded buffer")
    envelope = {"body": hashlib.sha256(body).hexdigest(), "context": context, "order": order}
    if raw_order is not None:
        envelope["raw_order"] = raw_order
    # Validate candidate-only values now as well, including nonfinite numbers.
    _json(envelope)
    return body, envelope


def _restore_fields(shared, context, order):
    if (not isinstance(shared, dict) or not isinstance(context, dict) or not isinstance(order, list)
            or any(not isinstance(k, str) for k in order) or len(set(order)) != len(order)
            or set(order) != set(shared) | set(context)):
        raise ValueError("Malformed model-store field layout")
    result = {}
    for key in order:
        if key == "paf" and key in shared and key in context:
            parts, token = shared[key], context[key]
            if (not isinstance(parts, list) or len(parts) != 3
                    or not isinstance(token, dict) or set(token) != {"paf_query_token"}
                    or not isinstance(token["paf_query_token"], str)
                    or parts[0] not in {None, "#PAF\t", "##PAF\t"} or not isinstance(parts[1], str)
                    or parts[2] not in {None, "\n", "\r\n"}):
                raise ValueError("Malformed shared PAF binding")
            result[key] = (parts[0] or "") + token["paf_query_token"] + parts[1] + (parts[2] or "")
        elif key in shared and key in context:
            raise ValueError("Overlapping model-store fields")
        else:
            result[key] = shared[key] if key in shared else context[key]
    return result


def _unpack(envelope, body):
    if (not isinstance(envelope, dict) or set(envelope) - {"body", "context", "order", "raw_order", "ordinal"}
            or not _SHA.fullmatch(str(envelope.get("body", "")))):
        raise ValueError("Malformed model-store envelope")
    # Copies must be independent: callers routinely modify CDS/evidence/status.
    shared = json.loads(body)
    context = json.loads(_json(envelope["context"]))
    if "raw_order" in envelope:
        if "raw_prediction" not in shared or "raw_prediction" not in context:
            raise ValueError("Missing raw-prediction field layout")
        context["raw_prediction"] = _restore_fields(
            shared.pop("raw_prediction"), context["raw_prediction"], envelope["raw_order"])
    return _restore_fields(shared, context, envelope["order"])


def _database(path):
    db = sqlite3.connect(path)
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("CREATE TABLE bodies (sha TEXT PRIMARY KEY, data BLOB NOT NULL) WITHOUT ROWID")
    return db


def _insert_body(db, body, level, codec, cache=None):
    sha = hashlib.sha256(body).hexdigest()
    if cache is not None and sha in cache.cache:
        if cache.cache[sha] != body:
            raise ValueError("Model-store body hash collision")
        cache.cache.move_to_end(sha)
        return sha
    found = db.execute("SELECT data FROM bodies WHERE sha=?", (sha,)).fetchone()
    if found is None:
        db.execute("INSERT INTO bodies VALUES (?,?)", (sha, _compress(body, codec, level)))
    elif _decompress(found[0], codec) != body:
        raise ValueError("Model-store body hash collision")
    if cache is not None:
        cache.remember(sha, body)
    return sha


def _overlay(path):
    db = _database(path)
    db.execute("CREATE TABLE entries (ordinal INTEGER PRIMARY KEY, envelope BLOB NOT NULL)")
    return db


def _add_overlay(db, ordinal, body, envelope, level, codec, cache=None):
    _insert_body(db, body, level, codec, cache)
    db.execute("INSERT INTO entries VALUES (?,?)", (ordinal, _compress(_json(envelope), codec, level)))


def write_model_store(directory, models, *, partial_models=None, revisions=None,
                      shard_bytes=8 * 1024 * 1024, codec="auto", compression_level=1):
    """Atomically publish into ``directory/model_store`` without replacing data.

    ``partial_models=None`` records references for ``partial_evidence.partial``.
    Explicit partial records must be exact members of models. Revision records
    can differ from models and are stored in an independent, small overlay.
    Publication is immutable; the enclosing stage subsequently hashes all files
    into its normal receipt. No receipt is fabricated by this writer.
    """
    codec, codec_version = _codec(codec)
    if type(compression_level) is not int or not 0 <= compression_level <= (9 if codec == "gzip" else 22):
        raise ValueError("Unsupported model-store compression")
    if type(shard_bytes) is not int or shard_bytes < 1:
        raise ValueError("Invalid model-store shard bound")
    directory = Path(directory)
    target = directory if directory.name == "model_store" else directory / "model_store"
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        raise ValueError("Model-store publication already exists")
    tmp = Path(tempfile.mkdtemp(prefix=".model-store-", dir=target.parent))
    dbs, shard_handle = [], None
    try:
        body_db = _database(tmp / "bodies.sqlite")
        lookup_db = None
        if partial_models is not None:
            # Needed only to resolve an explicit, possibly reordered legacy
            # partial stream. Never keep all-candidate lookup metadata in the
            # published DNA-body database or build it for normal producers.
            lookup_db = sqlite3.connect(tmp / ".lookup.sqlite")
            lookup_db.execute("PRAGMA journal_mode=OFF")
            lookup_db.execute("PRAGMA synchronous=OFF")
            lookup_db.execute("CREATE TABLE model_refs (fingerprint BLOB, ordinal INTEGER PRIMARY KEY, shard INTEGER, line INTEGER)")
            lookup_db.execute("CREATE INDEX fingerprints ON model_refs(fingerprint)")
        accepted_db = _overlay(tmp / "accepted.sqlite")
        revision_db = _overlay(tmp / "revision.sqlite")
        partial_db = sqlite3.connect(tmp / "partial.sqlite")
        partial_db.execute("CREATE TABLE refs (position INTEGER PRIMARY KEY, ordinal INTEGER, shard INTEGER, line INTEGER)")
        partial_db.execute("CREATE INDEX locations ON refs(shard,line)")
        dbs = [body_db, accepted_db, revision_db, partial_db]
        if lookup_db is not None:
            dbs.append(lookup_db)
        body_cache = _Bodies(body_db, codec)
        accepted_cache = _Bodies(accepted_db, codec)
        revision_cache = _Bodies(revision_db, codec)
        shards, shard_size, line, count, accepted_count, partial_count = [], 0, 0, 0, 0, 0
        partial_ordered, last_partial_ordinal = True, -1
        for ordinal, record in enumerate(models):
            body, envelope = _pack(record)
            if lookup_db is not None:
                fingerprint = hashlib.sha256(_json(envelope)).digest()
            envelope["ordinal"] = ordinal
            packed = _json(envelope) + b"\n"
            if len(packed) > MAX_RECORD_BYTES:
                raise ValueError("Model-store record exceeds bounded buffer")
            if shard_handle is None or (shard_size and shard_size + len(packed) > shard_bytes):
                if shard_handle is not None:
                    shard_handle.close()
                extension = "gz" if codec == "gzip" else "zst"
                shard = f"models-{len(shards):06d}.jsonl.{extension}"
                shards.append({"path": shard, "first": ordinal, "count": 0})
                shard_handle = _writer(tmp / shard, codec, compression_level)
                shard_size, line = 0, 0
            shard_handle.write(packed)
            _insert_body(body_db, body, compression_level, codec, body_cache)
            if lookup_db is not None:
                lookup_db.execute("INSERT INTO model_refs VALUES (?,?,?,?)",
                                  (fingerprint, ordinal, len(shards) - 1, line))
            if record.get("status") == "accepted":
                _add_overlay(accepted_db, ordinal, body, envelope, compression_level, codec, accepted_cache)
                accepted_count += 1
            if partial_models is None and isinstance(record.get("partial_evidence"), dict) and record["partial_evidence"].get("partial"):
                partial_db.execute("INSERT INTO refs VALUES (?,?,?,?)", (partial_count, ordinal, len(shards) - 1, line))
                partial_count += 1
            shard_size += len(packed)
            line += 1
            shards[-1]["count"] += 1
            count += 1
        if shard_handle is not None:
            shard_handle.close()
            shard_handle = None
        if partial_models is not None:
            for record in partial_models:
                body, envelope = _pack(record)
                fingerprint = hashlib.sha256(_json(envelope)).digest()
                match = lookup_db.execute("SELECT ordinal,shard,line FROM model_refs WHERE fingerprint=? ORDER BY ordinal LIMIT 1", (fingerprint,)).fetchone()
                if match is None:
                    raise ValueError("Partial record is not an exact model-store member")
                # Verify equal body bytes, not only equal SHA, before referencing.
                _insert_body(body_db, body, compression_level, codec, body_cache)
                if match[0] < last_partial_ordinal:
                    partial_ordered = False
                last_partial_ordinal = match[0]
                partial_db.execute("INSERT INTO refs VALUES (?,?,?,?)", (partial_count, *match))
                partial_count += 1
            lookup_db.close()
            dbs.remove(lookup_db)
            (tmp / ".lookup.sqlite").unlink()
        revision_count = 0
        for ordinal, record in enumerate(revisions or ()):
            body, envelope = _pack(record)
            envelope["ordinal"] = ordinal
            _add_overlay(revision_db, ordinal, body, envelope, compression_level, codec, revision_cache)
            revision_count += 1
        body_count = body_db.execute("SELECT COUNT(*) FROM bodies").fetchone()[0]
        partial_shards = [shards[r[0]]["path"] for r in partial_db.execute("SELECT DISTINCT shard FROM refs ORDER BY shard")]
        for db in dbs:
            db.commit()
            db.close()
        dbs = []
        files = {p.name: _digest(p) for p in sorted(tmp.iterdir())}
        scopes = {
            "models": ["bodies.sqlite", *(s["path"] for s in shards)],
            "accepted": ["accepted.sqlite"],
            "revision": ["revision.sqlite"],
            "partial": ["partial.sqlite", *(["bodies.sqlite", *partial_shards] if partial_count else [])],
        }
        manifest = {"schema": 1, "format": FORMAT, "codec": codec, "body_codec": codec,
                    "codec_version": codec_version,
                    "compression_level": compression_level, "files": files, "scopes": scopes,
                    "shards": shards, "partial_ordered": partial_ordered,
                    "counts": {"models": count, "accepted": accepted_count,
                                                 "partial": partial_count, "revision": revision_count,
                                                 "bodies": body_count}}
        (tmp / "manifest.json").write_bytes(_json(manifest) + b"\n")
        for path in tmp.iterdir():
            with path.open("rb") as handle:
                os.fsync(handle.fileno())
        os.rename(tmp, target)
        _sync_directory(target.parent)
        return manifest
    except BaseException:
        if shard_handle is not None:
            shard_handle.close()
        for db in dbs:
            db.close()
        shutil.rmtree(tmp, ignore_errors=True)
        raise


def _manifest(store):
    manifest, sha = _snapshot(store / "manifest.json")
    if (not isinstance(manifest, dict) or type(manifest.get("schema")) is not int
            or manifest.get("schema") != 1 or manifest.get("format") != FORMAT
            or manifest.get("codec") not in {"gzip", "zstd"} or manifest.get("body_codec") != manifest.get("codec")
            or not isinstance(manifest.get("files"), dict) or set(manifest.get("scopes", {})) != KINDS
            or not isinstance(manifest.get("counts"), dict) or not isinstance(manifest.get("shards"), list)):
        raise ValueError("Unknown or malformed model-store manifest")
    for member, digest in manifest["files"].items():
        _safe_member(member)
        if not isinstance(digest, str) or not _SHA.fullmatch(digest):
            raise ValueError("Malformed model-store digest")
    for kind, members in manifest["scopes"].items():
        if (not isinstance(members, list) or len(members) != len(set(members))
                or any(p not in manifest["files"] for p in members)):
            raise ValueError("Malformed model-store scope")
        if type(manifest["counts"].get(kind)) is not int or manifest["counts"][kind] < 0:
            raise ValueError("Malformed model-store count")
    _codec(manifest["codec"])
    ordinal = 0
    shard_names = []
    for shard in manifest["shards"]:
        if (not isinstance(shard, dict) or set(shard) != {"path", "first", "count"}
                or shard["path"] not in manifest["files"] or shard["path"] in shard_names
                or type(shard["first"]) is not int or shard["first"] != ordinal
                or type(shard["count"]) is not int or shard["count"] < 1):
            raise ValueError("Malformed model-store shard index")
        shard_names.append(shard["path"])
        ordinal += shard["count"]
    if (ordinal != manifest["counts"]["models"] or manifest["scopes"]["models"] != ["bodies.sqlite", *shard_names]
            or manifest["scopes"]["accepted"] != ["accepted.sqlite"]
            or manifest["scopes"]["revision"] != ["revision.sqlite"]
            or not manifest["scopes"]["partial"] or manifest["scopes"]["partial"][0] != "partial.sqlite"
            or any(p not in {"partial.sqlite", "bodies.sqlite", *shard_names} for p in manifest["scopes"]["partial"])):
        raise ValueError("Malformed model-store scope dependencies")
    return manifest, sha


def _receipt(root):
    path = root / "receipt.json"
    if not path.is_file() or path.is_symlink():
        raise ValueError("Model-store producer receipt is required")
    receipt, sha = _snapshot(path)
    if not isinstance(receipt, dict) or not isinstance(receipt.get("key"), dict) or not isinstance(receipt.get("files"), dict):
        raise ValueError("Malformed model-store producer receipt")
    for member, digest in receipt["files"].items():
        _safe_member(member)
        if not isinstance(digest, str) or not _SHA.fullmatch(digest):
            raise ValueError("Malformed producer receipt digest")
    return receipt, sha


def frozen_model_store_key(directory, *, kind="models"):
    """Cheap receipt-bound key; hash selected contents at verification/iteration."""
    if kind not in KINDS:
        raise ValueError("Unknown model-store kind")
    root, store = _layout(directory)
    receipt, receipt_sha = _receipt(root)
    if (store / "manifest.json").is_file():
        manifest, manifest_sha = _manifest(store)
        prefix = store.relative_to(root).as_posix()
        files = {f"{prefix}/manifest.json": manifest_sha}
        files.update({f"{prefix}/{p}": manifest["files"][p] for p in manifest["scopes"][kind]})
        if any(receipt["files"].get(p) != sha for p, sha in files.items()):
            raise ValueError("Model-store manifest/members are not bound to producer receipt")
        return {"schema": 1, "format": FORMAT, "kind": kind, "manifest_sha256": manifest_sha,
                "producer_receipt_sha256": receipt_sha, "files": files, "absent": []}
    member = {"models": "models.json", "accepted": "models.json", "partial": "partial_models.json",
              "revision": "revision_candidates.json"}[kind]
    files, absent = {}, []
    if member in receipt["files"]:
        files[member] = receipt["files"][member]
    elif kind in {"partial", "revision"} and not (root / member).exists() and not (root / member).is_symlink():
        absent.append(member)
    else:
        raise ValueError("Legacy model file is not bound to producer receipt")
    return {"schema": 1, "format": "legacy_json_array", "kind": kind, "manifest_sha256": None,
            "producer_receipt_sha256": receipt_sha, "files": files,
            "absent": ["model_store/manifest.json", *absent]}


def _check_key(directory, key, kind=None):
    if (not isinstance(key, dict) or type(key.get("schema")) is not int
            or key.get("schema") != 1 or key.get("kind") not in KINDS
            or (kind is not None and key["kind"] != kind)):
        raise ValueError("Malformed frozen model-store key")
    if frozen_model_store_key(directory, kind=key["kind"]) != key:
        raise ValueError("Model-store publication generation changed")


def verify_model_store_key(directory, key):
    """Return True after hashing every consumed member; reject changed generations."""
    _check_key(directory, key)
    root, _ = _layout(directory)
    for member, expected in key["files"].items():
        if _digest(_member(root, member)) != expected:
            raise ValueError(f"Model-store checksum changed: {member}")
    _check_key(directory, key)
    return True


def _readonly(path):
    return sqlite3.connect(Path(path).as_uri() + "?mode=ro&immutable=1", uri=True)


class _Bodies:
    def __init__(self, db, codec):
        self.db = db
        self.codec = codec
        self.cache = OrderedDict()
        self.bytes = 0

    def read(self, sha):
        if sha in self.cache:
            self.cache.move_to_end(sha)
            return self.cache[sha]
        row = self.db.execute("SELECT data FROM bodies WHERE sha=?", (sha,)).fetchone()
        if row is None:
            raise ValueError("Missing model-store body")
        body = _decompress(row[0], self.codec)
        if len(body) > MAX_RECORD_BYTES or hashlib.sha256(body).hexdigest() != sha:
            raise ValueError("Malformed model-store body")
        self.remember(sha, body)
        return body

    def remember(self, sha, body):
        if len(body) <= 8 * 1024 * 1024:
            self.cache[sha] = body
            self.bytes += len(body)
            while len(self.cache) > 2048 or self.bytes > 16 * 1024 * 1024:
                _, old = self.cache.popitem(last=False)
                self.bytes -= len(old)


class _HashReader:
    def __init__(self, handle):
        self.handle = handle
        self.hash = hashlib.sha256()

    def read(self, size=-1):
        value = self.handle.read(size)
        self.hash.update(value)
        return value

    def readinto(self, buffer):
        value = self.read(len(buffer))
        buffer[:len(value)] = value
        return len(value)

    def close(self):
        self.handle.close()


def _lines(path, codec, expected=None):
    with Path(path).open("rb") as raw:
        reader = _HashReader(raw)
        if codec == "gzip":
            handle = gzip.GzipFile(fileobj=reader, mode="rb")
        else:
            import zstandard
            handle = io.BufferedReader(zstandard.ZstdDecompressor().stream_reader(reader))
        with handle:
            while line := handle.readline(MAX_RECORD_BYTES + 1):
                if len(line) > MAX_RECORD_BYTES or not line.endswith(b"\n"):
                    raise ValueError("Malformed or oversized model-store JSONL record")
                try:
                    yield json.loads(line)
                except (UnicodeDecodeError, json.JSONDecodeError) as error:
                    raise ValueError("Malformed model-store JSONL record") from error
        if expected is not None and reader.hash.hexdigest() != expected:
            raise ValueError("Model-store shard checksum changed")


def _stat(path):
    st = path.stat()
    return st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns


def _iter(directory, kind, *, statuses=None, verify=True, frozen_key=None):
    root, store = _layout(directory)
    key = frozen_key or (frozen_model_store_key(directory, kind=kind) if verify else None)
    if key is not None:
        _check_key(directory, key, kind)
    paths = {p: _member(root, p) for p in key["files"]} if key is not None else {}
    stats = {p: _stat(path) for p, path in paths.items()}
    if not (store / "manifest.json").is_file():
        # Lazy import avoids a cycle with prediction-cache's new-store reader.
        try:
            from rescue_prediction_cache import stream_json_array
        except ImportError:
            from .rescue_prediction_cache import stream_json_array
        member = {"models": "models.json", "accepted": "models.json", "partial": "partial_models.json",
                  "revision": "revision_candidates.json"}[kind]
        if (root / member).is_file():
            hasher = hashlib.sha256()
            for record in stream_json_array(root / member, hasher=hasher):
                if not isinstance(record, dict):
                    raise ValueError("Legacy model record is not an object")
                if ((kind != "accepted" or record.get("status") == "accepted")
                        and (statuses is None or record.get("status") in statuses)):
                    yield record
            if key is not None and hasher.hexdigest() != key["files"][member]:
                raise ValueError("Legacy model checksum changed")
    else:
        manifest, _ = _manifest(store)
        prefix = store.relative_to(root).as_posix()
        if key is not None:
            for member, path in paths.items():
                if not member.endswith((".jsonl.gz", ".jsonl.zst")) and _digest(path) != key["files"][member]:
                    raise ValueError(f"Model-store checksum changed: {member}")
        if kind in {"accepted", "revision"}:
            db = _readonly(store / f"{kind}.sqlite")
            try:
                bodies = _Bodies(db, manifest["body_codec"])
                count = 0
                for ordinal, blob in db.execute("SELECT ordinal,envelope FROM entries ORDER BY ordinal"):
                    envelope = json.loads(_decompress(blob, manifest["codec"]))
                    if envelope.get("ordinal") != ordinal:
                        raise ValueError("Model-store overlay order changed")
                    record = _unpack(envelope, bodies.read(envelope["body"]))
                    if kind == "accepted" and record.get("status") != "accepted":
                        raise ValueError("Non-primary record in accepted overlay")
                    count += 1
                    if statuses is None or record.get("status") in statuses:
                        yield record
                if count != manifest["counts"][kind]:
                    raise ValueError("Model-store overlay count changed")
            finally:
                db.close()
        elif kind == "partial":
            yield from _partials(root, store, manifest, key)
        else:
            db = _readonly(store / "bodies.sqlite")
            try:
                bodies, ordinal = _Bodies(db, manifest["body_codec"]), 0
                for shard in manifest["shards"]:
                    count = 0
                    relative = f"{prefix}/{shard['path']}"
                    path = _member(root, relative)
                    expected = key["files"][relative] if key is not None else None
                    if shard["first"] != ordinal:
                        raise ValueError("Model-store shard order changed")
                    for envelope in _lines(path, manifest["codec"], expected):
                        if envelope.get("ordinal") != ordinal:
                            raise ValueError("Model-store candidate order changed")
                        record = _unpack(envelope, bodies.read(envelope["body"]))
                        ordinal += 1
                        count += 1
                        if statuses is None or record.get("status") in statuses:
                            yield record
                    if count != shard["count"]:
                        raise ValueError("Model-store shard count changed")
                if ordinal != manifest["counts"]["models"]:
                    raise ValueError("Model-store candidate count changed")
            finally:
                db.close()
    if key is not None:
        _check_key(directory, key, kind)
        if any(_stat(path) != stats[member] for member, path in paths.items()):
            raise ValueError("Model-store member changed during consumption")


def _partials(root, store, manifest, key):
    refs = _readonly(store / "partial.sqlite")
    body_db = None
    try:
        count = refs.execute("SELECT COUNT(*) FROM refs").fetchone()[0]
        if count != manifest["counts"]["partial"]:
            raise ValueError("Model-store partial count changed")
        if not count:
            return
        body_db = _readonly(store / "bodies.sqlite")
        bodies = _Bodies(body_db, manifest["body_codec"])
        if manifest.get("partial_ordered") is True:
            position = 0
            for found_position, envelope in _partial_entries(root, store, manifest, key, refs):
                if found_position != position:
                    raise ValueError("Model-store partial order changed")
                position += 1
                yield _unpack(envelope, bodies.read(envelope["body"]))
            if position != count:
                raise ValueError("Missing model-store partial reference")
            return
        # Explicit partial input may reorder or repeat models. Spool only compact
        # envelopes, on temporary disk, to recover exact partial order boundedly.
        with tempfile.TemporaryDirectory(prefix="gg-model-partial-") as tmp:
            spool = sqlite3.connect(Path(tmp) / "partial.sqlite")
            try:
                spool.execute("CREATE TABLE entries (position INTEGER PRIMARY KEY,envelope BLOB)")
                for position, envelope in _partial_entries(root, store, manifest, key, refs):
                    spool.execute("INSERT INTO entries VALUES (?,?)", (position, _json(envelope)))
                if spool.execute("SELECT COUNT(*) FROM entries").fetchone()[0] != count:
                    raise ValueError("Missing model-store partial reference")
                for (_, blob) in spool.execute("SELECT position,envelope FROM entries ORDER BY position"):
                    envelope = json.loads(blob)
                    yield _unpack(envelope, bodies.read(envelope["body"]))
            finally:
                spool.close()
    finally:
        refs.close()
        if body_db is not None:
            body_db.close()


def _partial_entries(root, store, manifest, key, refs):
    prefix = store.relative_to(root).as_posix()
    for (shard_index,) in refs.execute("SELECT DISTINCT shard FROM refs ORDER BY shard"):
        if type(shard_index) is not int or not 0 <= shard_index < len(manifest["shards"]):
            raise ValueError("Malformed model-store partial shard reference")
        shard = manifest["shards"][shard_index]["path"]
        relative = f"{prefix}/{_safe_member(shard)}"
        expected = key["files"][relative] if key is not None else None
        references = refs.execute(
            "SELECT line,position,ordinal FROM refs WHERE shard=? ORDER BY line,position", (shard_index,))
        reference = next(references, None)
        count = 0
        for line, envelope in enumerate(_lines(_member(root, relative), manifest["codec"], expected)):
            count += 1
            if (not isinstance(envelope, dict) or type(envelope.get("ordinal")) is not int
                    or envelope["ordinal"] != manifest["shards"][shard_index]["first"] + line):
                raise ValueError("Model-store partial shard order changed")
            if reference is not None and (type(reference[0]) is not int or reference[0] < line):
                raise ValueError("Malformed model-store partial line reference")
            while reference is not None and reference[0] == line:
                _, position, ordinal = reference
                if (type(position) is not int or position < 0 or type(ordinal) is not int
                        or envelope["ordinal"] != ordinal):
                    raise ValueError("Model-store partial reference changed")
                yield position, envelope
                reference = next(references, None)
        if count != manifest["shards"][shard_index]["count"]:
            raise ValueError("Model-store partial shard count changed")
        if reference is not None:
            raise ValueError("Missing model-store partial line reference")


def iter_models(directory, *, statuses=None, verify=True, frozen_key=None):
    return _iter(directory, "models", statuses=statuses, verify=verify, frozen_key=frozen_key)


def iter_accepted_models(directory, *, verify=True, frozen_key=None):
    return _iter(directory, "accepted", verify=verify, frozen_key=frozen_key)


def iter_partial_models(directory, *, verify=True, frozen_key=None):
    return _iter(directory, "partial", verify=verify, frozen_key=frozen_key)


def iter_revision_models(directory, *, verify=True, frozen_key=None):
    return _iter(directory, "revision", verify=verify, frozen_key=frozen_key)


def export_legacy_models(source, destination, *, kind="models", verify=True, frozen_key=None, directory=False):
    """Atomically export an exact ordered legacy array without loading it all."""
    if directory:
        return _export_directory(source, destination, verify=verify)
    if kind not in KINDS:
        raise ValueError("Unknown model-store kind")
    destination = Path(destination)
    root, store = _layout(source)
    if destination.resolve().is_relative_to(store):
        raise ValueError("Legacy export cannot overwrite the immutable model store")
    if destination.resolve() == root / "receipt.json":
        raise ValueError("Legacy export cannot overwrite the producer receipt")
    key = frozen_key or (frozen_model_store_key(source, kind=kind) if verify else None)
    if key is not None and destination.resolve() in {_member(root, p).resolve() for p in key["files"]}:
        raise ValueError("Legacy export cannot overwrite consumed source members")
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{destination.name}-", dir=destination.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(b"[")
            first = True
            for record in _iter(source, kind, verify=verify, frozen_key=key):
                if not first:
                    handle.write(b",")
                handle.write(_json(record))
                first = False
            handle.write(b"]\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, destination)
        _sync_directory(destination.parent)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _export_directory(source, destination, *, verify):
    destination = Path(destination)
    root, publication = _layout(source)
    if destination.resolve() == root or destination.resolve().is_relative_to(publication):
        raise ValueError("Legacy export cannot overwrite source publication")
    if destination.exists() or destination.is_symlink():
        raise ValueError("Legacy export destination already exists")
    destination.parent.mkdir(parents=True, exist_ok=True)
    keys = {kind: frozen_model_store_key(source, kind=kind) for kind in ("models", "partial", "revision")}
    tmp = Path(tempfile.mkdtemp(prefix=".model-export-", dir=destination.parent))
    try:
        names = {"models": "models.json", "partial": "partial_models.json", "revision": "revision_candidates.json"}
        for kind, name in names.items():
            export_legacy_models(source, tmp / name, kind=kind, verify=verify, frozen_key=keys[kind])
        # Source checksums and metadata are fenced after the last exported row.
        for key in keys.values():
            verify_model_store_key(source, key)
        export_receipt = {"schema": 1, "format": "gg_rescue_legacy_export_v1", "source": str(root),
                          "source_keys": keys, "files": {name: _digest(tmp / name) for name in names.values()}}
        (tmp / "export_receipt.json").write_bytes(_json(export_receipt) + b"\n")
        with (tmp / "export_receipt.json").open("rb") as handle:
            os.fsync(handle.fileno())
        os.rename(tmp, destination)
        _sync_directory(destination.parent)
        return export_receipt
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
