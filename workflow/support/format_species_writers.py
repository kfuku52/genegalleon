import bz2
import gzip
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

COMMON_REPLACEMENTS = (
    ("evm_27.model.", ""),
    ("evm.model.", ""),
    ("Oropetium_20150105_", ""),
)


def output_gzip_compresslevel():
    raw = os.environ.get("GG_INPUT_GZIP_LEVEL", "").strip()
    if raw == "":
        return 9
    try:
        level = int(raw)
    except ValueError as exc:
        raise ValueError("GG_INPUT_GZIP_LEVEL must be an integer between 1 and 9.") from exc
    if level < 1 or level > 9:
        raise ValueError("GG_INPUT_GZIP_LEVEL must be between 1 and 9.")
    return level


def output_compression_threads():
    raw = os.environ.get("GG_TASK_CPUS", "").strip()
    if raw == "":
        return 1
    try:
        value = int(raw)
    except ValueError:
        return 1
    return max(1, value)


def open_text(path, mode, errors="strict"):
    path = Path(path)
    if path.name.endswith(".gz"):
        kwargs = {"encoding": "utf-8", "errors": errors}
        if any(flag in mode for flag in ("w", "a", "x")):
            kwargs["compresslevel"] = output_gzip_compresslevel()
        return gzip.open(path, mode, **kwargs)
    if path.name.endswith(".bz2"):
        return bz2.open(path, mode, encoding="utf-8", errors=errors)
    return open(path, mode, encoding="utf-8", errors=errors)


def open_binary(path, mode="rb"):
    path = Path(path)
    if path.name.endswith(".gz"):
        return gzip.open(path, mode)
    if path.name.endswith(".bz2"):
        return bz2.open(path, mode)
    return open(path, mode)


def inspect_invalid_utf8(path, max_line_numbers=20):
    invalid_bytes = 0
    invalid_sequences = 0
    affected_lines = 0
    line_numbers = []
    with open_binary(path, "rb") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            offset = 0
            line_invalid = False
            while offset < len(raw_line):
                try:
                    raw_line[offset:].decode("utf-8")
                    break
                except UnicodeDecodeError as exc:
                    width = max(1, int(exc.end) - int(exc.start))
                    invalid_bytes += width
                    invalid_sequences += 1
                    line_invalid = True
                    offset += max(int(exc.end), int(exc.start) + 1)
            if line_invalid:
                affected_lines += 1
                if len(line_numbers) < max(0, int(max_line_numbers)):
                    line_numbers.append(line_number)
    return {
        "invalid_utf8_bytes": invalid_bytes,
        "invalid_utf8_sequences": invalid_sequences,
        "invalid_utf8_line_count": affected_lines,
        "invalid_utf8_lines": line_numbers,
    }


def make_temporary_output_path(output_path):
    base_name = output_path.name
    suffix = output_path.suffix
    if suffix != "" and base_name.endswith(suffix):
        base_name = base_name[: -len(suffix)]
    return output_path.parent / ".{}.tmp.{}.{}{}".format(
        base_name,
        os.getpid(),
        time.time_ns(),
        suffix,
    )


def write_text_output_via_command(output_path, writer, command_builder, output_via_stdout):
    tmp_output = make_temporary_output_path(output_path)
    stderr_text = ""
    proc = None
    stdout_handle = None
    try:
        stdout_target = subprocess.DEVNULL
        if output_via_stdout:
            stdout_handle = open(tmp_output, "wb")
            stdout_target = stdout_handle
        command = command_builder(tmp_output)
        proc = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stdout=stdout_target,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )
        try:
            if proc.stdin is None:
                raise RuntimeError("Failed to open stdin for: {}".format(" ".join(command)))
            writer(proc.stdin)
            proc.stdin.close()
            proc.stdin = None
            if proc.stderr is not None:
                stderr_text = proc.stderr.read()
            return_code = proc.wait()
        except Exception:
            if proc.stdin is not None:
                try:
                    proc.stdin.close()
                except Exception:
                    pass
            proc.kill()
            proc.wait()
            if proc.stderr is not None:
                try:
                    stderr_text = proc.stderr.read()
                except Exception:
                    stderr_text = ""
            raise
        finally:
            if stdout_handle is not None:
                stdout_handle.close()
                stdout_handle = None

        if return_code != 0:
            raise RuntimeError(
                "Command failed while writing '{}': {}".format(
                    output_path,
                    stderr_text.strip() or "exit code {}".format(return_code),
                )
            )
        if not tmp_output.exists() or tmp_output.stat().st_size == 0:
            raise RuntimeError("Command produced no output for '{}'".format(output_path))
        tmp_output.replace(output_path)
    finally:
        if proc is not None and proc.stderr is not None:
            proc.stderr.close()
        if stdout_handle is not None:
            stdout_handle.close()
        if tmp_output.exists():
            tmp_output.unlink()


def write_fasta_record(handle, record_id, sequence, width=80):
    handle.write(">{}\n".format(record_id))
    for idx in range(0, len(sequence), width):
        handle.write(sequence[idx : idx + width] + "\n")


def write_fasta_records_gzip(output_path, records):
    seqkit_path = shutil.which("seqkit")
    if seqkit_path is None:
        with open_text(output_path, "wt") as handle:
            for record_id, sequence in records:
                write_fasta_record(handle, record_id, sequence)
        return

    def writer(handle):
        for record_id, sequence in records:
            write_fasta_record(handle, record_id, sequence)

    write_text_output_via_command(
        output_path,
        writer,
        lambda tmp_output: [
            seqkit_path,
            "seq",
            "--threads",
            str(output_compression_threads()),
            "-o",
            str(tmp_output),
            "-",
        ],
        output_via_stdout=False,
    )


def apply_common_replacements(text):
    out = text
    for old, new in COMMON_REPLACEMENTS:
        out = out.replace(old, new)
    return out


def write_gff_gzip(input_path, output_path):
    pigz_path = shutil.which("pigz")
    line_count = 0

    if pigz_path is None:
        with open_text(input_path, "rt", errors="replace") as fin, open_text(output_path, "wt") as fout:
            for line in fin:
                fout.write(apply_common_replacements(line))
                line_count += 1
        return line_count

    def writer(handle):
        nonlocal line_count
        with open_text(input_path, "rt", errors="replace") as fin:
            for line in fin:
                handle.write(apply_common_replacements(line))
                line_count += 1

    write_text_output_via_command(
        output_path,
        writer,
        lambda _tmp_output: [
            pigz_path,
            "-p",
            str(output_compression_threads()),
            "-c",
        ],
        output_via_stdout=True,
    )
    return line_count


def write_gff_lines_gzip(output_path, lines):
    pigz_path = shutil.which("pigz")
    line_count = 0
    feature_count = 0

    if pigz_path is None:
        with open_text(output_path, "wt") as fout:
            for line in lines:
                normalized = apply_common_replacements(line)
                fout.write(normalized)
                line_count += 1
                if normalized.strip() != "" and not normalized.startswith("#"):
                    feature_count += 1
        return line_count, feature_count

    def writer(handle):
        nonlocal line_count, feature_count
        for line in lines:
            normalized = apply_common_replacements(line)
            handle.write(normalized)
            line_count += 1
            if normalized.strip() != "" and not normalized.startswith("#"):
                feature_count += 1

    write_text_output_via_command(
        output_path,
        writer,
        lambda _tmp_output: [
            pigz_path,
            "-p",
            str(output_compression_threads()),
            "-c",
        ],
        output_via_stdout=True,
    )
    return line_count, feature_count


def write_text_lines(output_path, lines):
    tmp_output = make_temporary_output_path(output_path)
    try:
        with open_text(tmp_output, "wt") as handle:
            for line in lines:
                handle.write(str(line))
        tmp_output.replace(output_path)
    except Exception:
        try:
            tmp_output.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            pass
        raise


def write_fasta_records(output_path, records):
    tmp_output = make_temporary_output_path(output_path)
    try:
        with open_text(tmp_output, "wt") as handle:
            for header, sequence in records:
                handle.write(">{}\n{}\n".format(str(header or "").strip(), re.sub(r"\s+", "", str(sequence or ""))))
        tmp_output.replace(output_path)
    except Exception:
        try:
            tmp_output.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            pass
        raise
