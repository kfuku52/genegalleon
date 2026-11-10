"""Supervise shell commands started by workflow integration tests."""

import os
import signal
import subprocess
from contextlib import suppress


def run_workflow_command(command, *, cwd, env, timeout):
    """Capture the command output and terminate its process group on timeout."""
    with subprocess.Popen(
        command, cwd=cwd, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True,
    ) as process:
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            # The task scratch launcher forwards TERM to its own child session.
            with suppress(ProcessLookupError):
                os.killpg(process.pid, signal.SIGTERM)
            try:
                process.communicate(timeout=2)
            except (subprocess.TimeoutExpired, UnicodeError):
                pass
            finally:
                with suppress(ProcessLookupError):
                    os.killpg(process.pid, signal.SIGKILL)
            raise
        return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
