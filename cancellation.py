"""Cooperative cancellation for Aurora pipeline subprocesses."""
from __future__ import annotations

import os
import signal
import subprocess
import threading
import time

_CANCEL_EVENT = threading.Event()
_ACTIVE_LOCK = threading.Lock()
_ACTIVE_PROCESSES = set()


class PipelineCancelled(RuntimeError):
    """Raised when the user cancels the active Aurora run."""


def reset_cancellation():
    _CANCEL_EVENT.clear()


def cancellation_requested():
    return _CANCEL_EVENT.is_set()


def check_cancelled():
    if _CANCEL_EVENT.is_set():
        raise PipelineCancelled("Pipeline run cancelled by user.")


def _terminate_process(proc, grace_seconds=2.0):
    if proc is None or proc.poll() is not None:
        return
    try:
        if os.name == "posix":
            os.killpg(proc.pid, signal.SIGTERM)
        else:
            proc.terminate()
    except (ProcessLookupError, PermissionError, OSError):
        try:
            proc.terminate()
        except Exception:
            pass

    deadline = time.monotonic() + grace_seconds
    while proc.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)

    if proc.poll() is None:
        try:
            if os.name == "posix":
                os.killpg(proc.pid, signal.SIGKILL)
            else:
                proc.kill()
        except (ProcessLookupError, PermissionError, OSError):
            try:
                proc.kill()
            except Exception:
                pass


def request_cancellation():
    """Signal cancellation and terminate all currently tracked child processes."""
    _CANCEL_EVENT.set()
    with _ACTIVE_LOCK:
        processes = list(_ACTIVE_PROCESSES)
    for proc in processes:
        _terminate_process(proc)


def run_cancellable(
    cmd,
    *,
    check=False,
    stdout=None,
    stderr=None,
    capture_output=False,
    text=False,
    env=None,
    timeout=None,
    **kwargs,
):
    """subprocess.run-like helper that can be interrupted from the GUI.

    Child tools are started in their own process group on POSIX so cancellation
    also stops descendants launched by Medaka or shell wrappers.
    """
    check_cancelled()

    if capture_output:
        if stdout is not None or stderr is not None:
            raise ValueError("stdout/stderr may not be used with capture_output")
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE

    proc = subprocess.Popen(
        cmd,
        stdout=stdout,
        stderr=stderr,
        text=text,
        env=env,
        start_new_session=(os.name == "posix"),
        **kwargs,
    )
    with _ACTIVE_LOCK:
        _ACTIVE_PROCESSES.add(proc)

    started = time.monotonic()
    out = err = None
    try:
        while True:
            if cancellation_requested():
                _terminate_process(proc)
                raise PipelineCancelled("Pipeline run cancelled by user.")

            if timeout is not None and (time.monotonic() - started) >= timeout:
                _terminate_process(proc)
                raise subprocess.TimeoutExpired(cmd, timeout, output=out, stderr=err)

            try:
                out, err = proc.communicate(timeout=0.2)
                if cancellation_requested():
                    raise PipelineCancelled("Pipeline run cancelled by user.")
                break
            except subprocess.TimeoutExpired:
                continue
    finally:
        with _ACTIVE_LOCK:
            _ACTIVE_PROCESSES.discard(proc)

    completed = subprocess.CompletedProcess(cmd, proc.returncode, out, err)
    if check and proc.returncode:
        raise subprocess.CalledProcessError(proc.returncode, cmd, output=out, stderr=err)
    return completed
