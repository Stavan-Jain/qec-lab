"""Distance backends for the web UI, built to outlive Tandem's next version.

Tandem is a moving target: it is a fork of MaxCDCL that the lab keeps
adding caller-verified hints to. So this module does not hardcode a
flag list. It **asks the binary** what it supports by parsing
`tandem --help` (a MiniSat-style option dump, complete with types,
ranges and defaults) and renders whatever it finds. A new flag appears
in the UI the moment you rebuild the solver.

What *is* hardcoded is the part a machine cannot infer: which flags
carry a soundness obligation. `FLAG_NOTES` below attaches a
human-readable blurb and, where one exists, the name of a premise from
`analysis.premises()` that must hold before the flag may be passed.
The check runs server-side on every request, so a stale browser tab or
a hand-rolled POST cannot talk the UI into an unsound run.

Adding support for a new Tandem flag is therefore either:
  * nothing at all — it shows up under "other solver options"; or
  * one `FLAG_NOTES` entry, if it needs a label or a premise.
"""

from __future__ import annotations

import itertools
import re
import subprocess
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np

from ..checks import CheckMatrices
from ..maxsat_distance import (
    decode_witness,
    parse_solver_output,
    verify_witness,
    write_wcnf,
)
from ..sat_distance import x_distance

LAB_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TANDEM = (
    LAB_ROOT / "third_party" / "maxcdcl" / "MaxCDCL" / "code" / "simp" / "tandem"
)

# Sections of `--help` that describe things the UI should not offer:
# HELP prints and exits; the positional proof/dimacs outputs are ours.
_SKIP_SECTIONS = {"HELP OPTIONS"}
_SKIP_FLAGS = {"-drup-file", "-dimacs", "-drup"}


# ------------------------------------------------------- option discovery

_SECTION_RE = re.compile(r"^([A-Z][A-Z /-]*OPTIONS):\s*$")
_BOOL_RE = re.compile(
    r"^\s*(-[\w-]+),\s*(-no-[\w-]+)\s*\(default:\s*(on|off)\)"
)
_VALUE_RE = re.compile(
    r"^\s*(-[\w-]+)\s*=\s*<(\w+)>"
    r"(?:\s*[\[(]\s*(.*?)\s*[\])])?"
    r"(?:\s*\(default:\s*(.*?)\))?\s*$"
)


@dataclass(frozen=True, slots=True)
class SolverOption:
    """One option as the binary itself describes it."""

    flag: str
    kind: str                    # 'bool' | 'int32' | 'double' | 'string'
    section: str
    default: str | None = None
    domain: str | None = None    # '1 .. 64', '0 .. inf', …
    negation: str | None = None  # '-no-pre' for booleans

    def to_json(self) -> dict[str, Any]:
        note = FLAG_NOTES.get(self.flag)
        return {
            "flag": self.flag,
            "kind": self.kind,
            "section": self.section,
            "default": self.default,
            "domain": self.domain,
            "label": note.label if note else self.flag.lstrip("-"),
            "blurb": note.blurb if note else None,
            "requires": note.requires if note else None,
            "soundness": note.soundness if note else "safe",
            "suggested": note.suggested if note else None,
            "featured": bool(note and note.featured),
        }


@dataclass(frozen=True, slots=True)
class FlagNote:
    """The human/soundness overlay on a discovered flag.

    `requires` names a key of `analysis.premises()`. When set, the
    server refuses the flag unless that premise holds for the code being
    solved. When `soundness == 'caller'` and `requires` is None, the
    obligation cannot be machine-checked and the UI demands an explicit
    acknowledgement instead.
    """

    label: str
    blurb: str
    soundness: str = "safe"      # 'safe' | 'bias' | 'caller'
    requires: str | None = None
    suggested: Any = None
    featured: bool = False


FLAG_NOTES: dict[str, FlagNote] = {
    "-cost-step": FlagNote(
        label="Cost step (weight parity)",
        blurb=(
            "Declares every feasible cost congruent mod N, letting the "
            "solver tighten its bound by N−1 after each improving model. "
            "N = 2 is exactly the coset weight-parity theorem, and is "
            "worth ~2–3× on BB instances."
        ),
        soundness="caller",
        requires="coset_parity_even",
        suggested=2,
        featured=True,
    ),
    "-init-lb": FlagNote(
        label="Certified floor",
        blurb=(
            "A lower bound on the optimum. The search stops the moment an "
            "incumbent reaches it, deleting the proof phase entirely when "
            "the floor is tight. Only pass analytically or kernel-certified "
            "floors — a wrong value silently returns a wrong distance."
        ),
        soundness="caller",
        featured=True,
    ),
    "-phase-file": FlagNote(
        label="Initial phases",
        blurb=(
            "Path to a 0/1 string over DIMACS vars used as initial branching "
            "phases — typically a known witness. Pure bias; cannot affect "
            "correctness."
        ),
        soundness="bias",
    ),
    "-prime-vars": FlagNote(
        label="Prime variables",
        blurb=(
            "Comma-separated 1-based DIMACS vars to branch on first with "
            "positive phase (the translation-symmetry anchor qubits). "
            "Pure bias; cannot affect correctness."
        ),
        soundness="bias",
    ),
    "-cpu-lim": FlagNote(
        label="CPU limit (s)",
        blurb="Solver-side CPU cap. Gives up rather than answering.",
        featured=True,
    ),
    "-mem-lim": FlagNote(
        label="Memory limit (MB)",
        blurb="Solver-side memory cap.",
    ),
    "-verb": FlagNote(
        label="Verbosity",
        blurb="0 quiet, 1 progress lines, 2 full statistics.",
    ),
}


def parse_help(text: str) -> list[SolverOption]:
    """Turn a MiniSat-family `--help` dump into structured options."""
    out: list[SolverOption] = []
    section = "OPTIONS"
    for line in text.splitlines():
        sec = _SECTION_RE.match(line.strip())
        if sec:
            section = sec.group(1)
            continue
        if section in _SKIP_SECTIONS:
            continue
        m = _BOOL_RE.match(line)
        if m:
            flag, neg, dflt = m.groups()
            if flag not in _SKIP_FLAGS:
                out.append(SolverOption(
                    flag=flag, kind="bool", section=section,
                    default=dflt, negation=neg,
                ))
            continue
        m = _VALUE_RE.match(line)
        if m:
            flag, kind, domain, dflt = m.groups()
            if flag not in _SKIP_FLAGS:
                out.append(SolverOption(
                    flag=flag, kind=kind, section=section,
                    default=dflt,
                    domain=re.sub(r"\s+", " ", domain) if domain else None,
                ))
    return out


@dataclass(frozen=True, slots=True)
class BackendInfo:
    """What the UI needs to know about the configured solver."""

    available: bool
    path: str | None
    options: list[SolverOption] = field(default_factory=list)
    fork: bool = False           # patched Tandem, or stock MaxCDCL?
    error: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "available": self.available,
            "path": self.path,
            "fork": self.fork,
            "error": self.error,
            "options": [o.to_json() for o in self.options],
        }


def probe(binary: str | Path | None) -> BackendInfo:
    """Ask the binary what it can do. Never raises."""
    if binary is None:
        return BackendInfo(
            available=False, path=None,
            error="No solver binary configured.",
        )
    path = Path(binary).expanduser()
    if not path.exists():
        return BackendInfo(
            available=False, path=str(path),
            error=(
                f"{path} does not exist. Build it with "
                "third_party/build_maxcdcl.sh, or pass --binary."
            ),
        )
    try:
        proc = subprocess.run(
            [str(path), "--help"], capture_output=True, text=True, timeout=30,
        )
    except Exception as e:  # unreadable, wrong arch, not executable …
        return BackendInfo(available=False, path=str(path), error=str(e))
    options = parse_help(proc.stdout + proc.stderr)
    if not options:
        return BackendInfo(
            available=False, path=str(path),
            error="Binary produced no recognisable --help output.",
        )
    fork = any(o.flag == "-cost-step" for o in options)
    return BackendInfo(
        available=True, path=str(path), options=options, fork=fork,
        error=None if fork else (
            "This binary has no -cost-step, so it is stock MaxCDCL rather "
            "than the Tandem fork. It will still solve; it just cannot "
            "take the lab's hints."
        ),
    )


# ------------------------------------------------------------ argv assembly


class FlagRejected(ValueError):
    """A requested flag failed its soundness precondition."""


def build_argv(
    binary: str | Path,
    wcnf: Path,
    flags: dict[str, Any],
    *,
    backend: BackendInfo,
    premises: dict[str, Any],
    acknowledged: Iterable[str] = (),
) -> list[str]:
    """Assemble the solver command line, enforcing every flag's premise.

    This is the choke point: the browser proposes, `build_argv` disposes.
    A flag whose `FlagNote.requires` premise is false is refused outright;
    a flag with an unverifiable caller obligation needs its name in
    `acknowledged`.
    """
    known = {o.flag: o for o in backend.options}
    ack = set(acknowledged)
    argv = [str(binary)]
    for flag, value in flags.items():
        if value in (None, "", False):
            continue
        if flag not in known:
            raise FlagRejected(
                f"{flag} is not an option of this solver build."
            )
        note = FLAG_NOTES.get(flag)
        if note is not None:
            if note.requires is not None:
                premise = premises.get(note.requires, {})
                if not premise.get("holds"):
                    raise FlagRejected(
                        f"{flag} needs the premise “{premise.get('label', note.requires)}”, "
                        f"which does not hold for this code "
                        f"({premise.get('detail', 'not checked')})."
                    )
            elif note.soundness == "caller" and flag not in ack:
                raise FlagRejected(
                    f"{flag} carries a caller obligation that cannot be "
                    "machine-checked; it must be acknowledged explicitly."
                )
        opt = known[flag]
        if opt.kind == "bool":
            argv.append(flag if value is True else str(opt.negation))
        else:
            argv.append(f"{flag}={value}")
    argv.append(str(wcnf))
    return argv


def method_string(flags: dict[str, Any], *, fork: bool) -> str:
    """The `d_method` tag this run would earn in the corpus.

    Matches the sweep scripts' vocabulary so a distance found here can be
    written back with the same provenance string it would have had from
    `ladder_sweep.py`.
    """
    base = "maxsat-tandem@mse23" if fork else "maxsat-maxcdcl@mse23"
    tags = []
    if str(flags.get("-cost-step", "")) == "2":
        tags.append("step2")
    for flag, value in sorted(flags.items()):
        if flag == "-cost-step" or value in (None, "", False):
            continue
        if flag in ("-cpu-lim", "-mem-lim", "-verb"):
            continue
        tags.append(f"{flag.lstrip('-')}={value}")
    return base + ("+" + "+".join(tags) if tags else "")


# ------------------------------------------------------------------- jobs


class Cancelled(Exception):
    """Raised inside a worker to unwind a cancelled run."""


@dataclass
class Event:
    seq: int
    kind: str
    payload: dict[str, Any]


class Job:
    """A running distance computation, streamable and cancellable."""

    _counter = itertools.count(1)

    def __init__(self, label: str, backend_name: str) -> None:
        self.id = f"job{next(self._counter)}"
        self.label = label
        self.backend = backend_name
        self.state = "running"          # running | done | error | cancelled
        self.events: list[Event] = []
        self.result: dict[str, Any] | None = None
        self.error: str | None = None
        self.started = time.perf_counter()
        self._seq = itertools.count()
        self._cond = threading.Condition()
        self._cancel = threading.Event()
        self._proc: subprocess.Popen | None = None
        self._log_budget = 4000

    # -- production side ------------------------------------------------

    def emit(self, kind: str, **payload: Any) -> None:
        payload.setdefault("elapsed", round(time.perf_counter() - self.started, 2))
        with self._cond:
            self.events.append(Event(next(self._seq), kind, payload))
            self._cond.notify_all()

    def log(self, line: str) -> None:
        if self._log_budget <= 0:
            return
        self._log_budget -= 1
        self.emit("log", line=line)

    def finish(self, result: dict[str, Any]) -> None:
        self.result = result
        self.state = "done"
        self.emit("done", **result)

    def fail(self, message: str) -> None:
        self.error = message
        if self.state != "cancelled":
            self.state = "error"
        self.emit("error", message=message)

    # -- control side ---------------------------------------------------

    @property
    def cancelled(self) -> bool:
        return self._cancel.is_set()

    def check_cancel(self) -> None:
        if self._cancel.is_set():
            raise Cancelled()

    def cancel(self) -> None:
        self._cancel.set()
        self.state = "cancelled"
        proc = self._proc
        if proc is not None and proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                proc.kill()
        self.emit("cancelled")

    # -- consumption side -----------------------------------------------

    def since(self, seq: int, timeout: float = 20.0) -> list[Event]:
        """Block for new events past `seq`; empty list on timeout."""
        with self._cond:
            if not any(e.seq > seq for e in self.events):
                self._cond.wait(timeout)
            return [e for e in self.events if e.seq > seq]

    @property
    def finished(self) -> bool:
        return self.state in ("done", "error", "cancelled")


# ------------------------------------------------------------ tandem runner


def _pump(stream, sink: Callable[[str], None]) -> None:
    try:
        for line in stream:
            sink(line.rstrip("\n"))
    except Exception:
        pass
    finally:
        try:
            stream.close()
        except Exception:
            pass


def run_tandem(
    job: Job,
    checks: CheckMatrices,
    binary: str | Path,
    argv_flags: dict[str, Any],
    *,
    backend: BackendInfo,
    premises: dict[str, Any],
    acknowledged: Iterable[str] = (),
    mode: str = "naive",
) -> None:
    """Solve with Tandem, streaming incumbents into `job`.

    Mirrors `maxsat_distance` — same encoder, same parser, same witness
    verification — but supervises the process itself so the incumbent
    (`o <cost>`) lines can be forwarded live instead of arriving in one
    lump at the end.
    """
    with tempfile.TemporaryDirectory(prefix="bb-ui-") as td:
        work = Path(td)
        job.emit("stage", stage="encoding",
                 detail=f"building the {mode} WCNF")
        qv = write_wcnf(checks, work / "instance.wcnf", mode=mode)
        job.check_cancel()

        argv = build_argv(
            binary, work / "instance.wcnf", argv_flags,
            backend=backend, premises=premises, acknowledged=acknowledged,
        )
        job.emit("stage", stage="solving", detail=" ".join(argv[1:-1]) or "(no flags)")

        out_lines: list[str] = []
        err_lines: list[str] = []

        def on_out(line: str) -> None:
            out_lines.append(line)
            if line.startswith("o "):
                try:
                    job.emit("incumbent", cost=int(line[2:].strip()))
                except ValueError:
                    job.log(line)
            elif line.startswith("v "):
                pass                      # the model; far too wide to log
            elif line.strip():
                job.log(line)

        t0 = time.perf_counter()
        proc = subprocess.Popen(
            argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, bufsize=1,
        )
        job._proc = proc
        # Cancel could have landed between build_argv and Popen, when there
        # was no process for it to kill; catch that here rather than let the
        # solver run to completion on a job nobody is waiting for.
        if job.cancelled and proc.poll() is None:
            proc.terminate()
        threads = [
            threading.Thread(target=_pump, args=(proc.stdout, on_out), daemon=True),
            threading.Thread(target=_pump, args=(proc.stderr, err_lines.append), daemon=True),
        ]
        for t in threads:
            t.start()
        proc.wait()
        for t in threads:
            t.join(timeout=5)
        dt = time.perf_counter() - t0
        job.check_cancel()

        stdout = "\n".join(out_lines)
        optimum, cost, bits = parse_solver_output(stdout)
        if not optimum or cost is None:
            raise RuntimeError(
                f"solver did not report OPTIMUM (exit {proc.returncode}); "
                f"tail:\n{stdout[-1200:]}\n{''.join(err_lines)[-400:]}"
            )
        if len(bits) < max(qv):
            raise RuntimeError(
                f"solver printed OPTIMUM but no usable v-line (cost={cost})"
            )

        job.emit("stage", stage="verifying",
                 detail="re-checking the witness against H_Z and the logicals")
        v = decode_witness(bits, qv)
        verify_witness(checks, v, cost)

        job.finish({
            "distance": int(cost),
            "witness_weight": int(v.sum()),
            "solver_seconds": round(dt, 2),
            "verified": True,
            "mode": mode,
            "method": method_string(argv_flags, fork=backend.fork),
            "command": " ".join(argv[1:-1]),
            "support": [int(i) for i in np.flatnonzero(v)],
        })


# --------------------------------------------------------- fallback runner


def run_sat_ladder(job: Job, checks: CheckMatrices) -> None:
    """CryptoMiniSat fallback: climb w = 1, 2, … until a logical appears.

    Slower than Tandem and it cannot be interrupted mid-rung, but it
    needs no external binary — so the UI still answers when the solver
    has not been built.
    """
    job.emit("stage", stage="solving",
             detail="iterated SAT (CryptoMiniSat), one rung per weight")

    def progress(w: int, sat: bool, seconds: float) -> None:
        job.emit("rung", weight=w, sat=sat, seconds=round(seconds, 2))
        job.check_cancel()

    t0 = time.perf_counter()
    res = x_distance(checks, progress=progress)
    dt = time.perf_counter() - t0

    job.emit("stage", stage="verifying", detail="re-checking the witness")
    verify_witness(checks, res.witness, int(res.witness.sum()))
    job.finish({
        "distance": int(res.distance),
        "witness_weight": int(res.witness.sum()),
        "solver_seconds": round(dt, 2),
        "verified": True,
        "mode": "sat-ladder",
        "method": "sat-cms-ladder",
        "command": "in-process CryptoMiniSat",
        "support": [int(i) for i in np.flatnonzero(res.witness)],
    })
