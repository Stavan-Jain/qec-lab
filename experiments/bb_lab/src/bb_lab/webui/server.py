"""The local HTTP server behind `bb-lab ui`.

Deliberately stdlib-only (`http.server` + threads): the lab's dependency
set is a solver stack, and a single-user localhost tool is not a reason
to add a web framework to it. Requests are JSON; solver progress is a
Server-Sent Events stream, which is the one browser primitive that gives
live output with no client library at all.
"""

from __future__ import annotations

import json
import threading
import traceback
from dataclasses import asdict
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

import yaml

from . import analysis, solver
from .analysis import CodeInputError
from .solver import FlagRejected, Job

STATIC = Path(__file__).resolve().parent / "static"
PRESETS_YAML = analysis.LAB_ROOT / "instances" / "bravyi_table.yaml"

_JOBS: dict[str, Job] = {}
_JOBS_LOCK = threading.Lock()

_CONTENT_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".svg": "image/svg+xml",
}


class Config:
    """Server-wide settings, set once by `serve()`."""

    binary: str | None = None
    backend: solver.BackendInfo | None = None

    @classmethod
    def refresh(cls) -> solver.BackendInfo:
        cls.backend = solver.probe(cls.binary)
        return cls.backend


def _presets() -> list[dict[str, Any]]:
    try:
        rows = yaml.safe_load(PRESETS_YAML.read_text())["instances"]
    except Exception:
        return []
    out = []
    for r in rows:
        g = r["group"]
        out.append({
            "code_id": r["code_id"],
            "name": r["display_name"],
            "orders": f"{g['ell']}x{g['m']}",
            "A": r["polynomials"]["A"],
            "B": r["polynomials"]["B"],
            "d": r["parameters"]["d"],
        })
    return out


class Handler(BaseHTTPRequestHandler):
    server_version = "bb-lab-ui"

    # Quiet by default: one line per request is noise for a UI that
    # polls, and the interesting logging is the solver's.
    def log_message(self, fmt: str, *args: Any) -> None:
        pass

    # ---------------------------------------------------------- plumbing

    def _send_json(self, obj: Any, status: int = 200) -> None:
        body = json.dumps(obj).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_error_json(self, message: str, status: int = 400, **extra: Any) -> None:
        self._send_json({"error": message, **extra}, status=status)

    def _body(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length") or 0)
        if not length:
            return {}
        return json.loads(self.rfile.read(length) or b"{}")

    def _send_static(self, name: str) -> None:
        path = (STATIC / name).resolve()
        if not path.is_file() or STATIC.resolve() not in path.parents:
            self.send_error(404)
            return
        body = path.read_bytes()
        self.send_response(200)
        self.send_header(
            "Content-Type",
            _CONTENT_TYPES.get(path.suffix, "application/octet-stream"),
        )
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    # ------------------------------------------------------------- routes

    def do_GET(self) -> None:  # noqa: N802  (stdlib naming)
        path = self.path.split("?", 1)[0]
        try:
            if path in ("/", "/index.html"):
                self._send_static("index.html")
            elif path.startswith("/static/"):
                self._send_static(path[len("/static/"):])
            elif path == "/api/backend":
                self._send_json(self._backend_payload())
            elif path == "/api/presets":
                self._send_json({"presets": _presets()})
            elif path.startswith("/api/solve/") and path.endswith("/events"):
                self._stream_events(path.split("/")[3])
            else:
                self.send_error(404)
        except BrokenPipeError:
            pass
        except Exception:
            traceback.print_exc()
            try:
                self._send_error_json("internal error", 500)
            except Exception:
                pass

    def do_POST(self) -> None:  # noqa: N802
        path = self.path.split("?", 1)[0]
        try:
            if path == "/api/analyse":
                self._analyse()
            elif path == "/api/solve":
                self._solve()
            elif path.startswith("/api/solve/") and path.endswith("/cancel"):
                self._cancel(path.split("/")[3])
            elif path == "/api/backend/refresh":
                Config.refresh()
                self._send_json(self._backend_payload())
            else:
                self.send_error(404)
        except CodeInputError as e:
            self._send_error_json(str(e))
        except FlagRejected as e:
            self._send_error_json(str(e), 422)
        except json.JSONDecodeError:
            self._send_error_json("malformed JSON body")
        except BrokenPipeError:
            pass
        except Exception as e:
            traceback.print_exc()
            self._send_error_json(f"{type(e).__name__}: {e}", 500)

    # --------------------------------------------------------- handlers

    def _backend_payload(self) -> dict[str, Any]:
        backend = Config.backend or Config.refresh()
        return {
            "tandem": backend.to_json(),
            "fallback": {
                "name": "CryptoMiniSat ladder",
                "available": True,
                "note": (
                    "In-process, no binary needed. Slower than Tandem and "
                    "only interruptible between weight rungs."
                ),
            },
        }

    def _analyse(self) -> None:
        body = self._body()
        report, _ = analysis.analyse(
            body.get("orders", ""), body.get("A", ""), body.get("B", ""),
        )
        self._send_json(asdict(report))

    def _solve(self) -> None:
        body = self._body()
        report, checks = analysis.analyse(
            body.get("orders", ""), body.get("A", ""), body.get("B", ""),
            lookup_corpus=False,
        )
        if report.k == 0:
            self._send_error_json(
                "This code has k = 0 — there are no logical qubits, so "
                "there is no distance to compute."
            )
            return

        backend_choice = body.get("backend", "tandem")
        flags = {k: v for k, v in (body.get("flags") or {}).items()}
        ack = body.get("acknowledged") or []
        mode = body.get("mode", "naive")
        if mode not in ("naive", "strengthened"):
            self._send_error_json(f"unknown encoding mode {mode!r}")
            return

        backend = Config.backend or Config.refresh()
        use_tandem = backend_choice == "tandem" and backend.available
        if backend_choice == "tandem" and not backend.available:
            self._send_error_json(
                backend.error or "solver binary unavailable", 422,
            )
            return

        if use_tandem:
            # Validate the flags before the job exists, so a rejection is
            # a clean 422 rather than an error event on a dead stream.
            solver.build_argv(
                backend.path, Path("/dev/null"), flags,
                backend=backend, premises=report.premises, acknowledged=ack,
            )

        label = f"[[{report.n},{report.k}]] {report.group_label}"
        job = Job(label, "tandem" if use_tandem else "sat-ladder")
        with _JOBS_LOCK:
            _JOBS[job.id] = job

        def work() -> None:
            try:
                if use_tandem:
                    solver.run_tandem(
                        job, checks, backend.path, flags,
                        backend=backend, premises=report.premises,
                        acknowledged=ack, mode=mode,
                    )
                else:
                    solver.run_sat_ladder(job, checks)
            except solver.Cancelled:
                job.state = "cancelled"
                job.emit("cancelled")
            except Exception as e:
                job.fail(f"{type(e).__name__}: {e}")

        threading.Thread(target=work, daemon=True).start()
        self._send_json({
            "job_id": job.id,
            "backend": job.backend,
            "label": label,
            "n": report.n,
            "k": report.k,
        })

    def _cancel(self, job_id: str) -> None:
        with _JOBS_LOCK:
            job = _JOBS.get(job_id)
        if job is None:
            self._send_error_json("no such job", 404)
            return
        job.cancel()
        self._send_json({"ok": True, "state": job.state})

    def _stream_events(self, job_id: str) -> None:
        with _JOBS_LOCK:
            job = _JOBS.get(job_id)
        if job is None:
            self.send_error(404)
            return
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Connection", "keep-alive")
        self.end_headers()
        seq = -1
        try:
            while True:
                events = job.since(seq, timeout=15.0)
                for e in events:
                    seq = max(seq, e.seq)
                    payload = json.dumps({"kind": e.kind, **e.payload})
                    self.wfile.write(
                        f"event: {e.kind}\ndata: {payload}\n\n".encode()
                    )
                self.wfile.flush()
                if not events:
                    self.wfile.write(b": keep-alive\n\n")   # hold the socket
                    self.wfile.flush()
                if job.finished and not job.since(seq, timeout=0.0):
                    break
        except (BrokenPipeError, ConnectionResetError):
            pass


def serve(
    *,
    host: str = "127.0.0.1",
    port: int = 8765,
    binary: str | Path | None = None,
) -> ThreadingHTTPServer:
    Config.binary = str(binary) if binary else None
    Config.refresh()
    httpd = ThreadingHTTPServer((host, port), Handler)
    httpd.daemon_threads = True
    return httpd
