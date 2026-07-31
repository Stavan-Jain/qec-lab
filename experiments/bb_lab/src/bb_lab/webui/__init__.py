"""Browser front end for the BB lab: type a code, get [[n, k, d]].

`bb-lab ui` starts a localhost server (see `server.serve`). The split is:

  analysis.py  — the instant half of a code's report (n, k, weights,
                 CSS guard, solver premises, corpus lookup)
  solver.py    — the slow half: Tandem, discovered from its own --help
                 so the UI tracks the fork as it evolves
  server.py    — stdlib HTTP + Server-Sent Events
  static/      — the page itself
"""

from __future__ import annotations

from .solver import DEFAULT_TANDEM, probe
from .server import serve

__all__ = ["serve", "probe", "DEFAULT_TANDEM"]
