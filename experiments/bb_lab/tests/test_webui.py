"""Tests for the browser UI's analysis, flag discovery and premise gating.

The distance itself is covered by the existing solver tests; what is new
here — and what would silently produce wrong physics if it broke — is the
layer that decides *which solver hints a code is allowed to be given*.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from bb_lab.webui import analysis, solver
from bb_lab.webui.analysis import CodeInputError

GROSS = ("12x6", "x^3 + y + y^2", "y^3 + x + x^2")
BB72 = ("6x6", "x^3 + y + y^2", "y^3 + x + x^2")
TORIC3 = ("3x3", "1 + x", "1 + y")


# ------------------------------------------------------------- group specs


@pytest.mark.parametrize(
    "spec, orders",
    [
        ("12x6", (12, 6)),
        ("12X6", (12, 6)),
        ("12 × 6", (12, 6)),
        ("12, 6", (12, 6)),
        ("[12, 6]", (12, 6)),
        ("6x6x3", (6, 6, 3)),
        ("9", (9,)),
        ([12, 6], (12, 6)),
    ],
)
def test_parse_orders_accepts_the_spellings_people_use(spec, orders):
    assert analysis.parse_orders(spec).orders == orders


@pytest.mark.parametrize(
    "spec", ["", "abc", "12x0", "12x-6", "3x3x3x3x3x3x3", "80x80"]
)
def test_parse_orders_rejects_nonsense_with_a_readable_message(spec):
    with pytest.raises(CodeInputError):
        analysis.parse_orders(spec)


def test_poly_errors_are_rephrased_not_leaked():
    G = analysis.parse_orders("12x6")
    with pytest.raises(CodeInputError, match="A:"):
        analysis.parse_poly("xy", G, "A")          # implicit product
    with pytest.raises(CodeInputError, match="empty"):
        analysis.parse_poly("   ", G, "B")
    with pytest.raises(CodeInputError, match="reduces to 0"):
        analysis.parse_poly("x + x", G, "A")       # cancels over F₂


# ----------------------------------------------------------------- report


@pytest.mark.parametrize(
    "spec, n, k, weight",
    [
        (BB72, 72, 12, 6),
        (GROSS, 144, 12, 6),
        (("15x3", "x^9 + y + y^2", "1 + x^2 + x^7"), 90, 8, 6),
        (("9x6", "x^3 + y + y^2", "y^3 + x + x^2"), 108, 8, 6),
        (TORIC3, 18, 2, 4),
    ],
)
def test_report_matches_the_published_parameters(spec, n, k, weight):
    report, _ = analysis.analyse(*spec, lookup_corpus=False)
    assert (report.n, report.k) == (n, k)
    assert report.check_weight == weight
    assert report.css_commutes
    assert not report.warnings


def test_check_weight_is_wtA_plus_wtB_and_degree_is_wtA():
    report, _ = analysis.analyse(*GROSS, lookup_corpus=False)
    assert report.check_weight == report.A_weight + report.B_weight
    # Each qubit sits in wt(A) X-checks (and as many Z-checks).
    assert report.qubit_degree == report.A_weight == report.B_weight
    assert report.num_checks == 2 * report.group_order
    assert report.n == 2 * report.group_order


def test_k_zero_is_reported_not_crashed():
    # A = 1 is a unit of F₂[G], so M_A = I and [M_A | M_B] has full rank
    # |G|: nothing is encoded, and the distance question is vacuous.
    report, _ = analysis.analyse("6x6", "1", "x", lookup_corpus=False)
    assert report.k == 0
    assert any("k = 0" in w for w in report.warnings)


def test_rank_three_group_is_accepted():
    report, _ = analysis.analyse(
        "3x3x2", "1 + x + z", "1 + y + z", lookup_corpus=False
    )
    assert report.n == 2 * 18
    assert report.css_commutes


# --------------------------------------------------------------- premises


def test_parity_premise_holds_for_bravyi_codes():
    for spec in (BB72, GROSS):
        report, _ = analysis.analyse(*spec, lookup_corpus=False)
        assert report.premises["coset_parity_even"]["holds"]


def test_parity_premise_fails_when_logicals_are_odd():
    # Toric code on an odd torus: loops have odd weight, so not every
    # feasible cost is even and -cost-step=2 would be unsound.
    report, _ = analysis.analyse(*TORIC3, lookup_corpus=False)
    assert not report.premises["coset_parity_even"]["holds"]


def test_premise_agrees_with_the_sweep_scripts_own_test():
    """The UI must compute the flag's precondition exactly as the batteries do."""
    from bb_lab.checks import bb_check_matrices
    from bb_lab.group import ZmZn
    from bb_lab.linalg import nullspace_f2, quotient_complement_basis
    from bb_lab.poly import Poly

    for orders, a, b in (BB72, GROSS, TORIC3):
        ell, m = (int(t) for t in orders.split("x"))
        G = ZmZn(ell, m)
        checks = bb_check_matrices(
            Poly.from_string(a, G), Poly.from_string(b, G)
        )
        # verbatim from scripts/tandem_verify.py
        V = quotient_complement_basis(checks.H_X, nullspace_f2(checks.H_Z))
        expected = (
            not any(int(r.sum()) % 2 for r in checks.H_X)
            and not any(int(v.sum()) % 2 for v in V)
        )
        assert analysis.premises(checks)["coset_parity_even"]["holds"] == expected


# ------------------------------------------------------- option discovery

HELP_SAMPLE = """\
USAGE: tandem [options] <input-file> <result-output-file>

CORE OPTIONS:

  -rnd-init, -no-rnd-init                 (default: off)

  -rinc         = <double> (   1 ..  inf) (default: 2)
  -rfirst       = <int32>  [   1 .. imax] (default: 100)

MAIN OPTIONS:

  -pre, -no-pre                           (default: on)

  -cost-step    = <int32>  [   1 ..   64] (default: 1)
  -init-lb      = <int32>  [   0 .. imax] (default: 0)

  -phase-file = <string>
  -prime-vars = <string>

HELP OPTIONS:

  --help        Print help message.
"""


def test_help_parser_reads_types_domains_and_defaults():
    opts = {o.flag: o for o in solver.parse_help(HELP_SAMPLE)}
    assert opts["-cost-step"].kind == "int32"
    assert opts["-cost-step"].domain == "1 .. 64"
    assert opts["-cost-step"].default == "1"
    assert opts["-pre"].kind == "bool"
    assert opts["-pre"].negation == "-no-pre"
    assert opts["-phase-file"].kind == "string"
    assert opts["-rinc"].domain == "1 .. inf"
    # HELP OPTIONS are not offered as solver flags.
    assert "--help" not in opts


def test_flag_payload_carries_what_the_page_renders():
    """The UI's compact rows key off these fields; keep them present."""
    opts = {o.flag: o.to_json() for o in solver.parse_help(HELP_SAMPLE)}
    step = opts["-cost-step"]
    assert step["label"] == "Cost step"          # short enough for one line
    assert step["suggested"] == 2                # prefilled value
    assert step["short"]                         # one-clause inline caveat
    assert step["blurb"] and len(step["blurb"]) > len(step["short"])  # hover text
    assert step["featured"] is True
    # -init-lb has no recommended value: the box stays empty until the user
    # supplies a floor they can certify.
    assert opts["-init-lb"]["suggested"] is None
    # An undecorated option still renders, just plainly.
    assert opts["-rfirst"]["short"] == ""
    assert opts["-rfirst"]["blurb"] is None


def test_unknown_flags_survive_discovery_without_code_changes():
    """A future Tandem flag must appear on its own — that is the whole point."""
    extended = HELP_SAMPLE.replace(
        "  -init-lb      = <int32>  [   0 .. imax] (default: 0)",
        "  -init-lb      = <int32>  [   0 .. imax] (default: 0)\n"
        "  -future-hint  = <int32>  [   0 ..   99] (default: 7)",
    )
    opts = {o.flag: o for o in solver.parse_help(extended)}
    assert opts["-future-hint"].default == "7"
    payload = opts["-future-hint"].to_json()
    # No FLAG_NOTES entry, so it is offered plainly and treated as safe.
    assert payload["soundness"] == "safe"
    assert payload["requires"] is None
    assert payload["featured"] is False


# ----------------------------------------------------------- flag gating


def _backend(help_text: str = HELP_SAMPLE) -> solver.BackendInfo:
    return solver.BackendInfo(
        available=True, path="/nonexistent/tandem",
        options=solver.parse_help(help_text), fork=True,
    )


HOLDS = {"coset_parity_even": {"holds": True, "label": "L", "detail": "D"}}
FAILS = {"coset_parity_even": {"holds": False, "label": "L", "detail": "D"}}


def test_flag_with_satisfied_premise_is_passed_through():
    argv = solver.build_argv(
        "tandem", Path("i.wcnf"), {"-cost-step": 2},
        backend=_backend(), premises=HOLDS,
    )
    assert argv == ["tandem", "-cost-step=2", "i.wcnf"]


def test_flag_with_violated_premise_is_refused():
    with pytest.raises(solver.FlagRejected, match="premise"):
        solver.build_argv(
            "tandem", Path("i.wcnf"), {"-cost-step": 2},
            backend=_backend(), premises=FAILS,
        )


def test_uncheckable_obligation_needs_explicit_acknowledgement():
    with pytest.raises(solver.FlagRejected, match="acknowledged"):
        solver.build_argv(
            "tandem", Path("i.wcnf"), {"-init-lb": 12},
            backend=_backend(), premises=HOLDS,
        )
    argv = solver.build_argv(
        "tandem", Path("i.wcnf"), {"-init-lb": 12},
        backend=_backend(), premises=HOLDS, acknowledged=["-init-lb"],
    )
    assert "-init-lb=12" in argv


def test_flags_the_binary_does_not_have_are_refused():
    with pytest.raises(solver.FlagRejected, match="not an option"):
        solver.build_argv(
            "tandem", Path("i.wcnf"), {"-cost-step": 2},
            backend=_backend(HELP_SAMPLE.replace("-cost-step", "-gone")),
            premises=HOLDS,
        )


def test_booleans_use_the_negated_form_when_switched_off():
    argv = solver.build_argv(
        "tandem", Path("i.wcnf"), {"-pre": True},
        backend=_backend(), premises=HOLDS,
    )
    assert "-pre" in argv
    # False/None/"" are treated as "leave it alone", not "negate it".
    argv = solver.build_argv(
        "tandem", Path("i.wcnf"), {"-pre": False, "-cost-step": None},
        backend=_backend(), premises=HOLDS,
    )
    assert argv == ["tandem", "i.wcnf"]


def test_method_string_matches_the_corpus_vocabulary():
    # scripts/merit_writeback.py: {1: 'maxsat-tandem@mse23+step2',
    #                              0: 'maxsat-tandem@mse23'}
    assert solver.method_string({"-cost-step": 2}, fork=True) == \
        "maxsat-tandem@mse23+step2"
    assert solver.method_string({}, fork=True) == "maxsat-tandem@mse23"
    assert solver.method_string({}, fork=False) == "maxsat-maxcdcl@mse23"
    # Resource caps are not provenance.
    assert solver.method_string({"-cpu-lim": 60}, fork=True) == \
        "maxsat-tandem@mse23"


# ----------------------------------------------------------------- probe


def test_probe_is_honest_about_a_missing_binary():
    info = solver.probe("/definitely/not/here/tandem")
    assert not info.available
    assert "does not exist" in info.error
    info = solver.probe(None)
    assert not info.available


@pytest.mark.skipif(
    not solver.DEFAULT_TANDEM.exists(),
    reason="Tandem not built (third_party/build_maxcdcl.sh)",
)
def test_probe_finds_the_fork_and_its_hint_flags():
    info = solver.probe(solver.DEFAULT_TANDEM)
    assert info.available and info.fork
    flags = {o.flag for o in info.options}
    assert {"-cost-step", "-init-lb", "-phase-file", "-prime-vars"} <= flags


# ------------------------------------------------------------ end to end


@pytest.fixture(scope="module")
def live_server():
    """The real server on an ephemeral port, no solver binary configured."""
    import threading

    from bb_lab.webui.server import serve

    httpd = serve(host="127.0.0.1", port=0, binary=None)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    host, port = httpd.server_address[:2]
    try:
        yield f"http://{host}:{port}"
    finally:
        httpd.shutdown()
        httpd.server_close()


def _get(base: str, path: str):
    import json
    import urllib.request

    with urllib.request.urlopen(base + path) as r:
        return json.load(r)


def _post(base: str, path: str, body: dict):
    import json
    import urllib.error
    import urllib.request

    req = urllib.request.Request(
        base + path, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req) as r:
            return r.status, json.load(r)
    except urllib.error.HTTPError as e:
        return e.code, json.load(e)


def test_index_and_static_assets_are_served(live_server):
    import urllib.request

    for path in ("/", "/static/app.js", "/static/style.css"):
        with urllib.request.urlopen(live_server + path) as r:
            assert r.status == 200
            assert r.read()


def test_static_route_refuses_path_traversal(live_server):
    import urllib.error
    import urllib.request

    with pytest.raises(urllib.error.HTTPError) as e:
        urllib.request.urlopen(live_server + "/static/../server.py")
    assert e.value.code == 404


def test_presets_carry_the_bravyi_table(live_server):
    names = {p["code_id"] for p in _get(live_server, "/api/presets")["presets"]}
    assert {"bb_72_12_6", "gross", "bb_288_12_18"} <= names


def test_analyse_endpoint_reports_gross(live_server):
    status, body = _post(live_server, "/api/analyse", {
        "orders": "12x6", "A": "x^3 + y + y^2", "B": "y^3 + x + x^2",
    })
    assert status == 200
    assert (body["n"], body["k"], body["check_weight"]) == (144, 12, 6)


def test_analyse_endpoint_returns_a_usable_error(live_server):
    status, body = _post(live_server, "/api/analyse",
                         {"orders": "nope", "A": "x", "B": "y"})
    assert status == 400
    assert "whole number" in body["error"]


def test_solve_endpoint_streams_to_a_verified_distance(live_server):
    """Full round trip on the fallback backend: no binary, real d."""
    import json
    import urllib.request

    status, job = _post(live_server, "/api/solve", {
        "orders": "6x6", "A": "x^3 + y + y^2", "B": "y^3 + x + x^2",
        "backend": "sat-ladder",
    })
    assert status == 200

    events = []
    with urllib.request.urlopen(
        f"{live_server}/api/solve/{job['job_id']}/events"
    ) as stream:
        for raw in stream:
            line = raw.decode().strip()
            if line.startswith("data: "):
                events.append(json.loads(line[6:]))
                if events[-1]["kind"] in ("done", "error", "cancelled"):
                    break

    assert [e["kind"] for e in events].count("rung") == 6   # w = 1 … 6
    done = events[-1]
    assert done["kind"] == "done"
    assert done["distance"] == 6            # published [[72,12,6]]
    assert done["witness_weight"] == 6
    assert done["verified"] is True


def test_solve_refuses_a_code_with_no_logicals(live_server):
    status, body = _post(live_server, "/api/solve",
                         {"orders": "6x6", "A": "1", "B": "x"})
    assert status == 400
    assert "k = 0" in body["error"]


def test_solve_refuses_tandem_flags_when_no_binary_is_configured(live_server):
    status, body = _post(live_server, "/api/solve", {
        "orders": "6x6", "A": "x^3 + y + y^2", "B": "y^3 + x + x^2",
        "backend": "tandem", "flags": {"-cost-step": 2},
    })
    assert status == 422
    assert "does not exist" in body["error"] or "configured" in body["error"]
