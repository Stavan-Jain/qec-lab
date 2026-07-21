"""Recursive deck-birth survey over the IBM qcode-discovery CSS catalog.

Parses the supplemental CSS tables of arXiv:2606.02418 (Tables I-III, all
verified codes at n = 144/288/360) from a pdftotext dump, validates every row
by recomputing k from the check matrices, then walks each distinct class's
x/y halving tower recording k of every order-2 deck quotient (x-half, y-half,
diagonal — general coset convolution, so non-product diagonal quotients are
measured too). The "birth rung" is the deepest frame reachable through
k-preserving x/y decks at which every deck jumps (or the odd core).

Usage (from experiments/bb_lab):
    pdftotext -layout -f 23 -l 42 <2606.02418.pdf> data/a19/ibm_supp.txt
    uv run scripts/a19_deck_survey_catalog.py data/a19/ibm_supp.txt \
        --json data/a19/ibm_catalog_css.json

Validation: reproduces all four A19 SS7 Bravyi-table verdicts (anchors run
first); rows whose recomputed k mismatches the table's k column are excluded
and reported.

Caveats (doctrine: A10-L1 / A19 SS7): literal presentations only — a JUMP is
decisive code-level only after an Aut+monomial orbit sweep; k-preservation
exhibited by one presentation is sufficient one-sided evidence. The descent
follows the first k-preserving x/y branch (halvings commute, but
path-dependence of k-preservation is not swept). Cells where only the
diagonal deck preserves k are classified "diag-survivor" (descent continues
through a non-product quotient, not followed).

2026-07-20 findings (research_log entry a19-wall-refutation-and-catalog-survey):
every k = 12 class has a unique deck-birth rung; [[288,12,16]] (class a) is a
second deck-nontrivial-top flagship; the (18,.)-family ([[288,8,20]] d = 20
exact, born (18,2)) is the cheapest known d >= 20 certification target.
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

BRAVYI_ANCHORS = [  # label, ell, m, A, B, expected k (A19 SS7 verdicts)
    ("BRAVYI [[72,12,6]]", 6, 6, "x^3+y+y^2", "y^3+x+x^2", 12),
    ("BRAVYI [[144,12,12]]", 12, 6, "x^3+y+y^2", "y^3+x+x^2", 12),
    ("BRAVYI [[288,12,18]]", 12, 12, "x^3+y^2+y^7", "y^3+x+x^2", 12),
    ("BRAVYI [[360,12,24]]", 30, 6, "x^9+y+y^2", "y^3+x^25+x^26", 12),
]


def parse_poly(s):
    """'y+y 2 +x3' (pdftotext superscripts) or 'x^3+y+y^2' -> F2 support set."""
    s = s.replace("^", "").replace(" ", "")
    supp = {}
    for term in s.split("+"):
        if not term:
            continue
        if term == "1":
            e = (0, 0)
        else:
            mo = re.fullmatch(r"(?:x(\d*))?(?:y(\d*))?", term)
            if mo is None or (mo.group(1) is None and mo.group(2) is None):
                return None
            ex = int(mo.group(1)) if mo.group(1) else (1 if "x" in term else 0)
            ey = int(mo.group(2)) if mo.group(2) else (1 if "y" in term else 0)
            e = (ex, ey)
        supp[e] = supp.get(e, 0) ^ 1
    return {e for e, c in supp.items() if c}


def rank_f2(M):
    M = M.copy().astype(np.uint8)
    r = 0
    for c in range(M.shape[1]):
        piv = np.nonzero(M[r:, c])[0]
        if piv.size == 0:
            continue
        p = r + piv[0]
        M[[r, p]] = M[[p, r]]
        rows = np.nonzero(M[:, c])[0]
        M[rows[rows != r]] ^= M[r]
        r += 1
        if r == M.shape[0]:
            break
    return r


def _hxhz(MA, MB):
    HX = np.concatenate([MA, MB], axis=1)
    HZ = np.concatenate([MB.T, MA.T], axis=1)
    return HX, HZ


def code_k(ell, m, sA, sB):
    elems = [(a, b) for a in range(ell) for b in range(m)]
    idx = {e: i for i, e in enumerate(elems)}

    def conv(supp):
        M = np.zeros((len(elems),) * 2, dtype=np.uint8)
        for j, (gx, gy) in enumerate(elems):
            for (ax, ay) in supp:
                M[idx[((gx + ax) % ell, (gy + ay) % m)], j] ^= 1
        return M

    HX, HZ = _hxhz(conv(sA), conv(sB))
    return 2 * ell * m - rank_f2(HX) - rank_f2(HZ)


def quotient_k(ell, m, sA, sB, t):
    """k of the quotient by the order-2 subgroup <t> (coset convolution)."""
    tx, ty = t
    seen, reps, ridx = set(), [], {}
    for a in range(ell):
        for b in range(m):
            c = frozenset({(a, b), ((a + tx) % ell, (b + ty) % m)})
            if c not in seen:
                seen.add(c)
                reps.append(min(c))
    for i, r in enumerate(reps):
        ridx[r] = i
        ridx[((r[0] + tx) % ell, (r[1] + ty) % m)] = i

    def qconv(supp):
        M = np.zeros((len(reps),) * 2, dtype=np.uint8)
        for j, (gx, gy) in enumerate(reps):
            for (ax, ay) in supp:
                M[ridx[((gx + ax) % ell, (gy + ay) % m)], j] ^= 1
        return M

    HX, HZ = _hxhz(qconv(sA), qconv(sB))
    return 2 * len(reps) - rank_f2(HX) - rank_f2(HZ)


def half_poly(supp, axis, half):
    out = {}
    for (ax, ay) in supp:
        e = (ax % half, ay) if axis == "x" else (ax, ay % half)
        out[e] = out.get(e, 0) ^ 1
    return {e for e, c in out.items() if c}


def descend(ell, m, sA, sB, k, depth=0):
    lines, decks = [], []
    if ell % 2 == 0:
        decks.append(("x", (ell // 2, 0)))
    if m % 2 == 0:
        decks.append(("y", (0, m // 2)))
    if ell % 2 == 0 and m % 2 == 0:
        decks.append(("xy", (ell // 2, m // 2)))
    if not decks:
        lines.append(f"{'  ' * depth}({ell},{m}) k={k}: odd core — k born at base")
        return (f"({ell},{m}) odd-core", lines)
    kq = {nm: quotient_k(ell, m, sA, sB, t) for nm, t in decks}
    lines.append("  " * depth + f"({ell},{m}) k={k}: " + ", ".join(
        f"{nm}->{kv}{'' if kv == k else ' JUMP(lit)'}" for nm, kv in kq.items()))
    for nm, _ in decks:
        if nm in ("x", "y") and kq[nm] == k:
            half = ell // 2 if nm == "x" else m // 2
            nell, nm_ = (half, m) if nm == "x" else (ell, half)
            birth, sub = descend(nell, nm_, half_poly(sA, nm, half),
                                 half_poly(sB, nm, half), k, depth + 1)
            lines.extend(sub)
            return (birth, lines)
    if any(v == k for v in kq.values()):
        lines.append("  " * depth + f"  => k={k} survives only the DIAGONAL "
                     f"deck at ({ell},{m}) — non-product descent not followed")
        return (f"({ell},{m}) diag-survivor", lines)
    lines.append("  " * depth +
                 f"  => k={k} BORN at ({ell},{m}) — all decks jump (lit)")
    return (f"({ell},{m}) deck-birth", lines)


def parse_tables(path):
    """Header-detected table sections; form-feed safe (\\f -> \\n)."""
    text = Path(path).read_text().replace("\f", "\n")
    section, rows = None, []
    for line in text.splitlines():
        mo = re.search(r"TABLE\s+\S+\s+All verified codes at n = (\d+)", line)
        if mo:
            n = int(mo.group(1))
            section = f"n{n}" if n in (144, 288, 360) else None  # CSS only
            continue
        if section is None:
            continue
        f = re.split(r"\s{2,}", line.strip())
        if len(f) < 6:
            continue
        fr = re.fullmatch(r"\((\d+),\s*(\d+)\)", f[1])
        if not fr:
            continue
        try:
            k_tab, d_tab = int(f[4]), int(f[5])
        except ValueError:
            continue
        ell, m = int(fr.group(1)), int(fr.group(2))
        sA, sB = parse_poly(f[2]), parse_poly(f[3])
        if sA is None or sB is None:
            print(f"  PARSE-FAIL {section} {f[0]}: A={f[2]!r} B={f[3]!r}")
            continue
        if 2 * ell * m != int(section[1:]):
            print(f"  FRAME/TABLE MISMATCH {section} {f[0]} ({ell},{m}) — "
                  f"row kept, reattributed to n{2*ell*m}")
            section_row = f"n{2*ell*m}"
        else:
            section_row = section
        rows.append({"table": section_row, "cls": f[0], "ell": ell, "m": m,
                     "A": f[2], "B": f[3], "sA": sA, "sB": sB,
                     "k": k_tab, "d": d_tab})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("supp_txt", help="pdftotext -layout dump of the "
                    "arXiv:2606.02418 supplemental catalog pages")
    ap.add_argument("--json", default=None,
                    help="write validated rows + birth verdicts as JSON")
    args = ap.parse_args()

    for label, ell, m, A, B, kexp in BRAVYI_ANCHORS:
        sA, sB = parse_poly(A), parse_poly(B)
        k = code_k(ell, m, sA, sB)
        birth, lines = descend(ell, m, sA, sB, k)
        status = "OK" if k == kexp else "ANCHOR-FAIL"
        print(f"== {label} k={k} (expect {kexp}) [{status}]  birth: {birth}")
        print("\n".join(lines))
        if k != kexp:
            sys.exit("anchor validation failed — aborting")

    rows = parse_tables(args.supp_txt)
    ok, bad = [], []
    for r in rows:
        kc = code_k(r["ell"], r["m"], r["sA"], r["sB"])
        (ok if kc == r["k"] else bad).append((r, kc))
    print(f"\nparsed {len(rows)} rows: {len(ok)} k-validated, "
          f"{len(bad)} mismatched (excluded)")
    for r, kc in bad:
        print(f"  MISMATCH {r['table']}/{r['cls']} ({r['ell']},{r['m']}) "
              f"table k={r['k']} computed {kc}  A={r['A']!r} B={r['B']!r}")

    out, seen = [], set()
    for r, _ in ok:
        key = (r["table"], r["cls"])
        if key in seen:
            continue
        seen.add(key)
        birth, lines = descend(r["ell"], r["m"], r["sA"], r["sB"], r["k"])
        print(f"== {r['table']}/{r['cls']} [[{2*r['ell']*r['m']},{r['k']},"
              f"{r['d']}]]  A={r['A']}  B={r['B']}")
        print("\n".join(lines))
        out.append({k2: v for k2, v in r.items() if k2 not in ("sA", "sB")}
                   | {"birth": birth})

    print("\n===== SUMMARY (distinct classes, first representation) =====")
    for r in out:
        print(f"{r['table']}/{r['cls']:<4s} k={r['k']:<3d} d={r['d']:<3d} "
              f"birth: {r['birth']}")
    if args.json:
        Path(args.json).write_text(json.dumps(out, indent=1))
        print(f"\nwrote {len(out)} classes -> {args.json}")


if __name__ == "__main__":
    main()
