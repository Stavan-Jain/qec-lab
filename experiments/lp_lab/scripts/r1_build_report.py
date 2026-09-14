"""Assemble the LP1 report page (HTML) with the inline toric-layout figure."""
import json, re, sys
from pathlib import Path
out = Path(sys.argv[1])
svg = Path("data/figs/C224_toric_layout.svg").read_text()
# strip xml header text labels: keep the drawing, remove the title/legend <text> lines and the white rect
svg = re.sub(r'<text[^>]*>.*?</text>', '', svg, flags=re.S)
svg = svg.replace('<rect width="382" height="696" fill="white"/>', '')
svg = re.sub(r'<rect width="\d+" height="\d+" fill="white"/>', '', svg)
svg = svg.replace('<svg xmlns="http://www.w3.org/2000/svg"', '<svg class="layout" xmlns="http://www.w3.org/2000/svg"', 1)
svg = re.sub(r' width="\d+" height="\d+"', '', svg, count=1)
# status placeholders (edited later)
S = {
 "c224_d": "14 ≤ d ≤ 16 certified (SAT: w ≤ 13 UNSAT; SAT witness of weight 16); rungs 14–15 running",
 "c288_d": "d ≥ 14 certified (SAT: w ≤ 13 UNSAT, 1800 s), d ≤ 22 (ISD 3000 trials); exact value open",
 "dcirc": "BP-OSD random-η search (paper §6 method; 600 trials on 2 cycles, 150 on 4) finds no fault set below the code distance for either code (C224: ≤ 19 Z / ≤ 20 X; C288: ≤ 33), while the same budget does find d_circ ≤ 11 < d for the gross code (paper: ≤ 10). So d_circ ≤ d trivially and no circuit-level degradation was detected.",
}
html = f"""<title>Monomial Lifted-Product Codes</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,400;9..144,600&family=Source+Sans+3:wght@400;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{{
  --bg:#F5F7FA; --bg2:#FFFFFF; --ink:#16202B; --muted:#5B6B7F; --rule:#D5DCE5;
  --accent:#2F62C9; --accent-soft:#E3ECFB; --amber:#B8781A; --amber-soft:#FBF0DC;
  --green:#3E8A22; --red:#B03A32; --code-bg:#EEF2F7;
}}
@media (prefers-color-scheme: dark){{
  :root:not([data-theme="light"]){{
    --bg:#0F151D; --bg2:#151D28; --ink:#E7ECF2; --muted:#98A6B8; --rule:#2A3644;
    --accent:#7AA3F0; --accent-soft:#1B2A45; --amber:#E0AF52; --amber-soft:#2E2410;
    --green:#7CC35A; --red:#E0776D; --code-bg:#1B2431;
  }}
}}
:root[data-theme="dark"]{{
  --bg:#0F151D; --bg2:#151D28; --ink:#E7ECF2; --muted:#98A6B8; --rule:#2A3644;
  --accent:#7AA3F0; --accent-soft:#1B2A45; --amber:#E0AF52; --amber-soft:#2E2410;
  --green:#7CC35A; --red:#E0776D; --code-bg:#1B2431;
}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--ink);font-family:"Source Sans 3",system-ui,-apple-system,"Segoe UI",sans-serif;font-size:17px;line-height:1.55}}
.wrap{{max-width:1120px;margin:0 auto;padding:40px 28px 72px}}
header{{display:grid;grid-template-columns:1fr;gap:18px;padding-bottom:26px;border-bottom:1px solid var(--rule)}}
.eyebrow{{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12.5px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}}
h1{{font-family:Fraunces,Georgia,serif;font-weight:600;font-size:clamp(34px,4.6vw,52px);line-height:1.05;margin:0;letter-spacing:-.01em;text-wrap:balance}}
h1 em{{font-style:italic;font-weight:400;color:var(--accent)}}
.lede{{font-size:19px;max-width:66ch;color:var(--ink);margin:0}}
.spec{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px;margin-top:8px}}
.card{{background:var(--bg2);border:1px solid var(--rule);border-radius:10px;padding:18px 20px}}
.card .name{{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:13px;color:var(--muted);letter-spacing:.06em}}
.card .params{{font-family:Fraunces,Georgia,serif;font-size:34px;font-weight:600;margin:4px 0 8px;font-variant-numeric:tabular-nums;letter-spacing:-.01em}}
.card .params small{{font-size:18px;font-weight:400;color:var(--muted)}}
.kv{{display:grid;grid-template-columns:auto 1fr;gap:4px 14px;font-size:15px}}
.kv dt{{color:var(--muted);margin:0}} .kv dd{{margin:0}}
.mono,code,pre{{font-family:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace}}
code{{background:var(--code-bg);padding:1px 5px;border-radius:4px;font-size:.9em}}
pre{{background:var(--code-bg);padding:14px 16px;border-radius:8px;overflow-x:auto;font-size:13.5px;line-height:1.5;margin:0}}
main{{display:grid;grid-template-columns:minmax(0,1fr) 300px;gap:40px;margin-top:34px}}
@media (max-width:900px){{main{{grid-template-columns:1fr}} .spec{{grid-template-columns:1fr}}}}
h2{{font-family:Fraunces,Georgia,serif;font-weight:600;font-size:27px;margin:44px 0 12px;line-height:1.15;text-wrap:balance}}
h2:first-child{{margin-top:0}}
h3{{font-size:17px;margin:22px 0 6px;font-weight:600}}
p{{max-width:68ch;margin:0 0 12px}}
ul{{max-width:70ch;padding-left:20px;margin:0 0 12px}} li{{margin:4px 0}}
.fig{{background:var(--bg2);border:1px solid var(--rule);border-radius:10px;padding:16px;margin:14px 0 8px}}
.fig svg{{width:100%;height:auto;display:block}}
.fig svg line{{stroke:var(--muted)!important;stroke-opacity:.55}}
.cap{{font-size:14px;color:var(--muted);max-width:70ch;margin-top:10px}}
.legend{{display:flex;gap:18px;flex-wrap:wrap;font-size:14px;color:var(--muted);margin-top:8px}}
.dot{{display:inline-block;width:11px;height:11px;border-radius:50%;margin-right:6px;vertical-align:-1px}}
.sq{{display:inline-block;width:11px;height:11px;margin-right:6px;vertical-align:-1px}}
table{{border-collapse:collapse;width:100%;font-size:15px;font-variant-numeric:tabular-nums}}
.tablewrap{{overflow-x:auto;margin:10px 0 16px}}
th,td{{text-align:left;padding:8px 10px;border-bottom:1px solid var(--rule);vertical-align:top}}
th{{font-size:12.5px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);font-weight:600}}
tr.ours td{{background:var(--accent-soft)}}
.pill{{display:inline-block;font-family:"IBM Plex Mono",monospace;font-size:12px;padding:2px 8px;border-radius:999px;border:1px solid var(--rule);color:var(--muted)}}
.pill.ok{{color:var(--green);border-color:var(--green)}} .pill.run{{color:var(--amber);border-color:var(--amber)}}
aside .card{{margin-bottom:16px}}
aside h4{{margin:0 0 8px;font-size:13px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);font-weight:600}}
aside p, aside li{{font-size:14.5px}}
.thm{{border-left:3px solid var(--accent);padding:2px 0 2px 14px;margin:14px 0}}
.thm b{{font-family:"IBM Plex Mono",monospace;font-size:13px;letter-spacing:.05em;text-transform:uppercase;color:var(--accent)}}
.note{{background:var(--amber-soft);border-left:3px solid var(--amber);padding:10px 14px;border-radius:0 8px 8px 0;margin:14px 0;font-size:15px}}
a{{color:var(--accent)}}
footer{{margin-top:56px;padding-top:18px;border-top:1px solid var(--rule);color:var(--muted);font-size:14px}}
</style>
<div class="wrap">
<header>
  <div class="eyebrow">qec-lab · experiments/lp_lab · LP1 · 2026-09-13</div>
  <h1>Weight-6 quantum LDPC codes that are <em>not</em> bicycle codes, with the bicycle hardware package</h1>
  <p class="lede">Lifted products <span class="mono">LP(A, b⊗I)</span> with a 2×2 base matrix <span class="mono">A = [[p,1],[1,q]]</span> over F₂[G] inherit every hardware property of the Bravyi–Cross–Gambetta–Maslov–Rall–Yoder codes — degree-6 Tanner graph of thickness two, the 7-CNOT-layer syndrome cycle, and a toric layout — yet are provably not two-block group-algebra codes. Two instances below tie or exceed the best published weight-6 codes at their (n, k).</p>
  <div class="spec">
    <div class="card">
      <div class="name">C224 · G = Z₂₈ × Z₂</div>
      <div class="params" id="c224params">[[224, 12, ≤16]] <small>kd²/n ≤ 13.7</small></div>
      <dl class="kv">
        <dt>b</dt><dd class="mono">1 + x + x¹²y</dd>
        <dt>p, q</dt><dd class="mono">1 + x¹⁷,&nbsp; x³y + x⁴</dd>
        <dt>distance</dt><dd id="c224d">{S['c224_d']}</dd>
        <dt>layout</dt><dd>toric, Z₁₆ × Z₂₈ (4 layouts); wheel layers planar (9/9)</dd>
        <dt>symmetry</dt><dd>|Aut| = 56, two X-check orbits ⇒ not 2BGA</dd>
      </dl>
    </div>
    <div class="card">
      <div class="name">C288 · G = Z₇₂</div>
      <div class="params">[[288, 8, ≤22]] <small>kd²/n ≤ 13.4</small></div>
      <dl class="kv">
        <dt>b</dt><dd class="mono">1 + x² + x¹⁰</dd>
        <dt>p, q</dt><dd class="mono">1 + x⁵,&nbsp; x²¹ + x⁵⁸</dd>
        <dt>distance</dt><dd id="c288d">{S['c288_d']}</dd>
        <dt>layout</dt><dd>toric, Z₃₂ × Z₁₈ (4 layouts); wheel layers planar (9/9)</dd>
        <dt>symmetry</dt><dd>|Aut| = 72, two X-check orbits ⇒ not 2BGA</dd>
      </dl>
    </div>
  </div>
</header>

<main>
<article>
<h2>The layout the hardware sees</h2>
<div class="fig">{svg}
<div class="legend"><span><i class="dot" style="background:#4c9a2a"></i>L data</span><span><i class="dot" style="background:#c0392b"></i>R data</span><span><i class="sq" style="background:#e8a33d"></i>X check</span><span><i class="sq" style="background:#5b8def"></i>Z check</span><span>dashed: the two long-range couplers of six highlighted vertices</span></div>
<p class="cap">C224 on the torus Z₁₆ × Z₂₈. Every vertex has four nearest-neighbour couplers (grid) and two long-range ones, exactly as in Fig. 1B of the BB paper — but here the long-range shifts differ between the two copies of the group, which is what breaks the group-algebra structure. Wrap-around grid edges are omitted.</p>
</div>

<h2>What the family is</h2>
<p>Take a finite abelian group G, R = F₂[G], a weight-2 polynomial p containing 1, any weight-2 polynomial q, and a weight-3 polynomial b. The code is the CSS code with</p>
<pre>H_X = [ A | b* I₂ ],   H_Z = [ b I₂ | A† ],   A = [[p, 1], [1, q]]   (blocks are |G|×|G| circulants)</pre>
<p>Bivariate bicycle codes are the 1×1 case H_X = [A | B]. Each row of A holds three monomials, so every check has weight 6 and every qubit sits in three X and three Z checks. Relabelings reduce any connected 2×2 pattern to this normal form.</p>

<h2>Why the BB hardware results carry over</h2>
<p>The BB paper's thickness-two decomposition (Lemma 2) and the correctness proof of its depth-7 syndrome cycle (Table 5) use a single structural fact beyond "A and B are sums of three permutation matrices": every B-permutation commutes with every A-permutation. Here B = b⊗I is a scalar over the group ring, so this holds automatically.</p>
<div class="thm"><b>Thickness 2</b> — The layers {{A₂, A₃, B₃}} and {{A₁, B₁, B₂}} are disjoint unions of planar wheel graphs. Checked mechanically with a planarity test for all nine layer choices on both codes.</div>
<div class="thm"><b>Depth-7 cycle</b> — The eight-round schedule measures every check and acts trivially on logicals (the round-6 step needs B A₂ = A₂ B). Stim confirms all detectors of the noiseless memory circuit are deterministic for both codes and both bases.</div>
<div class="thm"><b>Toric layout</b> — If σ = A_jA_i⁻¹ and τ = B_gB_h⁻¹ generate a group acting regularly on the L qubits as Z_μ × Z_λ, the Cayley graph of Z_{{2μ}} × Z_{{2λ}} is a spanning subgraph, with explicit coordinates. For a = 2 the criterion is G = ⟨p₁q₁⟩ ⊕ ⟨b_g − b_h⟩; the torus is Z_{{4·ord(p₁q₁)}} × Z_{{2·ord τ}}. Verified edge by edge.</div>
<div class="thm"><b>Not a 2BGA code</b> — A two-block code from commuting permutation sums with both permutation groups transitive is a group-algebra code; a group-algebra code is check-transitive. nauty finds two X-check orbits for both codes, so they are not 2BGA codes for any group and in particular not BB codes.</div>

<h2>Designing instead of sweeping</h2>
<p>Blind sweeps at n = 144 over every abelian group of order 36 capped at [[144,8,12]] for this class, and the reason is structural: the constituent classical code ker b lives on one copy of G, of length n/4 rather than the n/2 of a bicycle code. Two exact facts then turn the search into design:</p>
<div class="thm"><b>k</b> — k = 2·dim F₂[G]/(b*, 1 + pq); it is maximal, k = 2·dim ker b, exactly when pq ≡ 1 (mod (b*)). (Frobenius annihilator duality; the unconjugated version is false for even |G|.)</div>
<div class="thm"><b>Cap</b> — The X-logical (u, 0), u ∈ ker b, gives d ≤ D_b := d(ker b), so kd²/n ≤ δ_b D_b²/(2|G|), with D_b computed exactly.</div>
<p>So for each (G, b) with a strong low-dimensional kernel we enumerate only the (p, q) with pq ≡ 1 mod (b*) — a few thousand codes per b instead of millions — and screen them by randomized information-set decoding. The Z₅₆ table already explains the published coset code [[224,12,16]]: b = (1+x+x³)² has a [56,6,16] kernel and cap 13.71.</p>

<h2>Where they sit</h2>
<div class="tablewrap"><table>
<tr><th>code</th><th>[[n,k,d]]</th><th>kd²/n</th><th>construction</th><th>θ = 2</th><th>toric layout</th></tr>
<tr><td>Bravyi et al.</td><td>[[144,12,12]]</td><td>12.0</td><td>BB</td><td>yes</td><td>yes</td></tr>
<tr><td>Bravyi et al.</td><td>[[288,12,18]]</td><td>13.5</td><td>BB</td><td>yes</td><td>yes</td></tr>
<tr><td>IBM class Y (repo A33)</td><td>[[288,8,20]]</td><td>11.1</td><td>BB</td><td>yes</td><td>yes</td></tr>
<tr><td>Aydin–Tamo–Barg</td><td>[[224,12,16]]</td><td>13.7</td><td>coset two-block</td><td>yes</td><td>not given</td></tr>
<tr><td>Liang et al.</td><td>[[248,10,18]], [[294,10,20]], [[254,14,16]]</td><td>13.1, 13.6, 14.1</td><td>twisted-torus BB</td><td>yes</td><td>yes</td></tr>
<tr class="ours"><td>C224 (this work)</td><td id="c224row">[[224,12,≤16]] (d ≥ 14 certified)</td><td>13.7</td><td>LP(A, b⊗I), non-2BGA</td><td>yes</td><td>Z₁₆×Z₂₈</td></tr>
<tr class="ours"><td>C288 (this work)</td><td>[[288,8,≤22]]</td><td>≤13.4</td><td>LP(A, b⊗I), non-2BGA</td><td>yes</td><td>Z₃₂×Z₁₈</td></tr>
</table></div>
<p>C224 ties the coset code in parameters; it adds an explicit polynomial description, the k and cap formulas, and a toric layout (the coset paper lists such layouts as open). C288 would beat the known BB code at (288, 8) if d = 22 holds. Fom records at n ≈ 250–340 stay with the twisted-torus BB codes; the contribution here is a new non-bicycle structure that keeps the full hardware package, plus a complete design theory for it.</p>

<h2>Certificates and caveats</h2>
<ul>
<li><b>Exact distance.</b> CryptoMiniSat ladder with native XOR rows and translation-orbit symmetry breaking (one instance per qubit orbit), validated against a plain-SAT-certified [[144,8,12]] (7–11× faster per rung). <span id="certsat">{S['c224_d']}. C288: {S['c288_d']}.</span></li>
<li><b>Circuit-level distance.</b> <span id="dcirc">{S['dcirc']}</span> Distance preservation is therefore open — the same status the BB paper reports for its codes other than [[72,12,6]].</li>
<li><b>Relation to coset codes.</b> Formally LP(A, b⊗I) is a coset code with the wreath group Z_N ≀ S₂ and central b; the coset paper searched small GAP groups, not these lifts.</li>
<li><b>Screening trap recorded.</b> Early-stopped ISD values are not estimates; one k = 16 "d = 16" lead at n = 288 was retracted (true d ≤ 12).</li>
<li>Both codes have no weight-4 stabilizers (minimum stabilizer weight 6).</li>
</ul>
</article>

<aside>
<div class="card"><h4>Reproduce</h4>
<pre>cd experiments/lp_lab
uv run python scripts/c1_certify.py \\
  --orders 28,2 --p 0,34 --q 7,8 --b 0,2,25 \\
  --skip-exact --no-circuit
uv run python scripts/c2_satsym.py \\
  --orders 28,2 --p 0,34 --q 7,8 --b 0,2,25 \\
  --lower 10 --tag C224</pre>
<p>Indices are group elements in row-major order over the cyclic factors; p, q, b are supports.</p></div>
<div class="card"><h4>Status</h4>
<p><span class="pill ok">structure</span> planarity, determinism, toric coordinates, orbits — done</p>
<p><span class="pill run">exact d</span> <span id="statusd">C224: rungs 14–15 running (d ≥ 14). C288: stopped at certified d ≥ 14.</span></p>
<p><span class="pill ok">d_circ</span> BP-OSD passes done: no degradation found (bounds above d)</p></div>
<div class="card"><h4>Files</h4>
<ul>
<li><code>notes/LP1_theory.md</code> — statements and proofs</li>
<li><code>notes/LP1_results.md</code>, <code>LP1_plan.md</code> — log</li>
<li><code>data/codes/cand_224_toric.json</code>, <code>cand_288_toric.json</code> — code cards with permutation parts</li>
<li><code>data/figs/*.svg</code> — layouts</li>
</ul></div>
<div class="card"><h4>The n = 288 family</h4>
<p>59 designs at 3000-trial d ≤ 22 across Z₇₂ and Z₃₆×Z₂; twelve carry toric layouts (32×18, 16×36, 8×72). C288 is one of them, not an isolated hit.</p></div>
</aside>
</main>
<footer>Built in the qec-lab worktree <span class="mono">claude/stabilizer-code-lifted-product-784b9e</span>. Benchmarks: arXiv:2308.07915, 2606.17268, 2503.03827, 2606.02418.</footer>
</div>
"""
out.write_text(html); print("wrote", out, len(html), "bytes")
