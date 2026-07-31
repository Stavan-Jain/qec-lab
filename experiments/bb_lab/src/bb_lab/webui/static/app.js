"use strict";

const $ = (id) => document.getElementById(id);

const state = {
  backend: null,        // /api/backend payload
  report: null,         // last successful analysis
  choice: "tandem",     // which backend to run
  mode: "naive",
  flags: {},            // flag -> value (only what the user turned on)
  ack: new Set(),       // flags whose caller obligation was acknowledged
  job: null,
  source: null,
  solvedFor: null,      // signature of the code the solve card describes
};

// A distance belongs to one specific (G, A, B). Anything the solver said
// about a different code has to leave the screen the moment the inputs
// change, or the headline quietly attributes an old d to a new code.
const signature = () =>
  [$("orders").value, $("polyA").value, $("polyB").value].join("|");

// ───────────────────────────────────────────────────────── helpers

async function api(path, body) {
  const res = await fetch(path, {
    method: body === undefined ? "GET" : "POST",
    headers: body === undefined ? {} : { "Content-Type": "application/json" },
    body: body === undefined ? undefined : JSON.stringify(body),
  });
  const data = await res.json().catch(() => ({}));
  if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
  return data;
}

const debounce = (fn, ms) => {
  let t;
  return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); };
};

const el = (tag, cls, text) => {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (text !== undefined) n.textContent = text;
  return n;
};

const fmtSecs = (s) =>
  s < 60 ? `${s.toFixed(s < 10 ? 2 : 1)}s`
         : s < 3600 ? `${Math.floor(s / 60)}m ${Math.round(s % 60)}s`
                    : `${Math.floor(s / 3600)}h ${Math.round((s % 3600) / 60)}m`;

// ───────────────────────────────────────────────────────── backend

async function loadBackend() {
  state.backend = await api("/api/backend");
  const t = state.backend.tandem;
  const badge = $("backend-badge");

  if (t.available && t.fork) {
    badge.className = "badge badge-good";
    badge.textContent = `Tandem ready · ${t.options.length} solver options`;
    badge.title = t.path;
  } else if (t.available) {
    badge.className = "badge badge-warn";
    badge.textContent = "stock MaxCDCL (not the fork)";
    badge.title = t.error || t.path;
  } else {
    badge.className = "badge badge-warn";
    badge.textContent = "Tandem unavailable — using SAT fallback";
    badge.title = t.error || "";
    state.choice = "sat-ladder";
  }

  const sel = $("backend-select");
  sel.replaceChildren();
  const opts = [
    { id: "tandem", label: "Tandem", ok: t.available },
    { id: "sat-ladder", label: "SAT ladder", ok: true },
  ];
  for (const o of opts) {
    const b = el("button", state.choice === o.id ? "on" : "", o.label);
    b.type = "button";
    b.disabled = !o.ok;
    if (!o.ok) b.title = t.error || "";
    b.onclick = () => {
      state.choice = o.id;
      [...sel.children].forEach((c) => c.classList.toggle("on", c === b));
      $("mode-field").classList.toggle("hidden", o.id !== "tandem");
      renderFlags();
    };
    sel.appendChild(b);
  }
  $("mode-field").classList.toggle("hidden", state.choice !== "tandem");
  renderFlags();
}

async function loadPresets() {
  const { presets } = await api("/api/presets");
  const box = $("presets");
  box.replaceChildren();
  for (const p of presets) {
    const b = el("button", "", p.name.replace(" gross", " ✦"));
    b.type = "button";
    b.title = `${p.orders}  A = ${p.A}  B = ${p.B}`;
    b.onclick = () => {
      $("orders").value = p.orders;
      $("polyA").value = p.A;
      $("polyB").value = p.B;
      analyse();
    };
    box.appendChild(b);
  }
}

// ───────────────────────────────────────────────────────── flags

function renderFlags() {
  const box = $("flags");
  box.replaceChildren();
  if (state.choice !== "tandem") return;
  const t = state.backend?.tandem;
  if (!t?.available) return;

  const premises = state.report?.premises || {};
  const featured = t.options.filter((o) => o.featured);
  const rest = t.options.filter((o) => !o.featured);

  const list = el("div", "flaglist");
  for (const opt of featured) list.appendChild(featuredFlag(opt, premises));
  box.appendChild(list);

  if (rest.length) {
    const d = el("details", "advanced");
    d.appendChild(el("summary", "", `all ${t.options.length} solver options`));
    const grid = el("div", "advanced-grid");
    for (const opt of rest) grid.appendChild(advancedFlag(opt));
    d.appendChild(grid);
    box.appendChild(d);
  }
}

function featuredFlag(opt, premises) {
  const wrap = el("div", "flag");
  const head = el("div", "flag-head");
  const cb = el("input");
  cb.type = "checkbox";
  cb.checked = opt.flag in state.flags;

  head.appendChild(cb);
  head.appendChild(el("span", "flag-name", opt.label));
  head.appendChild(el("span", "flag-code", opt.flag));
  wrap.appendChild(head);
  if (opt.blurb) wrap.appendChild(el("div", "flag-blurb", opt.blurb));

  // A flag with a machine-checkable premise is only offered when the
  // premise holds for the code on screen.
  const premise = opt.requires ? premises[opt.requires] : null;
  const blocked = opt.requires && premise && !premise.holds;
  if (blocked) {
    wrap.classList.add("blocked");
    cb.checked = false;
    cb.disabled = true;
    delete state.flags[opt.flag];
    const r = el("div", "flag-block-reason");
    r.textContent = `Unavailable: “${premise.label}” is false here. ${premise.detail}`;
    wrap.appendChild(r);
    return wrap;
  }
  if (opt.requires && !premise) {
    cb.disabled = true;
    wrap.classList.add("blocked");
    wrap.appendChild(el("div", "flag-block-reason", "Enter a valid code to check this flag's premise."));
    return wrap;
  }

  let valueInput = null;
  if (opt.kind !== "bool") {
    const row = el("div", "flag-row");
    valueInput = el("input");
    valueInput.type = opt.kind === "string" ? "text" : "number";
    valueInput.value = state.flags[opt.flag] ?? opt.suggested ?? opt.default ?? "";
    valueInput.oninput = () => {
      if (cb.checked) state.flags[opt.flag] = valueInput.value;
    };
    row.appendChild(valueInput);
    if (opt.domain) row.appendChild(el("span", "flag-blurb", opt.domain));
    wrap.appendChild(row);
  }

  // No machine check exists for this obligation → demand an explicit ack.
  let ackBox = null;
  if (opt.soundness === "caller" && !opt.requires) {
    const lab = el("label", "flag-ack");
    ackBox = el("input");
    ackBox.type = "checkbox";
    ackBox.checked = state.ack.has(opt.flag);
    ackBox.onchange = () => {
      ackBox.checked ? state.ack.add(opt.flag) : state.ack.delete(opt.flag);
      sync();
    };
    lab.appendChild(ackBox);
    lab.appendChild(el("span", "", "I certify this value — an unsound one silently returns a wrong distance."));
    wrap.appendChild(lab);
  }

  const sync = () => {
    if (!cb.checked) { delete state.flags[opt.flag]; return; }
    state.flags[opt.flag] = opt.kind === "bool"
      ? true
      : (valueInput.value || opt.suggested || opt.default);
  };
  cb.onchange = sync;
  sync();
  return wrap;
}

function advancedFlag(opt) {
  const row = el("div", "advanced-row");
  const lab = el("label", "", opt.flag);
  row.appendChild(lab);
  if (opt.domain) row.appendChild(el("span", "dom", opt.domain));

  if (opt.kind === "bool") {
    const cb = el("input");
    cb.type = "checkbox";
    cb.checked = state.flags[opt.flag] === true;
    cb.onchange = () => {
      if (cb.checked) state.flags[opt.flag] = true;
      else delete state.flags[opt.flag];
    };
    row.appendChild(cb);
  } else {
    const inp = el("input");
    inp.type = opt.kind === "string" ? "text" : "number";
    inp.placeholder = opt.default ?? "";
    inp.value = state.flags[opt.flag] ?? "";
    inp.oninput = () => {
      if (inp.value === "") delete state.flags[opt.flag];
      else state.flags[opt.flag] = inp.value;
    };
    row.appendChild(inp);
  }
  return row;
}

// ───────────────────────────────────────────────────── analysis

// Analysis cost grows with |G| (~1 s at n = 288), so a fast request issued
// after a slow one can come back first. Stamp each request and drop any
// response that is not the newest, or the page paints an older code's report.
let analyseSeq = 0;

const analyse = debounce(async () => {
  const mine = ++analyseSeq;
  const body = {
    orders: $("orders").value,
    A: $("polyA").value,
    B: $("polyB").value,
  };
  try {
    const report = await api("/api/analyse", body);
    if (mine !== analyseSeq) return;
    state.report = report;
    $("input-error").classList.add("hidden");
    renderReport(report);
  } catch (e) {
    if (mine !== analyseSeq) return;
    state.report = null;
    const box = $("input-error");
    box.textContent = e.message;
    box.classList.remove("hidden");
    $("p-n").textContent = "·";
    $("p-k").textContent = "·";
    $("check-weight").textContent = "·";
  }
  renderFlags();
}, 220);

function renderReport(r) {
  $("group-hint").innerHTML =
    `<code>${r.group_label}</code> · |G| = ${r.group_order} · n = 2|G| = ${r.n}`;

  $("p-n").textContent = r.n;
  $("p-k").textContent = r.k;

  // d belongs to the solver, and only to the code it was solved for.
  const sig = signature();
  if (state.solvedFor !== null && state.solvedFor !== sig) {
    state.solvedFor = null;
    $("solve-card").hidden = true;
  }
  // A run in flight is answering a question the user has stopped asking.
  if (state.job && state.jobSig !== sig) {
    api(`/api/solve/${state.job}/cancel`, {}).catch(() => {});
  }
  if (!state.job) setDistance(null);

  const cw = $("check-weight");
  const sub = $("check-weight-sub");
  if (r.check_weight !== null) {
    cw.textContent = r.check_weight;
    sub.textContent = `wt(A) ${r.A_weight} + wt(B) ${r.B_weight}`;
  } else {
    cw.textContent = `${r.x_profile.row_min}–${r.x_profile.row_max}`;
    sub.textContent = "non-uniform rows";
  }

  const warn = $("warnings");
  warn.replaceChildren();
  for (const w of r.warnings) warn.appendChild(el("div", "notice", w));

  const qubitDeg = r.qubit_degree !== null
    ? { v: r.qubit_degree, note: "X-checks per qubit (Z alike)" }
    : { v: `${r.x_profile.col_min}–${r.x_profile.col_max}`, note: "non-uniform" };

  const stats = [
    ["Group", r.group_label, `|G| = ${r.group_order}`],
    ["Qubits n", r.n, "2|G|"],
    ["Logicals k", r.k, `n − rank H_X − rank H_Z`],
    ["Check weight", r.check_weight ?? `${r.x_profile.row_min}–${r.x_profile.row_max}`, "X and Z alike"],
    ["Qubit degree", qubitDeg.v, qubitDeg.note],
    ["Stabilizers", r.num_checks, `${r.group_order} X + ${r.group_order} Z`],
    ["rank H_X", r.rank_HX, ""],
    ["rank H_Z", r.rank_HZ, ""],
    ["Rate k/n", r.rate.toFixed(4), ""],
    ["CSS commute", r.css_commutes ? "✓" : "✗", "H_X · H_Zᵀ = 0", r.css_commutes ? "ok" : "bad"],
  ];
  const grid = $("stats");
  grid.replaceChildren();
  for (const [label, value, note, cls] of stats) {
    const s = el("div", "stat");
    s.appendChild(el("div", "stat-label", label));
    s.appendChild(el("div", `stat-value${cls ? " " + cls : ""}`, String(value)));
    if (note) s.appendChild(el("div", "stat-note", note));
    grid.appendChild(s);
  }

  renderPremises(r.premises);
  renderKnown(r.known);
  $("solve").disabled = r.k === 0;
}

function renderPremises(premises) {
  const card = $("premise-card");
  const body = $("premise-body");
  body.replaceChildren();
  const keys = Object.keys(premises || {});
  if (!keys.length) { card.hidden = true; return; }
  card.hidden = false;
  for (const key of keys) {
    const p = premises[key];
    const row = el("div", "premise");
    row.appendChild(el("span", `mark ${p.holds ? "yes" : "no"}`, p.holds ? "✓" : "✗"));
    const txt = el("div");
    txt.appendChild(el("div", "", p.label));
    txt.appendChild(el("div", "detail", p.detail));
    row.appendChild(txt);
    body.appendChild(row);
  }
}

function renderKnown(known) {
  const card = $("known-card");
  const body = $("known-body");
  if (!known || !known.found) { card.hidden = true; return; }
  card.hidden = false;
  body.replaceChildren();
  const dl = el("dl", "kv");
  const add = (k, v) => {
    if (v === null || v === undefined) return;
    dl.appendChild(el("dt", "", k));
    dl.appendChild(el("dd", "", String(v)));
  };
  add("code_id", known.code_id);
  add("d (exact)", known.d_exact);
  add("d bounds", known.d_exact === null ? `${known.d_lb ?? "?"} … ${known.d_ub ?? "?"}` : null);
  add("method", known.d_method);
  add("orbit size", known.orbit_size);
  body.appendChild(dl);
  body.appendChild(el(
    "div", "field-hint",
    "Matched on the canonical orbit representative under Aut(G) ⋉ G plus " +
    "block swap, so an equivalent presentation still finds this row. " +
    "Cited, not re-derived — press Compute distance to prove it here.",
  ));
  if (known.d_exact !== null && known.d_exact !== undefined) {
    setDistance(known.d_exact, "known");
  }
}

function setDistance(d, kind) {
  const n = $("p-d");
  n.className = "";
  if (d === null || d === undefined) {
    n.textContent = "?";
    n.classList.add("unknown");
    n.title = "";
    return;
  }
  n.textContent = d;
  // A value read out of the corpus and one this session actually proved
  // must not look the same — the second is evidence, the first is a claim
  // someone else's run made.
  if (kind === "known") {
    n.classList.add("cited");
    n.title = "from the corpus — cited, not solved in this session";
  } else {
    n.title = "solved and witness-verified in this session";
  }
}

// ─────────────────────────────────────────────────────── solving

async function solve() {
  if (!state.report) return;
  const sig = signature();
  state.solvedFor = null;
  const payload = {
    orders: $("orders").value,
    A: $("polyA").value,
    B: $("polyB").value,
    backend: state.choice,
    mode: state.mode,
    flags: state.choice === "tandem" ? state.flags : {},
    acknowledged: [...state.ack],
  };

  $("solve-card").hidden = false;
  $("incumbents").replaceChildren();
  $("log").textContent = "";
  $("solve-summary").textContent = "";
  setState("starting", "badge badge-run");

  let job;
  try {
    job = await api("/api/solve", payload);
  } catch (e) {
    setState("rejected", "badge badge-bad");
    $("solve-summary").textContent = e.message;
    return;
  }

  state.job = job.job_id;
  state.jobSig = sig;
  $("solve").disabled = true;
  $("cancel").classList.remove("hidden");
  $("solve-title").textContent =
    job.backend === "tandem" ? "Tandem" : "CryptoMiniSat ladder";
  const d = $("p-d");
  d.className = "pending";
  d.textContent = "…";

  const es = new EventSource(`/api/solve/${job.job_id}/events`);
  state.source = es;

  es.addEventListener("stage", (ev) => {
    const p = JSON.parse(ev.data);
    setState(p.stage, "badge badge-run");
    appendLog(`— ${p.stage}: ${p.detail}`);
  });

  es.addEventListener("incumbent", (ev) => {
    const p = JSON.parse(ev.data);
    addIncumbent(p.cost, p.elapsed);
    const node = $("p-d");
    node.className = "pending";
    node.textContent = `≤${p.cost}`;
  });

  es.addEventListener("rung", (ev) => {
    const p = JSON.parse(ev.data);
    addRung(p.weight, p.sat, p.seconds);
    const node = $("p-d");
    node.className = "pending";
    node.textContent = p.sat ? p.weight : `>${p.weight}`;
  });

  es.addEventListener("log", (ev) => appendLog(JSON.parse(ev.data).line));

  es.addEventListener("done", (ev) => {
    const p = JSON.parse(ev.data);
    setState("optimum", "badge badge-good");
    state.solvedFor = sig;
    setDistance(p.distance);
    const bits = [
      `d = <strong>${p.distance}</strong>`,
      `witness weight <strong>${p.witness_weight}</strong> (re-verified: in ker H_Z, not a stabilizer)`,
      `solver <strong>${fmtSecs(p.solver_seconds)}</strong>`,
      `method <strong>${p.method}</strong>`,
    ];
    $("solve-summary").innerHTML = bits.join(" · ");
    if (state.report) {
      const q = (state.report.k * p.distance ** 2) / state.report.n;
      $("solve-summary").innerHTML +=
        ` · merit q = k·d²/n = <strong>${q.toFixed(2)}</strong>`;
    }
    finish(es);
  });

  es.addEventListener("cancelled", () => {
    setState("cancelled", "badge badge-warn");
    setDistance(null);
    finish(es);
  });

  es.addEventListener("error", (ev) => {
    if (ev.data) {
      const p = JSON.parse(ev.data);
      setState("failed", "badge badge-bad");
      $("solve-summary").textContent = p.message;
      appendLog(p.message);
      setDistance(null);
      finish(es);
    }
  });
}

function finish(es) {
  es.close();
  state.source = null;
  state.job = null;
  $("solve").disabled = false;
  $("cancel").classList.add("hidden");
}

function setState(text, cls) {
  const b = $("solve-state");
  b.textContent = text;
  b.className = cls;
}

function addIncumbent(cost, elapsed) {
  const box = $("incumbents");
  [...box.children].forEach((c) => c.classList.remove("latest"));
  const n = el("span", "inc latest");
  n.appendChild(el("span", "", `≤ ${cost}`));
  n.appendChild(el("span", "t", `${elapsed}s`));
  box.appendChild(n);
}

function addRung(weight, sat, seconds) {
  const box = $("incumbents");
  [...box.children].forEach((c) => c.classList.remove("latest"));
  const n = el("span", "inc latest");
  n.appendChild(el("span", "", `w ≤ ${weight}: ${sat ? "SAT" : "UNSAT"}`));
  n.appendChild(el("span", "t", `${seconds}s`));
  box.appendChild(n);
}

function appendLog(line) {
  const log = $("log");
  log.textContent += (log.textContent ? "\n" : "") + line;
  log.scrollTop = log.scrollHeight;
}

// ──────────────────────────────────────────────────────── wiring

// Editing the code retracts the displayed d immediately, before the (possibly
// slow) analysis comes back. Otherwise the previous code's distance sits in
// the headline next to the new code's n and k.
function invalidateOnEdit() {
  if (state.job) return;              // a live run manages the headline itself
  if (state.solvedFor !== null && state.solvedFor === signature()) return;
  setDistance(null);
}

for (const id of ["orders", "polyA", "polyB"]) {
  $(id).addEventListener("input", () => { invalidateOnEdit(); analyse(); });
}

$("solve").onclick = solve;
$("cancel").onclick = async () => {
  if (state.job) await api(`/api/solve/${state.job}/cancel`, {});
};

for (const b of $("mode-select").children) {
  b.onclick = () => {
    state.mode = b.dataset.mode;
    [...$("mode-select").children].forEach((c) => c.classList.toggle("on", c === b));
  };
}

(async () => {
  await loadBackend();
  await loadPresets();
  analyse();
})();
