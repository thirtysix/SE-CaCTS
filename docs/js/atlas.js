/* atlas.js, browse specific super-enhancers per group, at each resolution.
   Enforces the resolution rule: CALLS at lineage/disease, RANKINGS ONLY at subtype. At cell-line level a line
   is compared four ways (vs all lines, and vs the lines of its subtype, primary disease and lineage), with calls
   where the line has two independent studies and enough relatives (data/lines/<key>.json). */
const Atlas = (() => {
  let manifest, level = "lineage", group = null, rows = [], cache = {};
  // cell-line level: the per-line comparison files, and the comparison on show
  let lineIdx = null, keyOf = {}, lineData = null, passSel = null;
  // the atlas names a line as scored ("C4-2B", "EA.hy926"); the line index keys it stripped ("C42B"), as DepMap does
  const lineKey = g => keyOf[String(g).toUpperCase().replace(/[^A-Z0-9]/g, "")];
  const CMPS = [["all", "All", "all lines"], ["lineage", "Lineage", "lineage"],          // broad to narrow
                ["disease", "Disease", "primary disease"], ["subtype", "Subtype", "subtype"]];
  let baseHead = null;                                   // the call-level table header, restored off the line level
  let sortKey = "rank", sortAsc = true, geneQuery = "", fdrMax = 1;
  // Hierarchy scoping. Every scored line has one lineage, primary disease and subtype (manifest.hierarchy),
  // so a level can be restricted to the groups inside a choice made at the levels above it.
  const ORDER = ["lineage", "disease", "subtype", "line"];
  let hier = [], scope = { lineage: null, disease: null, subtype: null };
  // the scoping can be switched off to list every group at a level at once (the original picker);
  // remembered per browser, and the page works the same when storage is unavailable
  let useScope = true;
  try { useScope = localStorage.getItem("secacts-atlas-scope") !== "off"; } catch (_) { /* default on */ }
  const inScope = (sc, upto) => hier.filter(h => ORDER.slice(0, upto).every(u => !sc[u] || h[u] === sc[u]));
  const linesOf = (lv, g) => hier.filter(h => h[lv] === g);
  const unique = (rows, lv) => { const s = new Set(rows.map(h => h[lv])); return s.size === 1 ? [...s][0] : null; };
  const majority = (rows, lv) => { const c = {}; rows.forEach(h => { c[h[lv]] = (c[h[lv]] || 0) + 1; });
    return Object.keys(c).sort((a, b) => c[b] - c[a])[0] || null; };
  const lab = (lv, g) => (manifest.levels[lv].groups[g] || {}).label || g;

  const levelDef = k => U.LEVELS.find(l => l.key === k);
  // analysis variant: applies to the CALL levels only; subtype rankings and cell-line comparisons use the main run
  const variantOf = k => levelDef(k).kind === "calls" ? U.getVariant() : "main";
  const vinfo = () => (manifest.variants || []).find(v => v.key === U.getVariant()) || null;
  const ck = k => `${k}|${variantOf(k)}`;                                // cache key
  const ginfo = (k, g) => {
    const v = variantOf(k), vi = v !== "main" && (manifest.variants || []).find(x => x.key === v);
    return (vi && vi.groups && vi.groups[k] && vi.groups[k][g]) || manifest.levels[k].groups[g] || {};
  };
  const fmtLen = bp => bp == null ? "" : (bp < 1000 ? `${bp} bp` : `${(bp / 1000).toFixed(1)} kb`);
  // The staged FDR is rounded to 4 dp, so report anything below 0.001 as a bound rather than as a
  // spuriously precise number (at cell-line level the degenerate null rounds many rows to 0).
  const fmtFdrVal = v => {
    if (v == null || v === "" || isNaN(+v)) return "n/a";
    v = +v;
    return v < 0.001 ? "<0.001" : v.toFixed(3);
  };

  async function loadLevel(k) {
    if (k === "line") {
      if (!lineIdx) {
        lineIdx = await DataLoader.loadJSON("data/lines/index.json");
        keyOf = Object.fromEntries(Object.entries(lineIdx).map(([key, v]) => [v.group, key]));
      }
      cache[k] = cache[k] || [];
      return cache[k];
    }
    const key = ck(k);
    if (!cache[key]) {
      const r = (await DataLoader.loadTSV(U.variantFile(levelDef(k).file, variantOf(k)))).rows;
      r.forEach(x => { x.len = (x.end || 0) - (x.start || 0); });   // SE span, for the Length column + sort
      cache[key] = r;
    }
    return cache[key];
  }

  function groupsFor(k) {
    const g = manifest.levels[k].groups, i = ORDER.indexOf(k);
    const allowed = useScope && hier.length && i > 0 ? new Set(inScope(scope, i).map(h => h[k])) : null;
    return Object.keys(g).filter(x => !allowed || allowed.has(x))
      .sort((a, b) => lab(k, a).localeCompare(lab(k, b)));                  // alphabetical, by display name
  }

  // Moving between levels keeps the context: drilling down makes the current group the scope (Ovary at
  // lineage level -> only the ovarian cell lines); going up selects the current group's parent.
  function switchLevel(nl) {
    const oi = ORDER.indexOf(level), ni = ORDER.indexOf(nl), anchor = group ? linesOf(level, group) : [];
    if (!useScope) group = null;                 // original behaviour: a fresh list at the new level
    else if (anchor.length) {
      if (ni > oi) {
        ORDER.slice(0, oi).forEach(u => { if (u in scope) scope[u] = unique(anchor, u) || scope[u]; });
        if (level in scope) scope[level] = group;
        ORDER.slice(oi + 1).forEach(u => { if (u in scope) scope[u] = null; });
        group = null;
      } else {
        group = majority(anchor, nl);
        ORDER.forEach((u, i) => { if (u in scope) scope[u] = i < ni ? unique(anchor, u) : null; });
      }
    }
    level = nl;
    onLevel();
  }

  // jump to a group at another level (breadcrumb), with the levels above it scoped to its parents
  function goTo(lv, g) {
    const rows = linesOf(lv, g), ni = ORDER.indexOf(lv);
    ORDER.forEach((u, i) => { if (u in scope) scope[u] = i < ni ? unique(rows, u) : null; });
    level = lv; group = g;
    onLevel();
  }

  function renderScope() {
    const idx = ORDER.indexOf(level), el = U.el("atlas-scope");
    if (!hier.length || idx === 0) { el.style.display = "none"; el.innerHTML = ""; return; }
    el.style.display = "flex";
    const toggle = `<label class="scope-toggle" title="on: list only the groups inside a chosen lineage / disease / subtype. Off: list every ${levelDef(level).label.toLowerCase()} group at once">
      <input type="checkbox" id="atlas-scope-on"${useScope ? " checked" : ""}> Narrow by parent group</label>`;
    const wire = () => {
      U.el("atlas-scope-on").onchange = e => {
        useScope = e.target.checked;
        try { localStorage.setItem("secacts-atlas-scope", useScope ? "on" : "off"); } catch (_) { /* not persisted */ }
        renderScope(); renderGroupPicker(); renderTable();
      };
    };
    if (!useScope) { el.innerHTML = toggle; wire(); return; }
    const sels = ORDER.slice(0, idx).map((u, i) => {
      const rows = inScope(scope, i), cnt = {};
      rows.forEach(h => { if (h[u]) cnt[h[u]] = (cnt[h[u]] || 0) + 1; });
      const lbl = levelDef(u).label.toLowerCase();
      return `<select class="sel scope-sel" data-u="${u}" title="show only the ${levelDef(level).label.toLowerCase()} groups inside one ${lbl}">
        <option value="">All ${lbl}s (${rows.length} lines)</option>
        ${Object.keys(cnt).sort((a, b) => a.localeCompare(b)).map(o =>
          `<option value="${U.esc(o)}"${scope[u] === o ? " selected" : ""}>${U.esc(o)} (${cnt[o]} line${cnt[o] > 1 ? "s" : ""})</option>`).join("")}
      </select>`;
    });
    const any = ORDER.slice(0, idx).some(u => scope[u]);
    el.innerHTML = toggle + `<span class="ctl-label" title="narrow this level to the groups inside a lineage, primary disease or subtype">Within</span>`
      + sels.join(`<span class="scope-sep">›</span>`)
      + (any ? `<button class="dl-btn scope-clear" title="show every group at this level">clear</button>` : "");
    wire();
    el.querySelectorAll(".scope-sel").forEach(sel => sel.onchange = () => {
      const u = sel.dataset.u, i = ORDER.indexOf(u);
      scope[u] = sel.value || null;
      // a narrower choice below this one survives only if it still fits
      ORDER.slice(i + 1, idx).forEach(d => {
        if (scope[d] && !inScope(scope, ORDER.indexOf(d)).some(h => h[d] === scope[d])) scope[d] = null;
      });
      // and a broader one above follows the new choice when it is unambiguous
      if (scope[u]) ORDER.slice(0, i).forEach(a => { const p = unique(linesOf(u, scope[u]), a); if (p) scope[a] = p; });
      renderScope(); renderGroupPicker(); renderTable();
    });
    const clr = el.querySelector(".scope-clear");
    if (clr) clr.onclick = () => { Object.keys(scope).forEach(k => { scope[k] = null; }); renderScope(); renderGroupPicker(); renderTable(); };
  }

  function renderControls() {
    U.el("atlas-level").innerHTML = U.LEVELS.map(l =>
      `<button data-lv="${l.key}" aria-selected="${l.key === level}" title="${l.kind === "calls"
        ? "specificity calls are supported at this resolution" : l.kind === "lines"
        ? "one cell line compared with all lines, and with the lines of its subtype, primary disease and lineage"
        : "rankings only: the panel does not support calls here (few lines per group)"}">${l.label}</button>`).join("");
    U.el("atlas-level").querySelectorAll("button").forEach(b =>
      b.onclick = () => { if (b.dataset.lv !== level) switchLevel(b.dataset.lv); });
  }

  function renderFdrControl() {
    // FDR threshold filter where calls exist (gotcha 72): lineage / disease, and the cell-line table's "vs all"
    // column. Shown but disabled at subtype, where the panel supports rankings only, so the bar never jumps.
    const kind = levelDef(level).kind, active = kind === "calls" || kind === "lines";
    U.el("atlas-fdr-wrap").style.display = "flex";
    if (!active) fdrMax = 1;
    const why = "subtype shows rankings only: most subtypes hold too few cell lines for a significance call";
    const opts = [[1, "All"], [0.10, "≤ 0.10"], [0.05, "≤ 0.05"], [0.01, "≤ 0.01"]];
    U.el("atlas-fdr").innerHTML = opts.map(([v, l]) =>
      `<button data-fdr="${v}" aria-selected="${v === fdrMax}"${active ? "" : " disabled"} title="${!active ? why : v === 1
        ? "show every call (all pass the ≤ 0.10 gate)" : `keep only super-enhancers at permutation FDR ${l}${kind === "lines" ? " vs all lines" : ""}`}">${l}</button>`).join("");
    U.el("atlas-fdr").querySelectorAll("button").forEach(b =>
      b.onclick = () => { fdrMax = +b.dataset.fdr; renderFdrControl(); renderTable(); });
  }

  let combo = null;
  function renderGroupPicker() {
    const gs = groupsFor(level), m = Object.fromEntries(Object.keys(manifest.levels[level].groups).map(g => [g, ginfo(level, g)]));
    if (!combo) combo = Combo.make(U.el("atlas-group"), key => { group = key; renderTable(); });
    // line level: sort and show by the display name (NIH:OVCAR-3), keyed by DepMap's stripped name
    const nm = g => m[g].label || g;
    const order = level === "line" ? gs.slice().sort((a, b) => nm(a).localeCompare(nm(b))) : gs;
    combo.setOptions(order.map(g => {
      const info = m[g], li = level === "line" && lineIdx ? lineIdx[lineKey(g)] : null;
      const tail = li ? (li.tested === false ? " · one study: ranking only" : ` · ${(li.n || {}).all ?? 0} specific vs all`)
                      : (info.n_calls != null ? ` · ${info.n_calls} calls` : "");
      const lbl = li ? `${nm(g)}${tail}` : `${nm(g)} (n=${info.n_lines}${tail})`;           // a line is n=1 by definition
      return { key: g, label: lbl, search: `${nm(g)} | ${g} | ${info.search || ""}` };
    }));
    U.el("atlas-group").placeholder = `${order.length} ${levelDef(level).label.toLowerCase()} group${order.length === 1 ? "" : "s"}; type to filter`;
    // lists are alphabetical, but a fresh view opens on the group with the most calls (not "Biliary Tract")
    if (!group || !gs.includes(group))
      group = level === "line" ? order[0] : gs.slice().sort((a, b) => (m[b].n_calls || 0) - (m[a].n_calls || 0))[0];
    combo.setValue(group);
  }

  function renderVariant() {
    const vs = manifest.variants || [], sel = U.el("atlas-variant"), wrap = U.el("atlas-variant-wrap");
    const isCalls = levelDef(level).kind === "calls";
    wrap.style.display = vs.length > 1 ? "inline-flex" : "none";
    // off the call levels the view exists for the default analysis only: show it, disabled; the choice is kept
    const shown = isCalls ? U.getVariant() : "main";
    sel.disabled = !isCalls;
    wrap.title = isCalls ? "which scoring run the lineage and disease calls come from"
      : `the ${levelDef(level).label.toLowerCase()} view is computed for the default analysis only`;
    sel.innerHTML = vs.map(v => `<option value="${v.key}"${v.key === shown ? " selected" : ""}>${U.esc(v.label)}</option>`).join("");
    // the description sits in a tooltip (and a chip in the summary bar), so switching never moves the page
    const v = vs.find(x => x.key === shown) || vs[0];
    U.el("atlas-variant-info").title = v ? `${v.label}: ${v.desc || ""}` : "";
  }
  // a short "analysis: ..." chip for the summary bar when a non-default run is on show
  const variantChip = () => {
    const v = variantOf(level) !== "main" && vinfo();
    return v ? `<span class="sep">·</span><span class="variant-chip" title="${U.esc(v.label)}: ${U.esc(v.desc || "")}">analysis: ${U.esc(v.label)}</span>` : "";
  };

  async function onLevel() {
    renderControls();
    renderVariant();
    renderScope();
    renderFdrControl();
    await loadLevel(level);
    renderGroupPicker();
    renderTable();
  }

  function sortRows(r) {
    const s = [...r].sort((a, b) => {
      let x = a[sortKey], y = b[sortKey];
      if (x == null || x === "") x = Infinity;
      if (y == null || y === "") y = Infinity;
      if (typeof x === "string") { x = x.toLowerCase(); y = String(y).toLowerCase(); }
      return x < y ? -1 : x > y ? 1 : 0;
    });
    return sortAsc ? s : s.reverse();
  }

  function currentRows() {
    const q = geneQuery.trim().toLowerCase();
    return cache[ck(level)].filter(r => r.group === group
      && (!q || String(r.gene || "").toLowerCase().includes(q))
      && (fdrMax >= 1 || (r.fdr != null && r.fdr <= fdrMax)));
  }

  // Mark each locus that overlaps a HIGHER-RANKED locus on the same chromosome, i.e. it tiles a
  // super-enhancer domain already listed above (kept a separate catalogue entry because reciprocal
  // overlap < 25%, aggregate.py). Computed from rank over the whole group so the "tiles #N" reference is
  // stable regardless of the current column sort. Returns {se_id -> rank of the domain it tiles}.
  function tilesOf(groupRows) {
    const byRank = [...groupRows].sort((a, b) => a.rank - b.rank);
    const cover = {};                                    // chrom -> [{start,end,rank}] already seen (higher rank)
    const out = {};
    for (const r of byRank) {
      if (!r.chrom) continue;
      const seen = cover[r.chrom] || (cover[r.chrom] = []);
      const hit = seen.find(c => r.start < c.end && r.end > c.start);   // overlaps an earlier (better-ranked) locus
      if (hit) out[r.se] = hit.rank;
      seen.push({ start: r.start, end: r.end, rank: r.rank });
    }
    return out;
  }

  // ---------------------------------------------------------------- cell-line level: four comparisons
  function cmpTitle(c, k) {
    if (k === "all") return `Specific vs. all ${Math.max(0, c.n_lines - 1)} other lines`;
    const w = CMPS.find(x => x[0] === k)[2];
    if (!c.stratum) return `Specific vs. lines in the same ${w}`;
    const n = Math.max(0, c.n_lines - 1);
    return `Specific vs. ${n} other line${n === 1 ? "" : "s"} in the same ${w} (${c.stratum})`;
  }
  const FLAG_TIP = { "chrY": "on chrY: presence depends on the line's sex, not its lineage (removed in the next release)",
                     "CN<0.3": "copy number below 0.3: dividing by a near-zero copy number inflates noise in deep deletions (removed in the next release)" };

  function crumbsHtml(lv, g) {
    const anchor = linesOf(lv, g);
    return ORDER.slice(0, ORDER.indexOf(lv)).map(u => {
      const p = unique(anchor, u);
      return p ? `<a class="crumb" href="#atlas" data-lv="${u}" data-g="${U.esc(p)}" title="open ${U.esc(p)} at ${levelDef(u).label.toLowerCase()} level">${U.esc(lab(u, p))}</a><span class="sep">›</span>` : "";
    }).join("");
  }

  function bindSort() {
    document.querySelectorAll("#atlas-table th[data-sort]").forEach(th =>
      th.onclick = () => {
        if (sortKey === th.dataset.sort) sortAsc = !sortAsc;
        else { sortKey = th.dataset.sort; sortAsc = !["cn_mean", "len", "n_called", "pass_ord"].includes(th.dataset.sort); }
        renderTable();
      });
  }
  function setHead(html) {
    const th = document.querySelector("#atlas-table thead");
    if (th.innerHTML !== html) { th.innerHTML = html; bindSort(); }
  }

  async function renderLine() {
    const g = group, key = lineKey(g);
    U.el("atlas-warn").style.display = "none";
    if (!["rank", "gene", "jsd", "cn_mean", "n_called", "len", "pass_ord"].includes(sortKey) && !sortKey.startsWith("fdr_")) { sortKey = "rank"; sortAsc = true; }
    if (!key) {
      U.el("atlas-desc").innerHTML = `<b>${U.esc(lab("line", g))}</b><span class="sep">·</span>no per-line data`;
      U.el("atlas-body").innerHTML = `<tr><td colspan="10" class="empty">no per-line data for this line</td></tr>`;
      rows = []; return;
    }
    if (!lineData || lineData.key !== key) {
      const d = await DataLoader.loadJSON(`data/lines/${key}.json`);
      if (group !== g || level !== "line") return;                   // the user moved on while it loaded
      lineData = d; passSel = null;                                  // a new line shows every combination
    }
    const L = lineData, C = L.comparisons, c0 = C.all, N = L.n_experiments;
    const COL = { all: "#e08214", lineage: "#7b3294", disease: "#c0392b", subtype: "#1b7837" };
    const tip = k => { const c = C[k];
      return `${cmpTitle(c, k)}${c.testable ? ` · ${c.n.toLocaleString()} called` : ` · not tested: ${c.reason || ""}`}${c.same_as ? ` · the same lines as the ${CMPS.find(x => x[0] === c.same_as)[2]} comparison` : ""}`; };
    setHead(`<tr>
      <th data-sort="rank" title="specificity rank (1 = most specific); one score, the same for all four comparisons. Click to sort.">Rank</th>
      <th data-sort="gene" title="nearest protein-coding gene. Click to sort.">Nearest gene</th>
      <th data-sort="pass_ord" title="which comparisons call the SE: A all lines, L same lineage, D same disease, S same subtype. Click to sort (broadest first); ▾ to filter.">Called vs ${U.passHeadBtn("atlas", passSel)}</th>
      <th data-sort="jsd" title="CaCTS score = Jensen–Shannon divergence. Lower = more specific. Click to sort.">JSD</th>
      ${CMPS.map(([k, w]) => `<th data-sort="fdr_${k}" style="color:${COL[k]}" title="${U.esc(tip(k))}. Click to sort (not-called rows last).">vs. ${w} FDR</th>`).join("")}
      <th data-sort="cn_mean" title="copy-number ratio at the SE in this line (1 = neutral; > 1.3 amplified). Click to sort.">Copy no.</th>
      <th title="genomic locus (GRCh38); links to UCSC">Locus</th>
      <th data-sort="n_called" title="experiments of this line whose SE calls cover the locus. Click to sort.">Called</th></tr>`);
    document.querySelectorAll("#atlas-table th[data-sort]").forEach(th => {
      th.classList.toggle("sorted-asc", th.dataset.sort === sortKey && sortAsc);
      th.classList.toggle("sorted-desc", th.dataset.sort === sortKey && !sortAsc);
    });

    const all = (L.rows || []).map(a => {
      const o = Object.fromEntries(L.cols.map((col, i) => [col, a[i]]));
      return { ...o, cn_mean: o.cn, n_called: o.called, len: (o.end || 0) - (o.start || 0), group: g, pass_ord: U.passOrd(o.pass) };
    });
    const q = geneQuery.trim().toLowerCase();
    rows = sortRows(all.filter(r => (!q || String(r.gene || "").toLowerCase().includes(q)) && (!passSel || passSel.has(r.pass || ""))
                                  && (fdrMax >= 1 || (r.fdr_all != null && r.fdr_all <= fdrMax))));
    const sets = CMPS.slice(1).map(([k, w, lvl]) => { const c = C[k];
      return c.testable ? `<span title="${U.esc(tip(k))}">${w.toLowerCase()} <b>${c.n.toLocaleString()}</b> <span class="muted-s">(vs ${Math.max(0, c.n_lines - 1)} in ${U.esc(c.stratum || "")})</span></span>`
                        : `<span title="${U.esc(tip(k))}">${w.toLowerCase()} <span class="muted-s">not tested</span></span>`; }).join(`<span class="sep">·</span>`);
    const status = c0.testable
      ? `<span class="mono">${c0.n.toLocaleString()}</span>&nbsp;specific vs all ${Math.max(0, c0.n_lines - 1)} other lines<span class="sep">·</span>of those, called within its ${sets}`
      : `<span class="rank-only-badge">not tested</span> <span class="cmp-status">${U.esc(c0.reason || "")}; the rows are the top of the ranking</span>`;
    U.el("atlas-desc").innerHTML = crumbsHtml("line", g) + `<b>${U.esc(L.name)}</b><span class="sep">·</span><span class="cmp-status">${status}</span>`;
    U.el("atlas-desc").querySelectorAll(".crumb").forEach(el => el.onclick = e => { e.preventDefault(); goTo(el.dataset.lv, el.dataset.g); });
    const fdrTd = (r, k) => { const c = C[k], v = r["fdr_" + k];
      if (!c.testable) return `<td class="mono num muted-s" title="not tested: ${U.esc(c.reason || "")}">–</td>`;
      return v == null ? `<td class="mono num muted-s" title="tested vs ${U.esc(c.level)}, not called">&gt; 0.10</td>`
                       : `<td class="mono num sig" title="called vs ${U.esc(c.level)} at FDR ≤ 0.10">${fmtFdrVal(v)}</td>`; };
    const tiles = tilesOf(all);
    const cnCell = v => v == null ? `<td class="mono">n/a</td>` : `<td class="mono ${U.cnClass(v)}">${(+v).toFixed(2)}</td>`;
    U.el("atlas-body").innerHTML = rows.map(r => `
      <tr>
        <td class="mono">${r.rank}</td>
        <td class="gene"><a class="gv-link" href="#line" data-locus="${r.chrom}:${Math.max(1, r.start - 25000)}-${r.end + 25000}" title="open ${U.esc(r.gene || r.se)} in the Genomic View (IGV) for ${U.esc(L.name)}">${U.esc(r.gene || "n/a")}</a><span class="th-sub">${r.dist_kb === 0 ? " overlaps" : " " + r.dist_kb + " kb"}</span>${r.flag ? ` <span class="flag-chip" title="possible artifact: ${U.esc(FLAG_TIP[r.flag] || r.flag)}">⚠ ${U.esc(r.flag)}</span>` : ""}</td>
        <td>${U.passBadge(r.pass)}</td>
        <td class="mono num">${U.fmtJsd(r.jsd)}</td>
        ${CMPS.map(([k]) => fdrTd(r, k)).join("")}
        ${cnCell(r.cn_mean)}
        <td class="coord"><a href="${U.ucsc(r.chrom, r.start, r.end)}" target="_blank" rel="noopener" title="${r.se}: open in the UCSC genome browser (GRCh38)">${r.chrom}:${(+r.start).toLocaleString()}–${(+r.end).toLocaleString()}</a>${tiles[r.se] != null
          ? `<span class="tiles-chip" title="same super-enhancer domain as row #${tiles[r.se]} above">↳ tiles #${tiles[r.se]}</span>` : ""}</td>
        <td class="mono num" title="experiments of ${U.esc(L.name)} whose SE calls cover this locus (of ${N}); 0 = specific signal, but not an SE in this line">${r.n_called}/${N}</td>
      </tr>`).join("") || `<tr><td colspan="11" class="empty">${!c0.testable && fdrMax < 1 ? "this line was not tested, so there is no FDR to filter on" : c0.testable ? (passSel || fdrMax < 1 ? "no rows pass the current filters" : "no specific super-enhancers") : "no rankings available"}</td></tr>`;
    U.el("atlas-body").querySelectorAll("a.gv-link").forEach(a => a.onclick = e => { e.preventDefault(); LineView.goto(key, a.dataset.locus); });
    U.wirePassHead("atlas", all, passSel, s => { passSel = s; renderLine(); });
    U.el("atlas-foot").innerHTML = (c0.testable
      ? `Top ${rows.length.toLocaleString()} of ${c0.n.toLocaleString()} super-enhancers specific vs all lines (permutation FDR ≤ 0.10), ranked by JSD. `
      : `Not tested (${U.esc(c0.reason || "")}); the rows are the top of the ranking, not calls. `)
      + `All four comparisons share this score and differ only in the null: all lines, or only the other lines of the same lineage, primary disease or subtype. Each FDR column says whether the super-enhancer also stands out there (<b>&gt; 0.10</b> = tested, not called; <b>–</b> = not tested; hover a header for the comparison set). Every relatives call is also a vs-all call, so this list holds them all. `
      + `<b>Called vs</b> repeats the called comparisons (${U.passBadge("ALDS")}); click it to sort, or ▾ to filter by combination: a lineage programme reads <b>A·L</b> without <b>D</b>. `
      + `A comparison is tested only with at least two independent studies and at least 4 lines in the group. <b>Called</b> = experiments of this line whose SE calls cover the locus. <span class="flag-chip" style="margin:0">⚠</span> marks known artifact classes (chrY; copy number below 0.3). Click a gene to open it in the Genomic View (IGV).`;
  }

  function renderTable() {
    if (level === "line") return renderLine();
    U.el("atlas-cmp").style.display = "none"; U.el("atlas-cmp").innerHTML = "";
    U.closePassHead();
    if (baseHead) setHead(baseHead);
    const def = levelDef(level), info = ginfo(level, group);
    const isCalls = def.kind === "calls";
    const groupAll = cache[ck(level)].filter(r => r.group === group);
    const total = groupAll.length;
    const tiles = tilesOf(groupAll);                     // se_id -> rank of the domain it tiles
    rows = sortRows(currentRows());

    // breadcrumb: the group's parents at the levels above, each a link to that group
    const anchor = linesOf(level, group);
    const crumbs = ORDER.slice(0, ORDER.indexOf(level)).map(u => {
      const p = unique(anchor, u);
      return p ? `<a class="crumb" href="#atlas" data-lv="${u}" data-g="${U.esc(p)}" title="open ${U.esc(p)} at ${levelDef(u).label.toLowerCase()} level">${U.esc(lab(u, p))}</a><span class="sep">›</span>` : "";
    }).join("");
    U.el("atlas-desc").innerHTML = crumbs +
      `<b>${U.esc(info.label || group)}</b><span class="sep">·</span><span class="mono">${info.n_lines} cell line${info.n_lines > 1 ? "s" : ""}</span>` +
      (isCalls ? `<span class="sep">·</span><span class="mono" title="catalogue entries; nested / tiling loci mean fewer independent SE domains, see the ↳ markers and About & methods">${info.n_calls.toLocaleString()} specific SEs</span><span>&nbsp;at permutation FDR ≤ 0.10</span>`
               : `<span class="rank-only-badge">rankings only</span>`) + variantChip();
    U.el("atlas-desc").querySelectorAll(".crumb").forEach(c =>
      c.onclick = e => { e.preventDefault(); goTo(c.dataset.lv, c.dataset.g); });
    U.el("atlas-warn").style.display = isCalls ? "none" : "block";

    const topN = 100;
    const shown = rows.slice(0, topN);
    const cnCell = v => {
      if (v == null) return `<td class="mono">n/a</td>`;
      const amp = +v > 1.3;
      return `<td class="mono ${U.cnClass(v)}"${amp ? ` title="amplified in this group's cell lines (mean copy-number ratio ${(+v).toFixed(2)}× vs neutral 1×), check the CN-ablation tab"` : ""}>${(+v).toFixed(2)}</td>`;
    };
    const concBadge = r => (r.conc === 1)
      ? ` <span class="conc-badge" title="Cross-layer support: ${U.esc(r.gene)} is itself a group-specific gene here (DepMap expression), and its expression tracks this SE's H3K27ac across cell lines (Spearman ρ = ${r.rho}). See the Concordance tab.">⇌ ${r.rho == null || r.rho === "" ? "" : (+r.rho).toFixed(2)}</span>`
      : (r.conc === 0 ? ` <span class="conc-no" title="the nearest gene is not itself group-specific in the expression layer">·</span>` : "");
    U.el("atlas-body").innerHTML = shown.map(r => `
      <tr>
        <td class="mono" title="specificity rank in this group (1 = most specific)">${r.rank}</td>
        <td class="gene">${U.esc(r.gene || "n/a")}<span class="th-sub" title="distance from the SE to this gene's body">${r.dist_kb === 0 ? " overlaps" : " " + r.dist_kb + " kb"}</span>${concBadge(r)}</td>
        <td class="mono num" title="Jensen–Shannon divergence specificity score (lower = more specific)">${U.fmtJsd(r.jsd)}</td>
        <td class="mono num${isCalls && r.fdr <= 0.10 ? " sig" : ""}" title="${!isCalls ? "FDR shown for reference, not a call at this resolution"
          : (r.fdr <= 0.10 ? "passes the permutation FDR ≤ 0.10 call" : "above the FDR ≤ 0.10 threshold")}">${fmtFdrVal(r.fdr)}</td>
        ${cnCell(r.cn_mean)}
        <td class="coord">${r.chrom ? `<a href="${U.ucsc(r.chrom, r.start, r.end)}" target="_blank" rel="noopener" title="${r.se}: open ${r.chrom}:${(+r.start).toLocaleString()}–${(+r.end).toLocaleString()} in the UCSC genome browser (GRCh38)">${r.chrom}:${(+r.start).toLocaleString()}–${(+r.end).toLocaleString()}</a>` : ", "}${tiles[r.se] != null
          ? `<span class="tiles-chip" title="Same super-enhancer domain as row #${tiles[r.se]} above, this locus overlaps it, but was kept a separate catalogue entry because their reciprocal overlap is under 25% (the cSEAdb-style union). Count SE domains, not rows.">↳ tiles #${tiles[r.se]}</span>` : ""}</td>
        <td class="mono num" title="length of the super-enhancer locus (end − start). Nested / overlapping entries near one gene tile a single SE domain.">${fmtLen(r.len)}</td>
        <td class="mono num" title="number of experiments in which this exact locus was itself called a super-enhancer, a low value (e.g. 1 beside a 74) marks a single-sample sub-peak of the broader domain above">${r.n_called == null ? "" : r.n_called}</td>
      </tr>`).join("") ||
      `<tr><td colspan="8" class="empty">no super-enhancers match the current filter</td></tr>`;

    const filtered = geneQuery.trim() || fdrMax < 1;
    const nestNote = ` Distinct rows near one gene often tile a single SE domain (<span class="tiles-chip" style="margin:0">↳</span>), compare the Locus and Length; a low <b>Called</b> count marks a single-sample sub-peak.`;
    U.el("atlas-foot").innerHTML = isCalls
      ? `Showing ${Math.min(topN, rows.length).toLocaleString()} of ${rows.length.toLocaleString()}${filtered ? ` matched (of <b>${total.toLocaleString()}</b> specific SEs in this group)` : ` <b>specific SEs</b> (permutation FDR ≤ 0.10)`}, ranked by JSD. Bold FDR passes the call threshold; the <span class="conc-badge" style="margin:0">⇌</span> marks a gene that is itself group-specific in expression.${nestNote} Download the current view below.`
      : `Top ${shown.length} SEs by specificity ranking${filtered ? ` (of ${rows.length.toLocaleString()} matched)` : ""}. <b>These are rankings, not calls</b>: at this resolution most groups have too few cell lines for the permutation null to support a significance threshold (see the note above).${nestNote}`;
    document.querySelectorAll("#atlas-table th[data-sort]").forEach(th => {
      th.classList.toggle("sorted-asc", th.dataset.sort === sortKey && sortAsc);
      th.classList.toggle("sorted-desc", th.dataset.sort === sortKey && !sortAsc);
    });
  }

  async function init() {
    manifest = await DataLoader.loadJSON("data/manifest.json");
    hier = (manifest.hierarchy || []).map(([line, lineage, disease, subtype]) => ({ line, lineage, disease, subtype }));
    baseHead = document.querySelector("#atlas-table thead").innerHTML;
    bindSort();
    U.el("atlas-variant").onchange = e => { U.setVariant(e.target.value); onLevel(); };
    const filterInput = U.el("atlas-filter");
    filterInput.addEventListener("input", () => { geneQuery = filterInput.value; renderTable(); });
    U.el("atlas-dl").onclick = () => U.downloadTSV(`SE-CaCTS.${level}.${group}${level !== "line" && variantOf(level) !== "main" ? "." + variantOf(level) : ""}.tsv`, [
      { label: "rank", key: "rank" }, { label: "se", key: "se" }, { label: "nearest_gene", key: "gene" },
      { label: "dist_kb", key: "dist_kb" }, { label: "jsd", key: "jsd" },
      ...(level === "line" ? [{ label: "fdr_vs_all", key: "fdr_all" }, { label: "fdr_vs_lineage", key: "fdr_lineage" },
                              { label: "fdr_vs_disease", key: "fdr_disease" }, { label: "fdr_vs_subtype", key: "fdr_subtype" },
                              { label: "called_in", key: "pass" }]
                           : [{ label: "fdr_permutation", key: "fdr" }]),
      { label: "cn_mean", key: "cn_mean" }, { label: "gene_concordant", key: "conc" },
      { label: "gene_expr_rho", key: "rho" }, { label: "chrom", key: "chrom" },
      { label: "start", key: "start" }, { label: "end", key: "end" },
      { label: "length_bp", key: "len" }, { label: "n_called_as_SE", key: "n_called" },
    ], rows);
    await onLevel();
  }
  return { init };
})();
