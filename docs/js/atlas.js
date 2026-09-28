/* atlas.js, browse specific super-enhancers per group, at each resolution.
   Enforces the resolution rule: CALLS at lineage/disease, RANKINGS ONLY at subtype/line. */
const Atlas = (() => {
  let manifest, level = "lineage", group = null, rows = [], cache = {};
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
  const fmtLen = bp => bp == null ? "" : (bp < 1000 ? `${bp} bp` : `${(bp / 1000).toFixed(1)} kb`);
  // The staged FDR is rounded to 4 dp, so report anything below 0.001 as a bound rather than as a
  // spuriously precise number (at cell-line level the degenerate null rounds many rows to 0).
  const fmtFdrVal = v => {
    if (v == null || v === "" || isNaN(+v)) return "n/a";
    v = +v;
    return v < 0.001 ? "<0.001" : v.toFixed(3);
  };

  async function loadLevel(k) {
    if (!cache[k]) {
      const r = (await DataLoader.loadTSV(levelDef(k).file)).rows;
      r.forEach(x => { x.len = (x.end || 0) - (x.start || 0); });   // SE span, for the Length column + sort
      cache[k] = r;
    }
    return cache[k];
  }

  function groupsFor(k) {
    const g = manifest.levels[k].groups, i = ORDER.indexOf(k);
    const allowed = useScope && hier.length && i > 0 ? new Set(inScope(scope, i).map(h => h[k])) : null;
    return Object.keys(g).filter(x => !allowed || allowed.has(x))
      .sort((a, b) => (g[b].n_calls || 0) - (g[a].n_calls || 0) || a.localeCompare(b));
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
        ? "specificity calls are supported at this resolution" : "rankings only: the panel does not support calls here (few lines per group)"}">${l.label}</button>`).join("");
    U.el("atlas-level").querySelectorAll("button").forEach(b =>
      b.onclick = () => { if (b.dataset.lv !== level) switchLevel(b.dataset.lv); });
  }

  function renderFdrControl() {
    // FDR threshold filter, only meaningful where calls exist (gotcha 72). Hidden at rankings-only levels.
    const isCalls = levelDef(level).kind === "calls";
    U.el("atlas-fdr-wrap").style.display = isCalls ? "flex" : "none";
    if (!isCalls) { fdrMax = 1; return; }
    const opts = [[1, "All"], [0.10, "≤ 0.10"], [0.05, "≤ 0.05"], [0.01, "≤ 0.01"]];
    U.el("atlas-fdr").innerHTML = opts.map(([v, l]) =>
      `<button data-fdr="${v}" aria-selected="${v === fdrMax}" title="${v === 1
        ? "show every call (all pass the ≤ 0.10 gate)" : `keep only super-enhancers at permutation FDR ${l}`}">${l}</button>`).join("");
    U.el("atlas-fdr").querySelectorAll("button").forEach(b =>
      b.onclick = () => { fdrMax = +b.dataset.fdr; renderFdrControl(); renderTable(); });
  }

  let combo = null;
  function renderGroupPicker() {
    const gs = groupsFor(level), m = manifest.levels[level].groups;
    if (!combo) combo = Combo.make(U.el("atlas-group"), key => { group = key; renderTable(); });
    // line level: sort and show by the display name (NIH:OVCAR-3), keyed by DepMap's stripped name
    const nm = g => m[g].label || g;
    const order = level === "line" ? gs.slice().sort((a, b) => nm(a).localeCompare(nm(b))) : gs;
    combo.setOptions(order.map(g => {
      const info = m[g];
      const tail = info.n_calls != null ? ` · ${info.n_calls} calls` : "";
      return { key: g, label: `${nm(g)} (n=${info.n_lines}${tail})`, search: `${nm(g)} | ${g} | ${info.search || ""}` };
    }));
    U.el("atlas-group").placeholder = `${order.length} ${levelDef(level).label.toLowerCase()} group${order.length === 1 ? "" : "s"}; type to filter`;
    if (!group || !gs.includes(group)) group = order[0];
    combo.setValue(group);
  }

  async function onLevel() {
    renderControls();
    renderScope();
    renderFdrControl();
    await loadLevel(level);
    renderGroupPicker();
    renderTable();
  }

  function sortRows(r) {
    const s = [...r].sort((a, b) => {
      let x = a[sortKey], y = b[sortKey];
      if (typeof x === "string") { x = x.toLowerCase(); y = (y || "").toLowerCase(); }
      return x < y ? -1 : x > y ? 1 : 0;
    });
    return sortAsc ? s : s.reverse();
  }

  function currentRows() {
    const q = geneQuery.trim().toLowerCase();
    return cache[level].filter(r => r.group === group
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

  function renderTable() {
    const def = levelDef(level), info = manifest.levels[level].groups[group];
    const isCalls = def.kind === "calls";
    const groupAll = cache[level].filter(r => r.group === group);
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
               : `<span class="rank-only-badge">rankings only</span>`);
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
    document.querySelectorAll("#atlas-table th[data-sort]").forEach(th =>
      th.onclick = () => {
        if (sortKey === th.dataset.sort) sortAsc = !sortAsc;
        else { sortKey = th.dataset.sort; sortAsc = !["cn_mean", "len", "n_called"].includes(th.dataset.sort); }
        renderTable();
      });
    const filterInput = U.el("atlas-filter");
    filterInput.addEventListener("input", () => { geneQuery = filterInput.value; renderTable(); });
    U.el("atlas-dl").onclick = () => U.downloadTSV(`SE-CaCTS.${level}.${group}.tsv`, [
      { label: "rank", key: "rank" }, { label: "se", key: "se" }, { label: "nearest_gene", key: "gene" },
      { label: "dist_kb", key: "dist_kb" }, { label: "jsd", key: "jsd" },
      { label: "fdr_permutation", key: "fdr" },
      { label: "cn_mean", key: "cn_mean" }, { label: "gene_concordant", key: "conc" },
      { label: "gene_expr_rho", key: "rho" }, { label: "chrom", key: "chrom" },
      { label: "start", key: "start" }, { label: "end", key: "end" },
      { label: "length_bp", key: "len" }, { label: "n_called_as_SE", key: "n_called" },
    ], rows);
    await onLevel();
  }
  return { init };
})();
