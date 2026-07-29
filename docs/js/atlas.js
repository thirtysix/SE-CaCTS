/* atlas.js — browse specific super-enhancers per group, at each resolution.
   Enforces the resolution rule: CALLS at lineage/disease, RANKINGS ONLY at subtype/line. */
const Atlas = (() => {
  let manifest, level = "lineage", group = null, rows = [], cache = {};
  let sortKey = "rank", sortAsc = true, geneQuery = "", fdrMax = 1;

  const levelDef = k => U.LEVELS.find(l => l.key === k);
  const fmtLen = bp => bp == null ? "" : (bp < 1000 ? `${bp} bp` : `${(bp / 1000).toFixed(1)} kb`);
  // The staged FDR is rounded to 4 dp, so report anything below 0.001 as a bound rather than as a
  // spuriously precise number (at cell-line level the degenerate null rounds many rows to 0).
  const fmtFdrVal = v => {
    if (v == null || v === "" || isNaN(+v)) return "—";
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
    const g = manifest.levels[k].groups;
    return Object.keys(g).sort((a, b) => (g[b].n_calls || 0) - (g[a].n_calls || 0) || a.localeCompare(b));
  }

  function renderControls() {
    U.el("atlas-level").innerHTML = U.LEVELS.map(l =>
      `<button data-lv="${l.key}" aria-selected="${l.key === level}" title="${l.kind === "calls"
        ? "specificity calls are supported at this resolution" : "rankings only — the panel does not support calls here (few lines per group)"}">${l.label}</button>`).join("");
    U.el("atlas-level").querySelectorAll("button").forEach(b =>
      b.onclick = () => { if (b.dataset.lv !== level) { level = b.dataset.lv; group = null; onLevel(); } });
  }

  function renderFdrControl() {
    // FDR threshold filter — only meaningful where calls exist (gotcha 72). Hidden at rankings-only levels.
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

  function renderGroupPicker() {
    const gs = groupsFor(level), m = manifest.levels[level].groups;
    const sel = U.el("atlas-group");
    sel.innerHTML = gs.map(g => {
      const info = m[g];
      const tail = info.n_calls != null ? ` — ${info.n_calls} calls` : "";
      return `<option value="${U.esc(g)}">${U.esc(g)} (n=${info.n_lines}${tail})</option>`;
    }).join("");
    if (!group || !gs.includes(group)) group = gs[0];
    sel.value = group;
    sel.onchange = () => { group = sel.value; renderTable(); };
  }

  async function onLevel() {
    renderControls();
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

  // Mark each locus that overlaps a HIGHER-RANKED locus on the same chromosome — i.e. it tiles a
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

    U.el("atlas-desc").innerHTML =
      `<b>${U.esc(group)}</b><span class="sep">·</span><span class="mono">${info.n_lines} cell line${info.n_lines > 1 ? "s" : ""}</span>` +
      (isCalls ? `<span class="sep">·</span><span class="mono" title="catalogue entries; nested / tiling loci mean fewer independent SE domains — see the ↳ markers and About & methods">${info.n_calls.toLocaleString()} specific SEs</span><span>&nbsp;at permutation FDR ≤ 0.10</span>`
               : `<span class="rank-only-badge">rankings only</span>`);
    U.el("atlas-warn").style.display = isCalls ? "none" : "block";

    const topN = 100;
    const shown = rows.slice(0, topN);
    const cnCell = v => {
      if (v == null) return `<td class="mono">—</td>`;
      const amp = +v > 1.3;
      return `<td class="mono ${U.cnClass(v)}"${amp ? ` title="amplified in this group's cell lines (mean copy-number ratio ${(+v).toFixed(2)}× vs neutral 1×) — check the CN-ablation tab"` : ""}>${(+v).toFixed(2)}</td>`;
    };
    const concBadge = r => (r.conc === 1)
      ? ` <span class="conc-badge" title="Cross-layer support: ${U.esc(r.gene)} is itself a group-specific gene here (DepMap expression), and its expression tracks this SE's H3K27ac across cell lines (Spearman ρ = ${r.rho}). See the Concordance tab.">⇌ ${r.rho == null || r.rho === "" ? "" : (+r.rho).toFixed(2)}</span>`
      : (r.conc === 0 ? ` <span class="conc-no" title="the nearest gene is not itself group-specific in the expression layer">·</span>` : "");
    U.el("atlas-body").innerHTML = shown.map(r => `
      <tr>
        <td class="mono" title="specificity rank in this group (1 = most specific)">${r.rank}</td>
        <td class="gene">${U.esc(r.gene || "—")}<span class="th-sub" title="distance from the SE to this gene's body">${r.dist_kb === 0 ? " overlaps" : " " + r.dist_kb + " kb"}</span>${concBadge(r)}</td>
        <td class="mono num" title="Jensen–Shannon divergence specificity score (lower = more specific)">${U.fmtJsd(r.jsd)}</td>
        <td class="mono num${isCalls && r.fdr <= 0.10 ? " sig" : ""}" title="${!isCalls ? "FDR shown for reference — not a call at this resolution"
          : (r.fdr <= 0.10 ? "passes the permutation FDR ≤ 0.10 call" : "above the FDR ≤ 0.10 threshold")}">${fmtFdrVal(r.fdr)}</td>
        ${cnCell(r.cn_mean)}
        <td class="coord">${r.chrom ? `<a href="${U.ucsc(r.chrom, r.start, r.end)}" target="_blank" rel="noopener" title="${r.se} — open ${r.chrom}:${(+r.start).toLocaleString()}–${(+r.end).toLocaleString()} in the UCSC genome browser (GRCh38)">${r.chrom}:${(+r.start).toLocaleString()}–${(+r.end).toLocaleString()}</a>` : "—"}${tiles[r.se] != null
          ? `<span class="tiles-chip" title="Same super-enhancer domain as row #${tiles[r.se]} above — this locus overlaps it, but was kept a separate catalogue entry because their reciprocal overlap is under 25% (the cSEAdb-style union). Count SE domains, not rows.">↳ tiles #${tiles[r.se]}</span>` : ""}</td>
        <td class="mono num" title="length of the super-enhancer locus (end − start). Nested / overlapping entries near one gene tile a single SE domain.">${fmtLen(r.len)}</td>
        <td class="mono num" title="number of experiments in which this exact locus was itself called a super-enhancer — a low value (e.g. 1 beside a 74) marks a single-sample sub-peak of the broader domain above">${r.n_called == null ? "" : r.n_called}</td>
      </tr>`).join("") ||
      `<tr><td colspan="8" class="empty">no super-enhancers match the current filter</td></tr>`;

    const filtered = geneQuery.trim() || fdrMax < 1;
    const nestNote = ` Distinct rows near one gene often tile a single SE domain (<span class="tiles-chip" style="margin:0">↳</span>) — compare the Locus and Length; a low <b>Called</b> count marks a single-sample sub-peak.`;
    U.el("atlas-foot").innerHTML = isCalls
      ? `Showing ${Math.min(topN, rows.length).toLocaleString()} of ${rows.length.toLocaleString()}${filtered ? ` matched (of <b>${total.toLocaleString()}</b> specific SEs in this group)` : ` <b>specific SEs</b> (permutation FDR ≤ 0.10)`}, ranked by JSD. Bold FDR passes the call threshold; the <span class="conc-badge" style="margin:0">⇌</span> marks a gene that is itself group-specific in expression.${nestNote} Download the current view below.`
      : `Top ${shown.length} SEs by specificity ranking${filtered ? ` (of ${rows.length.toLocaleString()} matched)` : ""}. <b>These are rankings, not calls</b> — at this resolution most groups have too few cell lines for the permutation null to support a significance threshold (see the note above).${nestNote}`;
    document.querySelectorAll("#atlas-table th[data-sort]").forEach(th => {
      th.classList.toggle("sorted-asc", th.dataset.sort === sortKey && sortAsc);
      th.classList.toggle("sorted-desc", th.dataset.sort === sortKey && !sortAsc);
    });
  }

  async function init() {
    manifest = await DataLoader.loadJSON("data/manifest.json");
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
