/* lineview.js, the Genomic View: one cell line in igv.js with the SEs it calls, the SEs specific to it in four
   comparisons (vs all lines, and vs the lines of its subtype, primary disease and lineage), its copy number as
   scored, and its H3K27ac coverage streamed from ChIP-Atlas.
   Data: data/lines/ (phase2/scripts/61_stage_lines.py). Coverage is ChIP-Atlas's own per-experiment bigWig,
   read by range request (their server allows it for this site); it is raw RPM, not the S3-normalised signal
   the atlas compares, so two studies' tracks are not directly comparable in height. */
const LineView = (() => {
  const DIR = "data/lines/";
  const DEFAULT_COVERAGE = 3;                  // coverage tracks shown on open (one per study where possible)
  let index = null, combo = null, browser = null, line = null, token = 0, pending = null;
  let passSel = null, sort = { k: "rank", asc: true };        // table filter (pass patterns) and sort

  // comparison -> track file tag, colour, level word
  const CMP = [                                // broad to narrow
    ["all", "vsall", "#e08214", "all lines"],
    ["lineage", "vslin", "#7b3294", "lineage"],
    ["disease", "vsdis", "#c0392b", "primary disease"],
    ["subtype", "vssub", "#1b7837", "subtype"],
  ];
  const COLORS = { called: "#7a7a7a", cn: "#2b6cb0", cnNeg: "#9b2c2c", cov: "#3d3d3d" };
  const CN_LABEL = { depmap_wgs: "DepMap WGS", depmap_mc_wes: "DepMap WES", cmp_wes: "CMP WES", ccle_snp6: "CCLE SNP6",
                     input_inferred: "inferred from ChIP input" };

  // "Specific vs. 40 other lines in the same primary disease (Invasive Breast Carcinoma)"
  function cmpTitle(c, key) {
    if (key === "all") return `Specific vs. all ${Math.max(0, c.n_lines - 1)} other lines`;
    const lvl = CMP.find(x => x[0] === key)[3];
    if (!c.stratum) return `Specific vs. lines in the same ${lvl}`;
    const k = Math.max(0, c.n_lines - 1);
    return `Specific vs. ${k} other line${k === 1 ? "" : "s"} in the same ${lvl} (${c.stratum})`;
  }
  const rowsOf = L => (L.rows || []).map(a => { const o = Object.fromEntries(L.cols.map((c, i) => [c, a[i]])); o.pass_ord = U.passOrd(o.pass); return o; });
  function sorted(rows) {
    const v = (r, k) => { const x = r[k]; return x == null || x === "" ? Infinity : x; };
    const out = [...rows].sort((a, b) => v(a, sort.k) < v(b, sort.k) ? -1 : v(a, sort.k) > v(b, sort.k) ? 1 : 0);
    return sort.asc ? out : out.reverse();
  }
  // FDR cell for one comparison: the value when called, "> 0.10" when tested but not called, "–" when not tested
  function fdrCell(L, r, k) {
    const c = L.comparisons[k], v = r["fdr_" + k];
    if (!c.testable) return `<td class="muted-s" title="not tested: ${U.esc(c.reason || "")}">–</td>`;
    return v == null ? `<td class="muted-s" title="tested, not called at FDR ≤ 0.10">&gt; 0.10</td>`
                     : `<td class="sig" title="called at FDR ≤ 0.10">${v < 0.001 ? "<0.001" : v.toFixed(3)}</td>`;
  }

  async function init() {
    index = await DataLoader.loadJSON(DIR + "index.json");
    const opts = Object.entries(index).map(([key, v]) => ({ key, label: v.name, search: `${v.search} | ${v.lineage}` }))
      .sort((a, b) => a.label.localeCompare(b.label));
    combo = Combo.make(U.el("lv-line"), key => show(key));
    combo.setOptions(opts);
    U.el("lv-count").textContent = `${opts.length} cell lines`;
    const first = (pending && index[pending.key] && { key: pending.key }) || opts.find(o => o.label === "NIH:OVCAR-3") || opts[0];
    if (first) { combo.setValue(first.key); show(first.key); }
  }

  // open a line at a locus from another tab (the SE Atlas gene links)
  function goto(key, locusStr) {
    pending = { key, locus: locusStr };
    if (location.hash !== "#line") location.hash = "line";
    if (index && combo) { combo.setValue(key); show(key); }
  }

  const locus = r => `${r.chrom}:${Math.max(1, r.start - 25000)}-${r.end + 25000}`;
  const studiesOf = exps => [...new Set(exps.map(e => e.study).filter(Boolean))];

  function pickCoverage(exps) {                // first experiment of each study, up to DEFAULT_COVERAGE
    const seen = new Set(), out = [];
    for (const e of exps) { if (!seen.has(e.study)) { seen.add(e.study); out.push(e.srx); } if (out.length >= DEFAULT_COVERAGE) break; }
    return new Set(out);
  }

  async function show(key) {
    const my = ++token;
    const raw = await DataLoader.loadJSON(`${DIR}${key}.json`);
    if (my !== token) return;
    // a copy: DataLoader caches the parsed JSON and the SE Atlas reads the same compact rows, so converting them
    // in place broke both views on reopening a line
    line = { ...raw, rows: rowsOf(raw) }; passSel = null; sort = { k: "rank", asc: true };
    renderSummary(); renderTables(); renderExperiments();
    await renderBrowser(my);
  }

  function renderSummary() {
    const L = line, st = L.n_studies ?? studiesOf(L.experiments).length;
    U.el("lv-summary").innerHTML = `
      <div class="lv-title"><h2>${U.esc(L.name)}</h2>
        <span class="muted-s">${U.esc(L.lineage || "no lineage")} › ${U.esc(L.disease || "no primary disease")} › ${U.esc(L.subtype || "no subtype")}</span></div>
      <div class="lv-kpis">
        <div><b>${L.experiments.length}</b><span>experiments</span></div>
        <div><b>${st}</b><span>studies</span></div>
        <div><b>${L.n_called.toLocaleString()}</b><span>SEs called</span></div>
        ${CMP.map(([k, , col, w]) => { const c = L.comparisons[k];
          return `<div title="${U.esc(cmpTitle(c, k))}"><b style="color:${col}">${c.testable ? c.n.toLocaleString() : "–"}</b><span>specific vs ${k === "all" ? "all" : "same " + w}</span></div>`; }).join("")}
      </div>
      <p class="muted-s">Copy number: ${U.esc(CN_LABEL[L.cn_source] || L.cn_source)}${L.has_cn ? "" : " (no track)"}.
        ${st < 2 ? " <b>One study only:</b> specific-SE calls need two independent studies, so the lists below are rankings, not calls." : ""}</p>`;
  }

  const FLAG_TIP = { "chrY": "on chrY: presence depends on the line's sex, not its lineage (removed in the next release)",
                     "CN<0.3": "copy number below 0.3 here: dividing by a near-zero copy number inflates noise in deep deletions (removed in the next release)" };
  const flagChip = f => f ? ` <span class="flag-chip" title="possible artifact: ${U.esc(FLAG_TIP[f] || f)}">⚠ ${U.esc(f)}</span>` : "";

  function renderTables() {
    const L = line, c0 = L.comparisons.all, nAll = c0.testable ? c0.n : 0;
    const rows = sorted(L.rows.filter(r => !passSel || passSel.has(r.pass || "")));
    const arrow = k => sort.k === k ? (sort.asc ? " ▲" : " ▼") : "";
    const head = CMP.map(([k, , col, w]) => { const c = L.comparisons[k];
      return `<th class="th-btn" data-k="fdr_${k}" style="color:${col}" title="${U.esc(cmpTitle(c, k))}${c.testable ? ` · ${c.n.toLocaleString()} called` : " · not tested: " + U.esc(c.reason || "")}${c.same_as ? " · the same lines as the " + U.esc(CMP.find(x => x[0] === c.same_as)[3]) + " comparison" : ""}">vs. ${k === "all" ? "All" : w === "primary disease" ? "Disease" : w[0].toUpperCase() + w.slice(1)} FDR${arrow("fdr_" + k)}</th>`; }).join("");
    const count = !c0.testable ? `not tested: ${U.esc(c0.reason || "")}; the rows are the top of the ranking`
      : (nAll > L.rows.length ? `top ${L.rows.length} of ${nAll.toLocaleString()} specific vs all lines` : `${nAll.toLocaleString()} specific vs all lines`)
        + (passSel ? ` · ${rows.length} shown` : "");
    const body = !rows.length ? `<p class="muted-s">${passSel ? "No rows with the selected combinations." : "None."}</p>` :
      `<table class="tbl lv-tbl" id="lv-t"><thead><tr><th class="th-btn" data-k="rank" title="specificity rank (1 = most specific); the same score for every comparison. Click to sort.">#${arrow("rank")}</th><th>Nearest gene</th><th class="th-btn" data-k="pass_ord" title="which comparisons call the SE: A all lines, L same lineage, D same disease, S same subtype. Click to sort (broadest first); ▾ to filter.">Called vs${arrow("pass_ord")} ${U.passHeadBtn("lv", passSel)}</th><th>kb</th><th>Locus</th>${head}<th title="experiments of this line whose SE calls cover the locus">Called in</th><th title="rank of the locus by this line's own signal among the SEs it calls (1 = strongest); blank = not an SE here">SE rank</th><th title="copy number at the SE in this line (ratio to the line median); amplified when > 1.3">CN</th></tr></thead><tbody>${
        rows.map((r, i) => `<tr data-i="${i}" title="show in the browser"><td>${r.rank}</td><td><b>${U.esc(r.gene)}</b>${flagChip(r.flag)}</td><td>${U.passBadge(r.pass)}</td><td>${r.dist_kb}</td>
          <td class="mono lv-locus" title="${((r.end - r.start) / 1000).toFixed(1)} kb">${r.chrom}:${(+r.start).toLocaleString()}–${(+r.end).toLocaleString()}</td>${CMP.map(([k]) => fdrCell(L, r, k)).join("")}
          <td class="${r.called ? "" : "muted-s"}">${r.called ? r.called + "/" + L.n_experiments : "not an SE here"}</td><td>${r.signal_rank ?? ""}</td>
          <td class="${U.cnClass(r.cn)}">${r.cn != null ? r.cn.toFixed(2) : ""}</td></tr>`).join("")}</tbody></table>`;
    U.el("lv-cmps").innerHTML = `<div class="card"><div class="card-h"><h3>Specific super-enhancers, one row each, with the FDR of every comparison</h3>
      <span class="muted-s">${count}</span></div><div class="card-b"><div class="scroll lv-scroll">${body}</div></div>
      <div class="card-note">Ranked by the specificity score, which is the same for all four comparisons; they differ in
      the null (all lines, or only the other lines of the same lineage, disease or subtype), so each column says whether
      the super-enhancer stands out in that comparison. <b>&gt; 0.10</b> = tested, not called; <b>–</b> = not tested
      (hover the column header for why). <b>Called vs</b> repeats the called comparisons (${U.passBadge("ALDS")}); click
      a header to sort, and ▾ on Called vs to filter by combination (a lineage programme reads A·L without D).</div></div>`;
    U.wirePassHead("lv", L.rows, passSel, s => { passSel = s; renderTables(); });
    U.el("lv-cmps").querySelectorAll("th.th-btn").forEach(th => th.onclick = e => {
      e.stopPropagation();
      const k = th.dataset.k;
      sort = sort.k === k ? { k, asc: !sort.asc } : { k, asc: k !== "pass_ord" };
      renderTables();
    });
    const t = U.el("lv-t");
    if (t) t.addEventListener("click", e => {
      const tr = e.target.closest("tr[data-i]");
      if (tr && browser) { browser.search(locus(rows[+tr.dataset.i])); U.el("lv-igv").scrollIntoView({ behavior: "smooth", block: "start" }); }
    });
  }

  function renderExperiments() {
    const on = pickCoverage(line.experiments);
    U.el("lv-exps").innerHTML = line.experiments.map(e =>
      `<label class="lv-exp" title="${U.esc(e.srx)} · ${U.esc(e.study)} · ${U.esc(e.layout)}"><input type="checkbox" data-srx="${e.srx}"${on.has(e.srx) ? " checked" : ""}>
       ${e.srx} <span class="muted-s">${U.esc(e.study)}${e.layout ? " · " + (e.layout === "PAIRED" ? "PE" : "SE") : ""}</span></label>`).join("");
    U.el("lv-exps").onchange = async ev => {
      const cb = ev.target.closest("input[data-srx]");
      if (!cb || !browser) return;
      const e = line.experiments.find(x => x.srx === cb.dataset.srx);
      if (cb.checked) await addTrack(coverageTrack(e));
      else browser.removeTrackByName(coverageName(e));
    };
  }

  const coverageName = e => `H3K27ac ${e.srx} (${e.study || "?"}), raw RPM`;
  const coverageTrack = e => ({ name: coverageName(e), type: "wig", format: "bigwig", url: e.bw, color: COLORS.cov,
                                autoscale: true, height: 60 });

  async function addTrack(cfg) {                // one bad track must not take the browser down with it
    try { await browser.loadTrack(cfg); return true; }
    catch (err) { console.warn("track failed", cfg.name, err); return false; }
  }

  async function renderBrowser(my) {
    const host = U.el("lv-igv"), L = line;
    if (browser) { try { igv.removeBrowser(browser); } catch (_) {} browser = null; }
    host.innerHTML = "";
    const want = pending && pending.key === L.key ? pending.locus : null;
    pending = null;
    const first = L.rows[0];
    const b = await igv.createBrowser(host, { genome: "hg38", locus: want || (first ? locus(first) : "MYC"), tracks: [] });
    if (my !== token) { try { igv.removeBrowser(b); } catch (_) {} return; }   // a newer show() owns the host
    browser = b;
    const f = ext => `${DIR}${L.key}.${ext}`;
    for (const [k, tag, col] of CMP) {
      const c = L.comparisons[k];
      if (!c.testable || !c.n) continue;
      await addTrack({ name: `${cmpTitle(c, k)}: ${c.n.toLocaleString()}`, type: "annotation", format: "bed", url: f(`${tag}.bed.gz`),
                       indexed: false, color: col, displayMode: k === "all" ? "COLLAPSED" : "EXPANDED", height: k === "all" ? 25 : 35 });
      if (my !== token) return;
    }
    await addTrack({ name: `SEs called in ${L.name} (${L.n_called})`, type: "annotation", format: "bed", url: f("called.bed.gz"),
                     indexed: false, color: COLORS.called, useScore: true, displayMode: "COLLAPSED" });
    if (L.has_cn) await addTrack({ name: `Copy number, log2 (${CN_LABEL[L.cn_source] || L.cn_source})`, type: "wig", format: "bedgraph",
                     url: f("cn.bedgraph.gz"), indexed: false, color: COLORS.cn, altColor: COLORS.cnNeg, min: -2, max: L.cn_max || 3, height: 50 });
    const on = pickCoverage(L.experiments);
    for (const e of L.experiments) if (on.has(e.srx)) { if (my !== token) return; await addTrack(coverageTrack(e)); }
  }

  return { init, goto };
})();
