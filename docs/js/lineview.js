/* lineview.js, one cell line in igv.js: the SEs it calls, the SEs specific to it (vs the whole panel and vs
   its relatives), its copy number as scored, and its H3K27ac coverage streamed from ChIP-Atlas.
   Data: data/lines/ (phase2/scripts/61_stage_lines.py). Coverage is ChIP-Atlas's own per-experiment bigWig,
   read by range request (their server allows it for this site); it is raw RPM, not the S3-normalised signal
   the atlas compares, so two studies' tracks are not directly comparable in height. */
const LineView = (() => {
  const DIR = "data/lines/";
  const DEFAULT_COVERAGE = 3;                  // coverage tracks shown on open (one per study where possible)
  let index = null, combo = null, browser = null, line = null, token = 0;

  const COLORS = { vsrel: "#c0392b", vsall: "#e08214", called: "#7a7a7a", cn: "#2b6cb0", cnNeg: "#9b2c2c", cov: "#3d3d3d" };
  const CN_LABEL = { depmap_wgs: "DepMap WGS", depmap_mc_wes: "DepMap WES", cmp_wes: "CMP WES", ccle_snp6: "CCLE SNP6",
                     input_inferred: "inferred from ChIP input" };

  async function init() {
    index = await DataLoader.loadJSON(DIR + "index.json");
    const opts = Object.entries(index).map(([key, v]) => ({ key, label: v.name, search: `${v.search} | ${v.lineage}` }))
      .sort((a, b) => a.label.localeCompare(b.label));
    combo = Combo.make(U.el("lv-line"), key => show(key));
    combo.setOptions(opts);
    U.el("lv-count").textContent = `${opts.length} lines in this prototype`;
    const first = opts.find(o => o.label === "NIH:OVCAR-3") || opts[0];
    if (first) { combo.setValue(first.key); show(first.key); }
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
    line = await DataLoader.loadJSON(`${DIR}${key}.json`);
    if (my !== token) return;
    renderSummary(); renderTables(); renderExperiments();
    await renderBrowser(my);
  }

  function renderSummary() {
    const L = line, st = studiesOf(L.experiments).length;
    const cmp = L.comparison ? `compared with ${L.comparison.n_lines - 1} other line${L.comparison.n_lines === 2 ? "" : "s"} in ${U.esc(L.comparison.stratum)}` : "no relatives comparison";
    U.el("lv-summary").innerHTML = `
      <div class="lv-title"><h2>${U.esc(L.name)}</h2>
        <span class="muted-s">${U.esc(L.lineage)} › ${U.esc(L.disease)} › ${U.esc(L.subtype)}</span></div>
      <div class="lv-kpis">
        <div><b>${L.experiments.length}</b><span>experiments</span></div>
        <div><b>${st}</b><span>studies</span></div>
        <div><b>${L.n_called.toLocaleString()}</b><span>SEs called</span></div>
        <div><b>${L.vsall.n.toLocaleString()}</b><span>specific vs all</span></div>
        <div><b>${L.vsrel.n.toLocaleString()}</b><span>specific vs relatives</span></div>
      </div>
      <p class="muted-s">Copy number: ${U.esc(CN_LABEL[L.cn_source] || L.cn_source)}${L.has_cn ? "" : " (no track)"} · "vs relatives": ${cmp}.
        ${st < 2 ? " <b>One study only:</b> specific-SE calls need two independent studies, so none are made." : ""}</p>`;
  }

  function table(rows, id) {
    if (!rows.length) return `<p class="muted-s">None at FDR ≤ 0.1.</p>`;
    return `<table class="tbl lv-tbl" id="${id}"><thead><tr><th>#</th><th>Nearest gene</th><th>kb</th><th>Locus</th><th title="copy number at the SE in this line (ratio to the line median); the same flag as the SE atlas: amplified when > 1.3">CN</th><th>FDR</th></tr></thead><tbody>${
      rows.map((r, i) => `<tr data-i="${i}" title="show in the browser"><td>${r.rank}</td><td><b>${U.esc(r.gene)}</b></td><td>${r.dist_kb}</td>
        <td class="mono">${r.chrom}:${(r.start / 1e6).toFixed(2)} Mb</td>
        <td class="${U.cnClass(r.cn)}">${r.cn != null ? r.cn.toFixed(2) : ""}</td><td>${U.fmtFdr ? U.fmtFdr(r.fdr) : r.fdr}</td></tr>`).join("")}</tbody></table>`;
  }

  function renderTables() {
    const L = line;
    U.el("lv-vsrel").innerHTML = table(L.vsrel.top, "lv-t-rel");
    U.el("lv-vsall").innerHTML = table(L.vsall.top, "lv-t-all");
    U.el("lv-vsrel-n").textContent = L.vsrel.n > L.vsrel.top.length ? `top ${L.vsrel.top.length} of ${L.vsrel.n}` : `${L.vsrel.n}`;
    U.el("lv-vsall-n").textContent = L.vsall.n > L.vsall.top.length ? `top ${L.vsall.top.length} of ${L.vsall.n}` : `${L.vsall.n}`;
    for (const [id, rows] of [["lv-t-rel", L.vsrel.top], ["lv-t-all", L.vsall.top]]) {
      const t = U.el(id);
      if (t) t.addEventListener("click", e => {
        const tr = e.target.closest("tr[data-i]");
        if (tr && browser) { browser.search(locus(rows[+tr.dataset.i])); U.el("lv-igv").scrollIntoView({ behavior: "smooth", block: "start" }); }
      });
    }
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
    const start = L.vsrel.top[0] || L.vsall.top[0];
    browser = await igv.createBrowser(host, { genome: "hg38", locus: start ? locus(start) : "MYC", tracks: [] });
    if (my !== token) return;
    const f = ext => `${DIR}${L.key}.${ext}`;
    await addTrack({ name: `Specific vs relatives (${L.vsrel.n})`, type: "annotation", format: "bed", url: f("vsrel.bed.gz"),
                     indexed: false, color: COLORS.vsrel, displayMode: "EXPANDED", height: 40 });
    await addTrack({ name: `Specific vs all lines (${L.vsall.n})`, type: "annotation", format: "bed", url: f("vsall.bed.gz"),
                     indexed: false, color: COLORS.vsall, displayMode: "COLLAPSED" });
    await addTrack({ name: `SEs called in ${L.name} (${L.n_called})`, type: "annotation", format: "bed", url: f("called.bed.gz"),
                     indexed: false, color: COLORS.called, useScore: true, displayMode: "COLLAPSED" });
    if (L.has_cn) await addTrack({ name: `Copy number, log2 (${CN_LABEL[L.cn_source] || L.cn_source})`, type: "wig", format: "bedgraph",
                     url: f("cn.bedgraph.gz"), indexed: false, color: COLORS.cn, altColor: COLORS.cnNeg, min: -2, max: L.cn_max || 3, height: 50 });
    const on = pickCoverage(L.experiments);
    for (const e of L.experiments) if (on.has(e.srx)) { if (my !== token) return; await addTrack(coverageTrack(e)); }
  }

  return { init };
})();
