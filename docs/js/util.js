/* util.js, shared helpers for the SE-CaCTS dashboard. */
const U = (() => {
  const el = id => document.getElementById(id);
  const esc = s => String(s == null ? "" : s).replace(/[&<>"]/g, c =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;" }[c]));

  // levels: which are CALLS (panel-supported) vs RANKINGS-ONLY (gotcha 72)
  const LEVELS = [
    { key: "lineage", label: "Lineage", kind: "calls", file: "data/calls_lineage.tsv" },
    { key: "disease", label: "Primary disease", kind: "calls", file: "data/calls_disease.tsv" },
    { key: "subtype", label: "Subtype", kind: "calls", file: "data/calls_subtype.tsv" },   // calls from v3.1
    { key: "line", label: "Cell line", kind: "lines", file: null },   // per-line comparisons, data/lines/<key>.json
  ];

  // v3.1 copy-number labels (fused calling, FINDINGS §55). cns: how the two statistics agree (r CN-robust: called with
  // and without correction; u CN-unmasked: only with it; a / g: only WITHOUT it, amplicon-driven at CN >= 2 or
  // gain-dependent below). sel: how the group's own experiments call the SE (g / A / H: only through gain, amplification
  // or high-level amplification). Never a merged FDR: the "only without correction" calls keep their own list.
  const SRC = { depmap_wgs: "DepMap WGS", cmp_wes: "Cell Model Passports WES", depmap_mc_wes: "DepMap WES",
                ccle_snp6: "CCLE SNP6", input_inferred: "inferred from ChIP input" };
  const srcText = s => String(s || "").split(",").filter(Boolean).map(x => SRC[x] || x).join(", ");
  const CNS = {
    u: ["CN-unmasked", "cns-u", "specific only WITH copy-number correction: without it, amplification in other lines hides this specificity"],
    a: ["amplicon-driven", "cns-a", "passes the FDR only WITHOUT copy-number correction, at a locus these lines carry at copy number 2 or more: the specificity comes from extra copies"],
    g: ["gain-dependent", "cns-g", "passes the FDR only WITHOUT copy-number correction, at a low-level gain (copy number below 2) in these lines"],
  };
  const SELT = { g: ["gain", "low-level gain (copy number below 2)"], A: ["amplified", "amplification (copy number 2-3)"],
                 H: ["high-level", "high-level amplification (copy number 3 or more)"] };
  const cnChips = (r, who) => {
    const src = r.cn_src ? `. Copy number from ${srcText(r.cn_src)}` : "", cn = r.cn_mean != null && r.cn_mean !== "" ? ` (copy number ${(+r.cn_mean).toFixed(2)})` : "";
    let h = "";
    const c = CNS[r.cns];
    if (c) h += `<span class="cns-chip ${c[1]}" title="${esc(c[2] + cn + src)}">${c[0]}</span>`;
    const s = SELT[r.sel], k = r.n_exp ? ` ${r.n_amp}/${r.n_exp}` : "";
    if (s) h += `<span class="cns-chip cns-sel" title="${esc(`${who} call this a super-enhancer only through ${s[1]}: on copy-number-corrected signal it falls below the cutoff${r.n_exp ? ` (${r.n_amp} of the ${r.n_exp} experiments that call it)` : ""}${src}`)}">SE via ${s[0]}${k}</span>`;
    return h;
  };
  const CN_FILTERS = [["all", "All calls"], ["r", "CN-robust"], ["u", "CN-unmasked"], ["sel", "SE via gain / amplification"],
                      ["dep", "Only without correction"]];
  const cnPass = (r, f) => f === "all" || f === "dep" || (f === "sel" ? ["g", "A", "H"].includes(r.sel) : r.cns === f);

  // CN class from a group-mean copy-number ratio
  const cnClass = v => v == null || v === "" ? "" : (+v > 1.3 ? "cn-amp" : "cn-neu");
  const fmtFdr = v => (v == null || v === "" || isNaN(+v)) ? "n/a" : (+v).toFixed(3);
  const fmtJsd = v => (+v).toFixed(3);

  // a UCSC out-link for an SE locus
  const ucsc = (chrom, start, end) => chrom
    ? `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${chrom}:${Math.max(1, start - 2000)}-${end + 2000}`
    : "";

  const downloadTSV = (filename, cols, rows) => {
    const cell = v => v == null ? "" : String(v).replace(/[\t\r\n]/g, " ");
    const lines = [cols.map(c => c.label).join("\t")].concat(
      rows.map(r => cols.map(c => cell(c.get ? c.get(r) : r[c.key])).join("\t")));
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([lines.join("\n")], { type: "text/tab-separated-values" }));
    a.download = filename; document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(a.href);
  };

  // which of a line's four comparisons call this SE, broad to narrow: "ALD" -> A, L, D lit, S dim
  const PASS = [["A", "all lines"], ["L", "same lineage"], ["D", "same primary disease"], ["S", "same subtype"]];
  const passBadge = p => {
    p = p || "";
    const on = PASS.filter(([k]) => p.includes(k)).map(([, w]) => w);
    const tip = on.length ? `specific vs ${on.join(", ")}` : "not called in any of the four comparisons (ranking only)";
    return `<span class="pass-b" title="${esc(tip)}">${PASS.map(([k]) => `<i class="${p.includes(k) ? "on p" + k : ""}">${k}</i>`).join("")}</span>`;
  };

  // sort key for a pass pattern: A=8, L=4, D=2, S=1, so descending lists ALDS, ALD, ALS, AL, ... A, then none
  const passOrd = p => (p || "").split("").reduce((t, k) => t + ({ A: 8, L: 4, D: 2, S: 1 }[k] || 0), 0);
  // checkbox per pass pattern present in `rows` (with counts); `sel` = Set of shown patterns or null for all
  const passFilter = (el, rows, sel, onChange) => {
    const cnt = {};
    rows.forEach(r => { const p = r.pass || ""; cnt[p] = (cnt[p] || 0) + 1; });
    const pats = Object.keys(cnt).sort((a, b) => passOrd(b) - passOrd(a));
    const on = p => !sel || sel.has(p);
    el.innerHTML = pats.map(p => `<label class="pass-opt" title="${p ? "called vs " + p.split("").join(", ") : "not called in any comparison (rankings)"}"><input type="checkbox" data-p="${p}"${on(p) ? " checked" : ""}>${passBadge(p)}<span class="muted-s">${cnt[p]}</span></label>`).join("") +
      `<button class="dl-btn pass-all" title="show every combination">all</button>`;
    el.querySelectorAll("input[data-p]").forEach(cb => cb.onchange = () => {
      const s = new Set([...el.querySelectorAll("input[data-p]:checked")].map(x => x.dataset.p));
      onChange(s.size === pats.length ? null : s);
    });
    el.querySelector(".pass-all").onclick = () => onChange(null);
  };

  // the filter lives in the "Called vs" column header: a ▾ button opens it in a popover that stays open while
  // boxes are ticked (the table re-renders, then calls this again), and closes on an outside click or Esc
  let popFor = null;
  const passHeadBtn = (tableId, sel) => `<button class="pass-dd${sel ? " on" : ""}" data-t="${tableId}" title="filter by combination of called comparisons${sel ? " (filtered)" : ""}" aria-label="filter">▾</button>`;
  const wirePassHead = (tableId, rows, sel, onChange) => {
    let pop = document.getElementById("pass-pop");
    if (!pop) {
      pop = document.createElement("div"); pop.id = "pass-pop"; pop.className = "pass-pop"; document.body.appendChild(pop);
      document.addEventListener("click", e => { if (popFor && !pop.contains(e.target) && !e.target.closest(".pass-dd")) { popFor = null; pop.style.display = "none"; } });
      document.addEventListener("keydown", e => { if (e.key === "Escape" && popFor) { popFor = null; pop.style.display = "none"; } });
      window.addEventListener("hashchange", () => { popFor = null; pop.style.display = "none"; });
    }
    const btn = document.querySelector(`.pass-dd[data-t="${tableId}"]`);
    const place = () => {
      if (popFor !== tableId || !btn) { if (popFor === tableId || !popFor) pop.style.display = "none"; return; }
      passFilter(pop, rows, sel, onChange);
      const r = btn.getBoundingClientRect();
      pop.style.display = "flex";
      pop.style.left = `${Math.max(8, Math.min(r.left + window.scrollX, window.scrollX + document.documentElement.clientWidth - pop.offsetWidth - 8))}px`;
      pop.style.top = `${r.bottom + window.scrollY + 4}px`;
    };
    if (btn) btn.onclick = e => { e.stopPropagation(); popFor = popFor === tableId ? null : tableId; place(); };
    place();
  };

  // analysis variants (manifest.variants): "main" is the default run; others stage calls_<level>.<key>.tsv
  let variant = "main";
  try { variant = localStorage.getItem("secacts-variant") || "main"; } catch (_) { /* default */ }
  const getVariant = () => variant;
  const setVariant = v => { variant = v; try { localStorage.setItem("secacts-variant", v); } catch (_) { /* not kept */ } };
  const variantFile = (file, v) => (!v || v === "main") ? file : file.replace(/\.(tsv|json)$/, `.${v}.$1`);

  const closePassHead = () => { popFor = null; const p = document.getElementById("pass-pop"); if (p) p.style.display = "none"; };
  const LINKS = {
    chipatlas: "https://chip-atlas.org/", geo: "https://www.ncbi.nlm.nih.gov/geo/", sra: "https://www.ncbi.nlm.nih.gov/sra",
    depmap: "https://depmap.org/portal/", cmp: "https://cellmodelpassports.sanger.ac.uk/", ccle: "https://sites.broadinstitute.org/ccle/",
    cellosaurus: "https://www.cellosaurus.org/", ncit: "https://ncit.nci.nih.gov/", oncotree: "https://oncotree.mskcc.org/",
    ensembl: "https://apr2022.archive.ensembl.org/", igv: "https://github.com/igvteam/igv.js",
    rose: "http://younglab.wi.mit.edu/super_enhancer_code.html", rose2: "https://github.com/linlabbcm/rose2",
    s3norm: "https://github.com/guanjue/S3norm", cacts: "https://doi.org/10.1126/sciadv.abf6123",
    cactscode: "https://github.com/lawrenson-lab/CaCTS", pycacts: "https://github.com/thirtysix/pyCaCTS",
    cnrose: "https://github.com/thirtysix/SE-CaCTS/tree/main/cnrose", repo: "https://github.com/thirtysix/SE-CaCTS",
    depmap_ref: "https://doi.org/10.1016/j.cell.2017.06.010", cmp_ref: "https://doi.org/10.1093/nar/gky872",
    ccle_ref: "https://doi.org/10.1038/s41586-019-1186-3", cbioportal: "https://www.cbioportal.org/study/summary?id=ccle_broad_2019",
    encode: "https://www.encodeproject.org/", encode_ref: "https://pubmed.ncbi.nlm.nih.gov/32728249/",
  };
  const link = (k, text) => `<a href="${LINKS[k]}" target="_blank" rel="noopener">${text}</a>`;

  // ---- genes near an SE (data/se_genes.json, 75_stage_se_genes.py): every protein-coding gene within 100 kb of the
  // SE, nearest first, as [gene, kb] (kb 0 = overlaps); [gene, kb, 1] = the nearest gene when none lies within 100 kb.
  // data/expr_specific.json: per level and group, the genes that are group-specific in DepMap expression (⇌).
  // data/se_lnc.json: HGNC-named lncRNA genes within 100 kb, [name, kb], listed after the protein-coding genes.
  let genesP = null, SEG = {}, EXPR = {}, SEL = {};
  const loadGenes = () => genesP || (genesP = Promise.all([DataLoader.loadJSON("data/se_genes.json"),
      DataLoader.loadJSON("data/expr_specific.json").catch(() => ({})), DataLoader.loadJSON("data/se_lnc.json").catch(() => ({}))])
    .then(([g, e, l]) => { SEG = g; SEL = l; for (const lv in e) { EXPR[lv] = {}; for (const k in e[lv]) EXPR[lv][k] = new Set(e[lv][k]); } })
    .catch(() => { genesP = null; }));
  const kbTxt = x => x[1] === 0 ? "overlaps" : `${x[1]} kb`;
  // a COPY of the row with the gene list (gene / dist_kb = the nearest); cached rows are shared, never mutate them
  const withGenes = r => {
    const gl = SEG[r.se] || (r.gene ? [[r.gene, r.dist_kb, 1]] : []);
    const ll = SEL[r.se] || [];
    return { ...r, gene0: r.gene, genes: gl, lnc: ll, gene: gl.length ? gl[0][0] : "", dist_kb: gl.length ? gl[0][1] : "",
             genes_100kb: gl.filter(x => !x[2]).map(x => `${x[0]}:${x[1]}`).join(";"),
             lncrna_100kb: ll.map(x => `${x[0]}:${x[1]}`).join(";") };
  };
  const geneMatch = (r, q) => !q || (r.se || "").toLowerCase() === q || (r.genes || []).some(x => x[0].toLowerCase().includes(q))
                                 || (r.lnc || []).some(x => x[0].toLowerCase().includes(q));
  const exprSpec = (lv, grp, g) => !!(lv && EXPR[lv] && EXPR[lv][grp] && EXPR[lv][grp].has(g));
  // the gene cell: the nearest gene(s), then any gene matching the filter, then "+N" (hover lists all; click expands).
  // lv/grp name the group whose expression marks the genes (⇌); r.rho is the staged nearest-gene correlation.
  function geneCell(r, { q = "", lv = null, grp = r.group, link = null } = {}) {
    const gl = r.genes || [];
    if (!gl.length) return `<span class="muted-s">n/a</span>`;
    const outside = !!gl[0][2];
    const mark = g => !exprSpec(lv, grp, g) ? "" : ` <span class="conc-badge" title="${esc(g)} is itself specific to ${esc(grp)} in DepMap expression (cross-layer concordant)${g === r.gene0 && r.rho != null && r.rho !== "" ? `; its expression tracks this SE's H3K27ac across lines (Spearman ρ = ${r.rho})` : ""}. See the Concordance tab.">⇌${g === r.gene0 && r.rho != null && r.rho !== "" ? " " + (+r.rho).toFixed(2) : ""}</span>`;
    const one = (x, first) => `<span class="g1">${first && link ? link(x[0]) : esc(x[0])}${mark(x[0])}<span class="th-sub"> ${kbTxt(x)}</span></span>`;
    // inline: the nearest gene and any tied with it (several genes often overlap one SE), then filter matches
    const ll = r.lnc || [];
    const lncOne = x => `<span class="g1 lnc" title="lncRNA">${esc(x[0])}<span class="th-sub"> ${kbTxt(x)}</span></span>`;
    const lead = gl.filter((x, i) => i === 0 || (!x[2] && x[1] === gl[0][1]) || (q && x[0].toLowerCase().includes(q)));
    const leadL = q ? ll.filter(x => x[0].toLowerCase().includes(q)) : [];
    const rest = gl.filter(x => !lead.includes(x)), restL = ll.filter(x => !leadL.includes(x));
    const all = gl.map(x => `${x[0]} ${kbTxt(x)}${exprSpec(lv, grp, x[0]) ? " ⇌" : ""}`).join(", ");
    const allL = ll.length ? ` lncRNAs within 100 kb: ${ll.map(x => `${x[0]} ${kbTxt(x)}`).join(", ")}.` : "";
    const tip = (outside ? `no protein-coding gene within 100 kb; the nearest is ${all}.` : `protein-coding genes within 100 kb of this super-enhancer, nearest first: ${all}.`) + allL + " Proximity only, not a scored link.";
    // <wbr>: the items are nowrap and joined without whitespace, so a long tie (HOXA cluster) could not wrap
    const sep = `<span class="sep">·</span><wbr>`;
    const nMore = rest.length + restL.length;
    const more = [...rest.map(x => one(x, false)), ...restL.map(lncOne)];
    return `<span class="genes" title="${esc(tip)}">${lead.map((x, i) => one(x, i === 0)).join(sep)}${outside ? ` <span class="genes-out">(no protein-coding gene within 100 kb)</span>` : ""}${leadL.length ? sep + leadL.map(lncOne).join(sep) : ""}${nMore
      ? ` <button class="genes-more" type="button" title="show the other ${nMore} gene${nMore > 1 ? "s" : ""} within 100 kb${restL.length ? ` (${restL.length} lncRNA${restL.length > 1 ? "s" : ""})` : ""}">+${nMore}</button><span class="genes-rest">${sep}${more.join(sep)}</span>` : ""}</span>`;
  }
  // "+N" expands in place; capture phase so a click on it never also triggers the row (Genomic View rows jump IGV)
  document.addEventListener("click", e => {
    const b = e.target.closest && e.target.closest(".genes-more");
    if (b) { e.preventDefault(); e.stopPropagation(); b.closest(".genes").classList.add("open"); }
  }, true);

  return { el, esc, LEVELS, LINKS, link, cnClass, cnChips, CN_FILTERS, cnPass, srcText, fmtFdr, fmtJsd, ucsc, downloadTSV, passBadge, passOrd, passFilter, passHeadBtn, wirePassHead, closePassHead,
           getVariant, setVariant, variantFile, loadGenes, withGenes, geneMatch, geneCell, exprSpec };
})();
