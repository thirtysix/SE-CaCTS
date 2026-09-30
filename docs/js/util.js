/* util.js, shared helpers for the SE-CaCTS dashboard. */
const U = (() => {
  const el = id => document.getElementById(id);
  const esc = s => String(s == null ? "" : s).replace(/[&<>"]/g, c =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;" }[c]));

  // levels: which are CALLS (panel-supported) vs RANKINGS-ONLY (gotcha 72)
  const LEVELS = [
    { key: "lineage", label: "Lineage", kind: "calls", file: "data/calls_lineage.tsv" },
    { key: "disease", label: "Primary disease", kind: "calls", file: "data/calls_disease.tsv" },
    { key: "subtype", label: "Subtype", kind: "rankings", file: "data/rank_subtype.tsv" },
    { key: "line", label: "Cell line", kind: "lines", file: null },   // per-line comparisons, data/lines/<key>.json
  ];

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

  // which of a line's four comparisons call this SE: "ADL" -> A, D, L lit, S dim
  const PASS = [["A", "all lines"], ["S", "same subtype"], ["D", "same primary disease"], ["L", "same lineage"]];
  const passBadge = p => {
    p = p || "";
    const on = PASS.filter(([k]) => p.includes(k)).map(([, w]) => w);
    const tip = on.length ? `specific vs ${on.join(", ")}` : "not called in any of the four comparisons (ranking only)";
    return `<span class="pass-b" title="${esc(tip)}">${PASS.map(([k]) => `<i class="${p.includes(k) ? "on p" + k : ""}">${k}</i>`).join("")}</span>`;
  };

  // analysis variants (manifest.variants): "main" is the default run; others stage calls_<level>.<key>.tsv
  let variant = "main";
  try { variant = localStorage.getItem("secacts-variant") || "main"; } catch (_) { /* default */ }
  const getVariant = () => variant;
  const setVariant = v => { variant = v; try { localStorage.setItem("secacts-variant", v); } catch (_) { /* not kept */ } };
  const variantFile = (file, v) => (!v || v === "main") ? file : file.replace(/\.(tsv|json)$/, `.${v}.$1`);

  return { el, esc, LEVELS, cnClass, fmtFdr, fmtJsd, ucsc, downloadTSV, passBadge, getVariant, setVariant, variantFile };
})();
