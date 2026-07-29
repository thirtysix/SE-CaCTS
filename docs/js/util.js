/* util.js — shared helpers for the SE-CaCTS dashboard. */
const U = (() => {
  const el = id => document.getElementById(id);
  const esc = s => String(s == null ? "" : s).replace(/[&<>"]/g, c =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;" }[c]));

  // levels: which are CALLS (panel-supported) vs RANKINGS-ONLY (gotcha 72)
  const LEVELS = [
    { key: "lineage", label: "Lineage", kind: "calls", file: "data/calls_lineage.tsv" },
    { key: "disease", label: "Primary disease", kind: "calls", file: "data/calls_disease.tsv" },
    { key: "subtype", label: "Subtype", kind: "rankings", file: "data/rank_subtype.tsv" },
    { key: "line", label: "Cell line", kind: "rankings", file: "data/rank_line.tsv" },
  ];

  // CN class from a group-mean copy-number ratio
  const cnClass = v => v == null || v === "" ? "" : (+v > 1.3 ? "cn-amp" : "cn-neu");
  const fmtFdr = v => (v == null || v === "" || isNaN(+v)) ? "—" : (+v).toFixed(3);
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

  return { el, esc, LEVELS, cnClass, fmtFdr, fmtJsd, ucsc, downloadTSV };
})();
