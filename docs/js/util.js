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
  return { el, esc, LEVELS, cnClass, fmtFdr, fmtJsd, ucsc, downloadTSV, passBadge, passOrd, passFilter, passHeadBtn, wirePassHead, closePassHead,
           getVariant, setVariant, variantFile };
})();
