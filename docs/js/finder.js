/* finder.js, type a gene symbol, see every lineage / disease / subtype / cell line where an SE within 100 kb of it, or in
   Hi-C contact with it, is specific. Groups: data/gene_index.json (75_stage_se_genes.py: every gene within 100 kb of each
   call, plus the SE's Hi-C contact gene, h = 1). Cell lines: data/gene_lines/<first character>.json, loaded on demand. */
const Finder = (() => {
  let index = null, lines = null, token = 0;
  const LV = [["lineage", "Lineage"], ["disease", "Primary disease"], ["subtype", "Subtype"]];
  const LINE_SHOW = 30;                                   // cell-line chips shown before "+N more"
  const shards = {};
  const shardOf = g => /^[A-Z]/.test(g) ? g[0] : "_";
  const shard = ch => shards[ch] || (shards[ch] = DataLoader.loadJSON(`data/gene_lines/${ch}.json`).catch(() => ({})));

  const kbOf = (d, hic) => d == null ? "" : d === 0 ? (hic ? "TSS inside" : "overlaps") : `${d} kb`;

  function groupChip(g, h) {
    const conc = h.c === 1 ? "; and the gene is itself group-specific in expression (cross-layer concordant)" : "";
    const hicOnly = h.h && h.n === 0 && h.o == null && h.t == null;
    const rel = h.t === "lnc" ? " (a lncRNA; nearest is counted among protein-coding genes only)" : h.o ? ", the nearest gene (none lies within 100 kb)" : h.n ? ", the nearest gene" : ", not the nearest gene";
    const hic = h.h ? "; it is the super-enhancer's Hi-C contact gene (the highest contact among the genes with a TSS within 1 Mb, ENCODE Hi-C of 14 cancer lines)" : "";
    const where = hicOnly ? `is in Hi-C contact with, ${h.d} kb from its TSS,` : h.d === 0 ? "overlaps" : `lies ${h.d} kb from`;
    const d = kbOf(h.d, hicOnly);
    return `<a class="chip chip-link${h.h ? " hic" : ""}" href="#atlas" data-lv="${h.lv}" data-g="${U.esc(h.g)}" data-gene="${U.esc(g)}" title="${U.esc(g)} ${where} a super-enhancer specific to ${U.esc(h.g)} (${h.lv} level) at rank ${h.r}, permutation FDR ${h.fdr}${hicOnly ? "" : rel}${hic}${conc}. Click to open it in the SE atlas.">
      <b>${U.esc(h.g)}</b> <span class="r">#${h.r}</span>${d ? ` <span class="r" style="opacity:.75">· ${h.h && (hicOnly || h.o) ? "Hi-C " : ""}${d}</span>` : ""}${h.c === 1 ? ` <span class="cc" title="cross-layer concordant: the gene is itself group-specific in DepMap expression">⇌</span>` : ""}</a>`;
  }

  function lineChip(g, e) {
    const [li, rank, ps, d, fl] = e, [key, name, lineage, grp] = lines[li];
    const hicOnly = fl === 1;                               // the Hi-C contact gene beyond 100 kb (d = kb to its TSS)
    const what = fl & 2 ? " (a lncRNA)" : fl & 4 ? ", the nearest gene (none lies within 100 kb)" : fl & 8 ? ", the nearest gene" : "";
    const hic = fl & 1 ? "; it is the super-enhancer's Hi-C contact gene (ENCODE Hi-C of 14 cancer lines)" : "";
    return `<a class="chip chip-link${fl & 1 ? " hic" : ""}" href="#atlas" data-lv="line" data-g="${U.esc(grp)}" data-gene="${U.esc(g)}" title="${U.esc(g)} ${hicOnly ? `is in Hi-C contact with, ${d} kb from its TSS,` : d === 0 ? "overlaps" : `lies ${d} kb from`} a super-enhancer called in ${U.esc(name)} (${U.esc(lineage)}) at rank ${rank}${what}${hic}. Click to open the line in the SE atlas.">
      <b>${U.esc(name)}</b> <span class="r">#${rank}</span> ${U.passBadge(ps)}<span class="r" style="opacity:.75"> · ${(fl & 1) && (hicOnly || (fl & 4)) ? "Hi-C " : ""}${kbOf(d, hicOnly)}</span></a>`;
  }

  async function render(qRaw) {
    const my = ++token, q = qRaw.trim().toUpperCase();
    const box = U.el("finder-res");
    if (!q) { box.innerHTML = `<div class="empty">Type a gene symbol above.</div>`; return; }
    const L = await shard(shardOf(q));
    if (!lines) lines = await DataLoader.loadJSON("data/gene_lines/lines.json").catch(() => []);
    if (my !== token) return;                              // a newer keystroke is already rendering
    // exact match first, else prefix matches, over the group index and the cell-line index
    const known = g => index[g] || L[g];
    let genes = known(q) ? [q] : [...new Set([...Object.keys(index), ...Object.keys(L)])].filter(g => g.startsWith(q)).sort().slice(0, 12);
    if (!genes.length) {
      box.innerHTML = `<div class="empty">No specific super-enhancer lies within 100 kb of <b>${U.esc(q)}</b>, or is in Hi-C contact with it,
        at any level.</div>`;
      return;
    }
    const v = U.getVariant();
    box.innerHTML = genes.map(g => {
      const hits = index[g] || [];
      const byGroup = {};                                  // best (lowest) rank per group, keeping the level
      hits.forEach(h => { const k = h.lv + "|" + h.g; if (!byGroup[k] || h.r < byGroup[k].r || (h.r === byGroup[k].r && h.d < byGroup[k].d)) byGroup[k] = h; });
      const best = Object.values(byGroup).sort((a, b) => a.r - b.r);
      const sec = LV.map(([lv, label]) => {
        const c = best.filter(h => h.lv === lv);
        return c.length ? `<div class="fr-lv"><span class="fr-l">${label}<small>${c.length}</small></span><div class="chips">${c.map(h => groupChip(g, h)).join("")}</div></div>` : "";
      }).join("");
      const byLine = {};                                   // best (lowest) rank per cell line
      (L[g] || []).forEach(e => { if (!byLine[e[0]] || e[1] < byLine[e[0]][1]) byLine[e[0]] = e; });
      const ls = Object.values(byLine).sort((a, b) => a[1] - b[1]);
      const lineSec = ls.length ? `<div class="fr-lv"><span class="fr-l">Cell line<small>${ls.length}</small>${v !== "main" ? `<small title="per-line calls come from the default analysis only">default analysis</small>` : ""}</span><div class="chips">${
        ls.slice(0, LINE_SHOW).map(e => lineChip(g, e)).join("")}${ls.length > LINE_SHOW
        ? `<button class="fr-more" type="button" title="show the other ${ls.length - LINE_SHOW} cell lines">+${ls.length - LINE_SHOW}</button><span class="fr-rest">${ls.slice(LINE_SHOW).map(e => lineChip(g, e)).join("")}</span>` : ""}</div></div>` : "";
      const isLnc = hits.length && hits.every(h => h.t === "lnc");
      return `<div class="fr card"><div class="fr-head"><span class="sym">${U.esc(g)}${isLnc ? ` <span class="lnc-tag" title="long non-coding RNA (HGNC-named, GENCODE v36)">lncRNA</span>` : ""}</span></div>${
        sec || `<div class="fr-lv"><span class="fr-l">Groups</span><span class="muted-s">no group-level call</span></div>`}${lineSec}</div>`;
    }).join("");
  }

  const cache = {};
  async function load() {
    const v = U.getVariant(), f = U.variantFile("data/gene_index.json", v);
    index = cache[v] || (cache[v] = await DataLoader.loadJSON(f));
    const man = await DataLoader.loadJSON("data/manifest.json"), vs = man.variants || [];
    const sel = U.el("finder-variant");
    sel.parentElement.style.display = vs.length > 1 ? "inline-flex" : "none";
    sel.innerHTML = vs.map(x => `<option value="${x.key}"${x.key === v ? " selected" : ""}>${U.esc(x.label)}</option>`).join("");
    const vi = vs.find(x => x.key === v);
    U.el("finder-variant-info").title = vi ? `${vi.label}: ${vi.desc || ""}` : "";
  }

  async function init() {
    await load();
    const inp = U.el("finder-input");
    inp.addEventListener("input", () => render(inp.value));
    // a result opens its group (or cell line) in the SE atlas, filtered to the gene; "+N" shows the other lines
    U.el("finder-res").addEventListener("click", e => {
      const more = e.target.closest(".fr-more");
      if (more) { more.nextElementSibling.classList.add("open"); more.remove(); return; }
      const a = e.target.closest("a.chip-link");
      if (!a) return;
      e.preventDefault();
      Atlas.open(a.dataset.lv, a.dataset.g, a.dataset.gene);
      location.hash = "atlas";
    });
    U.el("finder-variant").onchange = async e => { U.setVariant(e.target.value); await load(); render(inp.value); };
    render(inp.value || "");
  }
  return { init };
})();
