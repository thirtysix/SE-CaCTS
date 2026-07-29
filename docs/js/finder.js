/* finder.js — type a gene symbol, see every lineage / disease where an SE near it is group-specific. */
const Finder = (() => {
  let index = null;

  function render(qRaw) {
    const q = qRaw.trim().toUpperCase();
    const box = U.el("finder-res");
    if (!q) { box.innerHTML = `<div class="empty">Type a gene symbol above.</div>`; return; }
    // exact match first, else prefix matches
    let genes = index[q] ? [q] : Object.keys(index).filter(g => g.startsWith(q)).sort().slice(0, 12);
    if (!genes.length) {
      box.innerHTML = `<div class="empty">No lineage- or disease-specific super-enhancer is near <b>${U.esc(q)}</b>.
        <span class="muted-s">(The finder covers the two call levels; a gene absent here may still rank at subtype / cell-line level in the SE atlas.)</span></div>`;
      return;
    }
    box.innerHTML = genes.map(g => {
      const hits = index[g];
      // best (lowest) rank per group, keeping the level
      const byGroup = {};
      hits.forEach(h => { const k = h.lv + "|" + h.g; if (!byGroup[k] || h.r < byGroup[k].r) byGroup[k] = h; });
      const chips = Object.values(byGroup).sort((a, b) => a.r - b.r).map(h => {
        const conc = h.c === 1 ? "; and the gene is itself group-specific in expression (cross-layer concordant)" : "";
        return `<span class="chip" title="${U.esc(g)} is the nearest gene to a super-enhancer specific to ${U.esc(h.g)} (${h.lv} level) at rank ${h.r}, permutation FDR ${h.fdr}${conc}">
          <b>${U.esc(h.g)}</b> <span class="r">#${h.r}</span>${h.c === 1 ? ` <span class="cc" title="cross-layer concordant: the gene is itself group-specific in DepMap expression">⇌</span>` : ""}</span>`;
      }).join("");
      return `<div class="fr card"><span class="sym">${U.esc(g)}</span><div class="chips">${chips}</div></div>`;
    }).join("");
  }

  async function init() {
    index = await DataLoader.loadJSON("data/gene_index.json");
    const inp = U.el("finder-input");
    inp.addEventListener("input", () => render(inp.value));
    render(inp.value || "");
  }
  return { init };
})();
