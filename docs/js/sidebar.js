/* sidebar.js, the release card: which atlas release is on screen, and what it contains.
   Rendered at boot (not by a tab's lazy init), so the sidebar is filled whichever tab a link opens. */
const Sidebar = (() => {
  async function init() {
    const m = await DataLoader.loadJSON("data/meta.json");
    const r = m.release || {};
    const fmt = n => (n == null ? "n/a" : Number(n).toLocaleString());
    const delta = r.delta_lines ? ` <span class="rel-delta" title="change from the previous release">${r.delta_lines > 0 ? "+" : "−"}${fmt(Math.abs(r.delta_lines))}</span>` : "";
    const inferred = (m.cn_sources || {})["inferred from ChIP input"] || 0;
    const rows = [
      [fmt(m.n_lines) + delta, "cell lines", `cancer cell lines with QC-passed H3K27ac and copy number (measured${inferred ? `, or inferred from ChIP input for ${fmt(inferred)}` : ""}); replicate experiments are collapsed to the line`],
      [fmt(m.n_samples), "H3K27ac samples", `untreated or control experiments passing the ≥2,000-peak QC gate, of ${fmt(m.n_pull)} pulled (ChIP-Atlas, plus SRA data it lacks, processed with its pipeline)`],
      [fmt(m.n_ses), "super-enhancer loci", "the union catalogue of super-enhancers across all scored samples"],
      [fmt(m.n_lineages), "lineages", "Oncotree lineages: specificity calls are supported at this level"],
      [fmt(m.n_diseases), "primary diseases", "Oncotree primary diseases: specificity calls are supported at this level"],
      [fmt(m.n_lineage_calls), "lineage calls", "super-enhancer × lineage tests passing the permutation FDR ≤ 0.10"],
      [fmt(m.n_disease_calls), "disease calls", "super-enhancer × primary-disease tests passing the permutation FDR ≤ 0.10"],
      ...(m.n_subtype_calls != null ? [[fmt(m.n_subtype_calls), "subtype calls", "super-enhancer × subtype tests passing the permutation FDR ≤ 0.10, in subtypes of two or more lines"]] : []),
    ];
    const cn = Object.entries(m.cn_sources || {}).filter(([, v]) => v > 0)
      .map(([k, v]) => `${v} ${U.esc(k)}`).join(" · ");
    U.el("snap").innerHTML = `
      <a class="rel" href="#about" data-goto="releases" title="release notes and the history of every release">
        <div class="rel-h"><span class="rel-v">release ${U.esc(r.version || "")}</span>
          <span class="rel-d">${U.esc(r.date || "")}</span></div>
        ${r.title ? `<div class="rel-t">${U.esc(r.title)}</div>` : ""}
      </a>
      <div class="snap-grid">${rows.map(([b, l, t]) => `<b title="${U.esc(t)}">${b}</b><span title="${U.esc(t)}">${l}</span>`).join("")}</div>
      ${cn ? `<div class="rel-cn" title="source of the copy number each line is corrected with">CN: ${cn}</div>` : ""}`;
    // the card links to the Releases section of the About tab
    const card = U.el("snap").querySelector(".rel");
    if (card) card.addEventListener("click", e => {
      e.preventDefault();
      TabManager.switchTab("about");
      const go = () => { const t = document.getElementById("releases"); if (t) t.scrollIntoView({ behavior: "smooth" }); else setTimeout(go, 80); };
      go();
    });
  }
  return { init };
})();
