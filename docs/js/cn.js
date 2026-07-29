/* cn.js — the call-based CN ablation: correction rescues real specificity and removes recurrent amplicons. */
const CN = (() => {
  async function init() {
    const abl = await DataLoader.loadJSON("data/cn_ablation.json");

    // the flow per level: uncorrected -> corrected, with rescued / amplicon-driven split
    U.el("cn-flow").innerHTML = abl.summary.map(s => `
      <div class="card">
        <div class="card-h"><h3>${U.esc(s.level)}</h3><span class="muted-s">permutation FDR ≤ 0.10</span></div>
        <div class="card-b">
          <div class="abl-flow" title="calls significant with CN correction but NOT without it — real specificity that amplicon variance in the null was masking">
            <span class="n">${s.rescued.toLocaleString()}</span>
            <span class="arw">◀ rescued by correction</span>
          </div>
          <div class="abl-grid">
            <div class="stat" title="super-enhancers passing the permutation FDR ≤ 0.10 with copy-number correction applied"><div class="k good">${s.corrected.toLocaleString()}</div><div class="l">specific SEs, CN-corrected</div></div>
            <div class="stat" title="calls that survive in both arms — significant with or without correction (copy-neutral specificity)"><div class="k">${(s.corrected - s.rescued).toLocaleString()}</div><div class="l">also called without correction</div></div>
            <div class="stat" title="calls significant WITHOUT correction but not WITH it — significant only because the locus is amplified. Median copy number ${s.amp_median_cn}× vs neutral 1×"><div class="k" style="color:var(--hot)">${s.amplicon_driven}</div><div class="l">amplicon-driven false calls removed<br><span class="muted-s">median CN ${s.amp_median_cn}×</span></div></div>
            <div class="stat" title="median copy-number ratio of the rescued calls — ≈1 means they sit at neutral copy number, so they are real specificity rather than amplicon artifacts"><div class="k">${s.resc_median_cn}</div><div class="l">median CN of rescued calls<br><span class="muted-s">(neutral = real specificity)</span></div></div>
          </div>
        </div>
      </div>`).join("");

    // the amplicon-driven calls — every one a named recurrent lineage amplicon
    U.el("cn-amp").innerHTML = abl.amplicon.map(a => `
      <div class="amp-row" title="${U.esc(a.gene)} is amplified at ${a.cn}× in ${U.esc(a.group)}; a super-enhancer there passed the FDR without CN correction and dropped after it — a recurrent amplicon, not lineage identity">
        <span class="g">${U.esc(a.gene)}</span>
        <span class="grp">${U.esc(a.group)} <span class="th-sub">(${a.level})</span></span>
        <span class="cnv">CN ${a.cn}×</span>
      </div>`).join("");
  }
  return { init };
})();
