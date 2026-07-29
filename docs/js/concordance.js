/* concordance.js — Phase-6 cross-layer check: specific SEs sit next to genes specific to the same group. */
const Concordance = (() => {
  async function init() {
    const c = await DataLoader.loadJSON("data/concordance.json");

    U.el("conc-summary").innerHTML = c.summary.map(s => `
      <div class="card">
        <div class="card-h"><h3>${U.esc(s.level)}</h3><span class="muted-s">${s.enrichment}× over background</span></div>
        <div class="card-b">
          <div class="abl-grid">
            <div class="stat" title="of all (specific SE, nearby gene) pairs, the fraction where the gene is itself specific to the SAME group in DepMap expression — the headline concordance"><div class="k good">${s.per_pair}%</div><div class="l">SE-proximal genes specific to the same group</div></div>
            <div class="stat" title="the rate at which any gene is specific to any group — the null expectation the ${s.per_pair}% is measured against (${s.enrichment}× enrichment)"><div class="k">${s.background}%</div><div class="l">background (all genes, all groups)</div></div>
            <div class="stat" title="the same test with each SE scored against a RANDOM other group — collapses to background, proving the effect is same-group-specific"><div class="k">${s.shuffled}%</div><div class="l">group-shuffled control</div></div>
            <div class="stat" title="the fraction of specific SEs that have at least one concordant group-specific gene within 100 kb"><div class="k">${s.per_se_any}%</div><div class="l">SEs with ≥1 concordant gene in 100 kb</div></div>
          </div>
        </div>
      </div>`).join("");

    // distance decay — the control that separates a local regulatory link from a lineage confound.
    // teal bar = concordance; a grey tick marks the group-shuffled chance level (stays flat ~2.5%).
    const maxC = Math.max(...c.distance.map(d => d.concordant));
    U.el("conc-dist").innerHTML = c.distance.map(d => `
      <div class="dbar" title="${d.n.toLocaleString()} SE–gene pairs; median Spearman ρ = ${d.rho}">
        <span class="dl">${d.bin}</span>
        <div class="track">
          <div class="fill" style="width:${(d.concordant / maxC * 100).toFixed(1)}%"></div>
          <div class="mark" style="left:${(d.shuffled / maxC * 100).toFixed(1)}%" title="group-shuffled control: ${d.shuffled}%"></div>
        </div>
        <span class="dv">${d.concordant}%</span>
      </div>`).join("") +
      `<div class="legend"><span><span class="sw acc"></span> concordant with a group-specific gene</span>` +
      `<span><span class="sw" style="width:3px;height:12px;border-radius:1px;background:var(--faint)"></span> group-shuffled control (chance)</span></div>`;
  }
  return { init };
})();
