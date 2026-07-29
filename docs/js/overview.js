/* overview.js — hero stats, the resolution caveat, and a guide to the tabs. */
const Overview = (() => {
  async function init() {
    const meta = await DataLoader.loadJSON("data/meta.json");

    U.el("snap").innerHTML = [
      [meta.n_ses.toLocaleString(), "super-enhancer loci", "union catalogue of SEs across the panel"],
      [meta.n_lines, "cancer cell lines", "DepMap cell lines with QC-passed H3K27ac"],
      [meta.n_samples.toLocaleString(), "H3K27ac samples", "individual experiments after the QC gate"],
      ["perm", "permutation-null FDR", "significance from a label-permutation null (B=1000)"],
    ].map(([b, l, t]) => `<div class="s" title="${U.esc(t)}"><b>${b}</b> ${l}</div>`).join("");

    U.el("ov-stats").innerHTML = [
      [meta.n_ses.toLocaleString(), "", "super-enhancer loci in the atlas", "the union catalogue of super-enhancers across all samples (≥25% reciprocal-overlap merge)"],
      [meta.n_lines, "", "cancer cell lines (282, DepMap-joined)", "distinct DepMap cell lines with H3K27ac data that passed QC; replicate experiments are collapsed to the line"],
      [meta.n_samples.toLocaleString(), "", "QC-passed H3K27ac experiments", "ChIP-Atlas H3K27ac experiments passing the ≥2,000-peak QC gate (of 2,916 pulled)"],
      ["0", "good", "false calls on shuffled labels<br>(the analytic null gave 6.05%)", "the calibration test: run the whole procedure on data whose group labels are shuffled, so nothing real exists. A working FDR calls ≈ nothing — the permutation null does; the normal-approximation null called 6.05%"],
    ].map(([k, cls, l, t]) => `<div class="stat" title="${U.esc(t)}"><div class="k ${cls}">${k}</div><div class="l">${l}</div></div>`).join("");

    // panel-at-a-glance: which resolutions carry CALLS vs rankings only
    U.el("ov-panel").innerHTML =
      `<span class="ps-cap" title="resolutions where a significance call is supported by the permutation null">Calls</span>` +
      [[meta.n_lineages, "lineages", "OncotreeLineage groups (e.g. Breast, Lung, Ovary/Fallopian Tube)"],
       [meta.n_diseases, "primary diseases", "OncotreePrimaryDisease groups (e.g. Invasive Breast Carcinoma)"]].map(([n, l, t]) =>
        `<div class="ps-item" title="${U.esc(t)}"><span class="ps-n">${n}</span><span class="ps-l">${l}</span></div>`).join("") +
      `<span class="ps-cap" style="margin-left:8px" title="finer resolutions — rankings are shown but significance is NOT called, because most groups have too few cell lines">Rankings only</span>` +
      [[meta.n_subtypes, "subtypes", "OncotreeSubtype groups — 29 of 75 contain a single cell line"],
       [meta.n_lines, "cell lines", "individual cell lines — permutation is degenerate for a single sample"]].map(([n, l, t]) =>
        `<div class="ps-item" title="${U.esc(t)}"><span class="ps-n">${n}</span><span class="ps-l">${l}</span></div>`).join("");

    U.el("ov-guide").innerHTML = [
      ["▦", "atlas", "SE atlas",
        "The core view. For any lineage or primary disease, the super-enhancers most specific to it — ranked by JSD, with the permutation FDR, the mean copy number at the locus, the nearest gene, the SE length and coordinates (each linked to the UCSC browser), and a ⇌ badge where the gene is <em>also</em> specific in expression. Filter by gene, tighten the FDR cutoff, or drop to subtype / cell-line rankings."],
      ["⊘", "cn", "CN ablation",
        "Why the copy-number layer earns its place. Scoring with vs without correction: the handful of recurrent-amplicon false calls it removes (MYCN, OTX2, FGFR2, ANO1) and the thousands of real, copy-neutral calls it <em>rescues</em> from amplicon noise in the null."],
      ["⇌", "concordance", "Concordance",
        "The cross-layer validation. Genes next to a group-specific super-enhancer are themselves group-specific in DepMap expression far above chance — and the concordance decays with SE→gene distance while a shuffled control stays flat, the signature of a local regulatory link."],
      ["⌕", "finder", "SE finder",
        "Look up a gene by symbol: every lineage and disease where a super-enhancer near it is group-specific, with its rank and cross-layer concordance mark."],
      ["ⓘ", "about", "About & methods",
        "How the atlas is built, the label-permutation FDR, the data sources, and an explicit list of what is deliberately <em>not</em> claimed."],
    ].map(([i, tab, t, s]) => `<a href="#${tab}" title="go to ${t}"><span class="gi">${i}</span><span class="gt"><b>${t}</b><small>${s}</small></span></a>`).join("");

    U.el("ov-credit").innerHTML = `
      <p class="cred-p">SE-CaCTS applies <b>CaCTS</b> (Cancer Core Transcription factor Specificity; Reddy
      <em>et al.</em>, Sci. Adv. 2021) — a Jensen–Shannon-divergence measure of how group-specifically a
      feature is active — to a new layer: the <b>super-enhancer</b>.</p>
      <p class="cred-p">It is the sister of <b>pyCaCTS</b>, which scores candidate master
      <b>transcription factors</b> from gene expression across the DepMap and TCGA panels. SE-CaCTS reuses
      pyCaCTS's JSD scoring engine and mirrors its dashboard design — the two are the same specificity idea
      on two layers of the same cells: the <b>genes</b> (pyCaCTS) and the <b>enhancers that drive them</b>
      (SE-CaCTS). The Concordance tab is where the two layers meet.</p>`;

    // Everything on this page is generated by the repo below — link it prominently, and point at the
    // two documents a reader should consult before quoting a number.
    const GH = "https://github.com/thirtysix/SE-CaCTS";
    const ghMark = `<svg viewBox="0 0 16 16" width="16" height="16" aria-hidden="true" fill="currentColor"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27s1.36.09 2 .27c1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.012 8.012 0 0 0 16 8c0-4.42-3.58-8-8-8Z"/></svg>`;
    U.el("ov-repo").innerHTML = `
      <a class="repo-cta" href="${GH}" target="_blank" rel="noopener"
         title="the full pipeline: SE calling, normalization, scoring, the permutation null, and this dashboard">
        ${ghMark}<span class="rc-t"><b>thirtysix/SE-CaCTS</b><small>the SE caller, the scoring pipeline, the permutation null, and this dashboard</small></span>
        <span class="rc-go">↗</span></a>
      <p class="cred-p" style="margin-top:14px">Open source and fully reproducible: the atlas artifacts,
      the per-group scores and the staging script that builds every file this page loads are all in the
      repository. Two documents are worth reading before quoting any number —
      <a href="${GH}/blob/main/RESULTS.md" target="_blank" rel="noopener"><b>RESULTS.md</b></a>, the claims
      document, which is deliberately narrower than the raw outputs; and
      <a href="${GH}/blob/main/README.md" target="_blank" rel="noopener"><b>README.md</b></a> for how the
      atlas was built. The two H3K27ac signal matrices exceed GitHub's file-size limit and are the one
      artifact not in the repository.</p>`;
  }
  return { init };
})();
