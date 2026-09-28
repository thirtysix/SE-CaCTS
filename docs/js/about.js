/* about.js — method, data provenance, and the explicit "what is NOT claimed" section (from RESULTS.md). */
const About = (() => {
  async function init() {
    const m = await DataLoader.loadJSON("data/meta.json");
    const rels = await DataLoader.loadJSON("data/releases.json").catch(() => []);
    const f = n => (n == null ? "—" : Number(n).toLocaleString());
    const relTable = !rels.length ? "" : `
      <h3 id="releases">Releases</h3>
      <p>Each release is a full rebuild and rescore; the dashboard always shows the newest. Counts are
      specific super-enhancers at permutation FDR ≤ 0.10.</p>
      <div class="card scroll" style="margin:0 0 12px"><table class="tbl">
        <thead><tr><th>release</th><th>date</th><th class="num">cell lines</th><th class="num">samples</th>
          <th class="num">SE loci</th><th class="num">lineage calls</th><th class="num">disease calls</th>
          <th>copy number</th></tr></thead>
        <tbody>${rels.slice().reverse().map(r => `<tr>
          <td><b>${U.esc(r.version)}</b></td><td class="mono" style="white-space:nowrap">${U.esc(r.date)}</td>
          <td class="num mono">${f(r.n_lines)}</td><td class="num mono">${f(r.n_samples)}</td>
          <td class="num mono">${f(r.n_ses)}</td><td class="num mono">${f(r.n_lineage_calls)}</td>
          <td class="num mono">${f(r.n_disease_calls)}</td>
          <td>${Object.entries(r.cn_sources || {}).filter(([, v]) => v > 0).map(([k, v]) => `${v} ${U.esc(k)}`).join(", ")}</td>
        </tr>`).join("")}</tbody></table></div>
      <ul>${rels.slice().reverse().map(r => `<li><b>${U.esc(r.version)}</b> (${U.esc(r.date)})${r.title ? ` — ${U.esc(r.title)}` : ""}.
        ${U.esc(r.notes || "")}</li>`).join("")}</ul>`;
    const src = Object.entries(m.cn_sources || {}).filter(([, v]) => v > 0);
    const cal = m.calibration || { analytic_shuffled_pct: 6.05, perm_shuffled_calls: 0 };
    U.el("about-body").innerHTML = `
      <h3>What this is</h3>
      <p>A reference atlas of <b>lineage-specific super-enhancers</b> across the cancer cell-line panel.
      Super-enhancers (SEs) are large clusters of active enhancers, marked by broad <b>H3K27ac</b>, that
      drive the genes defining a cell's identity. SE-CaCTS asks, for every SE, <em>how specific is it to one
      cancer lineage or disease?</em> — and reads out the SEs that most distinguish each group.</p>

      <h3>Pipeline, end to end</h3>
      <ul>
        <li><b>SE calling — <code>cnrose</code>.</b> Each H3K27ac experiment is called for super-enhancers
          with a bigWig-native, copy-number-aware reimplementation of ROSE, validated bit-for-bit against
          ROSE2. It stitches nearby enhancer peaks and applies the tangent-cutoff rule, but reads signal
          straight from the coverage track and can divide out copy number.</li>
        <li><b>Union atlas — reduce.</b> Per-sample SE calls are merged into one catalogue (≥ 25% reciprocal
          overlap) of <b>${m.n_ses.toLocaleString()} SE loci</b>, and signal is quantified for every locus in
          every sample, cross-study-normalized (S3norm) with a QC gate, and collapsed from
          ${m.n_samples.toLocaleString()} experiments to <b>${m.n_lines} cancer cell lines</b>.</li>
        <li><b>Specificity — CaCTS JSD.</b> For each Oncotree group, the per-group mean signal of every SE is
          scored by Jensen–Shannon divergence against a perfectly group-specific profile (via
          <code>pyCaCTS</code>). Lower score = more group-specific.</li>
        <li><b>Significance — a label-permutation FDR</b> (below).</li>
      </ul>

      <h3>The specificity call — a permutation FDR</h3>
      <p>For each group we ask which SEs are significantly more concentrated there than elsewhere. The null is
      <b>measured, not assumed</b>: shuffle which cell line carries which group label, recompute the JSD,
      repeat 1,000×, and take a Benjamini–Hochberg FDR against that empirical null. This replaced a
      normal-approximation null that <b>failed calibration outright</b> — run on shuffled labels, where
      nothing real exists to find, it called ${cal.analytic_shuffled_pct}% of tests "specific"; the permutation null
      made ${cal.perm_shuffled_calls.toLocaleString()} calls.</p>
      <div class="callout">Copy number is corrected at scoring time, from measured copy number for every line
      (${src.map(([k, v]) => `${k} for ${v}`).join(", ")}). Correction is not merely a
      penalty on amplified signal: at the group level it mainly <b>rescues</b> real, copy-neutral specificity
      that amplicon variance was masking in the permutation null — see the CN ablation tab.</div>

      <h3>Resolution — read this first</h3>
      <p>The panel supports specificity <b>calls only at the lineage and primary-disease levels</b>. At the
      subtype and single-cell-line levels the atlas shows <b>rankings only</b>: ${m.n_subtypes_single} of
      ${m.n_subtypes} subtypes contain a single cell line and ${m.n_subtypes_le4} contain ≤ 4, and because the permutation preserves group size, a random handful
      of lines is as "specific" as the real grouping. Rankings (which SE is most concentrated in a group)
      stay meaningful there; significance calls do not.</p>

      <h3>SE → gene links, and why one gene can appear many times</h3>
      <p>Each SE is labelled with its nearest protein-coding gene. That is a <b>locational label, not a scored
      regulatory association</b>. Two consequences to keep in mind when reading the atlas:</p>
      <ul>
        <li><b>Distinct rows near one gene usually tile a single SE domain.</b> A gene can appear several
          times in a group's top SEs — those are different loci (check the length and coordinates), but they
          often belong to one super-enhancer region: a large stitched domain plus the smaller peaks nested
          inside it (the catalogue keeps them separate when their reciprocal overlap is under 25%). Count SE
          <em>domains</em>, not rows.</li>
        <li><b>The nearest protein-coding gene can be a bystander.</b> The strongest lymphoid SEs sit at the
          immunoglobulin loci (IGH on chr14q32, IGL on chr22); their nearest annotated protein-coding genes
          (e.g. TMEM121) are not the functional target — the immunoglobulin genes themselves are not in the
          protein-coding annotation set.</li>
      </ul>
      <p>The <b>Concordance</b> tab is the aggregate cross-check on these links: SE-proximal genes are
      themselves group-specific in expression far above chance, and the effect decays with distance.</p>

      <h3>Data sources</h3>
      <ul>
        <li><b>H3K27ac ChIP-seq</b> — ChIP-Atlas (hg38); ${m.n_pull.toLocaleString()} experiments pulled,
          ${m.n_samples.toLocaleString()} passing the ≥ 2,000-peak QC gate.</li>
        <li><b>Cell-line annotation & copy number</b> — DepMap 2026q1: <code>Model.csv</code> (Oncotree
          lineage / disease / subtype), <code>OmicsCNGeneWGS.csv</code> (WGS gene-level copy number, the
          preferred source), <code>OmicsCNGeneMC_WES.csv</code> (WES), and the protein-coding expression
          matrix (used for the concordance layer). Cell Model Passports WES pureCN (2025) copy number for
          lines DepMap WGS does not cover. Lines outside DepMap take Oncotree labels via Cellosaurus → NCIt.</li>
        <li><b>Gene coordinates</b> — Ensembl GRCh38.106.</li>
        <li><b>Engines</b> — <code>cnrose</code> (SE calling), <code>pyCaCTS</code> (JSD specificity + the
          permutation null).</li>
      </ul>

      <h3>What is NOT claimed</h3>
      <ul>
        <li>Any specific-SE <b>count</b> from the analytic (normal-approximation) null.</li>
        <li>Subtype- or cell-line-level specificity <b>calls</b>.</li>
        <li>That a given SE regulates its nearest gene (proximity annotation only).</li>
        <li>That each row is an independent regulatory element (nested / tiling loci — see above).</li>
      </ul>

      ${relTable}

      <h3>Credit</h3>
      <p>Method: <b>CaCTS</b> — Reddy <em>et al.</em>, <em>Sci. Adv.</em> 2021 —
      <a href="https://github.com/lawrenson-lab/CaCTS" target="_blank" rel="noopener">lawrenson-lab/CaCTS</a>,
      here adapted from gene expression to the super-enhancer layer. SE calling reimplements ROSE2. Sister
      project: <b>pyCaCTS</b> (the master-transcription-factor atlas on the same JSD engine). Data: ChIP-Atlas
      H3K27ac, DepMap, Ensembl.</p>`;
  }
  return { init };
})();
