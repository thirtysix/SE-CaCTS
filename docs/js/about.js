/* about.js: what this is, the releases, data sources, the pipeline flow chart, then method details and the
   explicit "what is NOT claimed" section (from RESULTS.md). The flow chart is an inline SVG built from the
   release's own numbers (meta.json), styled with the page's theme tokens; its <style> carries light-theme
   fallbacks, so the same SVG exported on its own (phase2/figures/export_pipeline.py) is the poster figure. */
const About = (() => {
  const fmt = n => (n == null ? "n/a" : Number(n).toLocaleString());

  function pipelineSVG(m) {
    const cn = Object.entries(m.cn_sources || {}).filter(([, v]) => v > 0).map(([k, v]) => `${k} ${v}`).join(" · ");
    const inferred = (m.cn_sources || {})["inferred from ChIP input"] || 0;
    const measured = Object.entries(m.cn_sources || {}).filter(([k, v]) => v > 0 && k !== "inferred from ChIP input")
      .map(([k, v]) => `${k} ${v}`).join(" · ");
    const steps = [
      ["H3K27ac ChIP-seq", m.pull_desc || [`ChIP-Atlas, hg38: ${fmt(m.n_pull)} experiments on human cancer cell lines`]],
      ["Super-enhancer calling", ["cnrose: ROSE re-implemented on the bigWig coverage", "(calls identical to ROSE2)"]],
      ["Union catalogue", [`${fmt(m.n_ses)} super-enhancer loci; exact signal for every locus in every sample`]],
      ["Cross-study normalization", ["S3norm, behind a QC gate of 2,000 peaks:", `${fmt(m.n_samples)} experiments from ${fmt(m.n_lines)} cell lines`]],
      ["Copy-number correction", inferred ? [measured, `and ${inferred} lines with copy number inferred from ChIP input`]
                                          : ["measured copy number for every line:", cn]],
      ["Specificity score", ["CaCTS Jensen-Shannon divergence at each Oncotree level:", "lineage, primary disease, subtype, cell line"]],
      ["Significance", ["label-permutation FDR (1,000 shuffles, group sizes kept);", `${fmt((m.calibration || {}).perm_shuffled_calls)} calls on shuffled labels`]],
      ["Specific super-enhancers", [`${fmt(m.n_lineage_calls)} lineage and ${fmt(m.n_disease_calls)} disease calls (FDR ≤ 0.10),`, "checked by CN ablation and expression concordance"]],
    ];
    const inputs = [                                   // [target step index, title, subtitle]
      m.pull_desc ? [0, "ChIP-Atlas, SRA, GEO", "coverage, reads, sample metadata"] : [0, "ChIP-Atlas", "bigWig coverage + peaks"],
      (m.cn_sources || {})["inferred from ChIP input"] ? [4, "DepMap, CMP, CCLE", "WGS / WES / SNP6; ChIP input"]
                                                         : [4, "DepMap, Cell Model Passports", "WGS / WES copy number"],
      [5, "DepMap Model, Cellosaurus", "Oncotree labels per line"],
      [7, "DepMap RNA-seq", "expression, for validation"],
    ];
    const W = 780, BX = 16, BW = 470, IX = 540, IW = 224, BH = 62, GAP = 24, TOP = 12;
    const y = i => TOP + i * (BH + GAP);
    const H = y(steps.length - 1) + BH + TOP;
    const esc = U.esc;
    const box = (s, i) => {
      const [t, lines] = s, yy = y(i);
      return `<g>
        <rect class="pl-box" x="${BX}" y="${yy}" width="${BW}" height="${BH}" rx="10"/>
        <rect class="pl-bar" x="${BX}" y="${yy}" width="5" height="${BH}" rx="2"/>
        <circle class="pl-num" cx="${BX + 28}" cy="${yy + BH / 2}" r="13"/>
        <text class="pl-numt" x="${BX + 28}" y="${yy + BH / 2 + 4.5}" text-anchor="middle">${i + 1}</text>
        <text class="pl-t" x="${BX + 52}" y="${yy + (lines.length > 1 ? 21 : 25)}">${esc(t)}</text>
        ${lines.map((l, k) => `<text class="pl-s" x="${BX + 52}" y="${yy + (lines.length > 1 ? 38 : 44) + k * 16}">${esc(l)}</text>`).join("")}
      </g>`;
    };
    const arrow = i => `<line class="pl-arrow" x1="${BX + BW / 2}" y1="${y(i) + BH}" x2="${BX + BW / 2}" y2="${y(i + 1) - 3}" marker-end="url(#pl-head)"/>`;
    const input = ([i, t, s]) => {
      const yy = y(i) + 8, h = BH - 16;
      return `<g>
        <rect class="pl-in" x="${IX}" y="${yy}" width="${IW}" height="${h}" rx="8"/>
        <text class="pl-it" x="${IX + 12}" y="${yy + 19}">${esc(t)}</text>
        <text class="pl-is" x="${IX + 12}" y="${yy + 35}">${esc(s)}</text>
        <line class="pl-arrow" x1="${IX}" y1="${yy + h / 2}" x2="${BX + BW + 4}" y2="${yy + h / 2}" marker-end="url(#pl-head)"/>
      </g>`;
    };
    return `<svg id="pipeline" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="100%" role="img"
      aria-label="SE-CaCTS pipeline, from ChIP-Atlas H3K27ac to specific super-enhancers">
      <style>
        #pipeline text{font-family:var(--pl-font, "Liberation Sans", Arial, Helvetica, sans-serif)}
        #pipeline .pl-box{fill:var(--surface,#ffffff); stroke:var(--border,#d4dcdf); stroke-width:1}
        #pipeline .pl-bar{fill:var(--accent,#008a7e)}
        #pipeline .pl-num{fill:var(--accent-soft,#e1f1ef)}
        #pipeline .pl-numt{fill:var(--accent-deep,#0a615c); font-size:13px; font-weight:700}
        #pipeline .pl-t{fill:var(--text,#12222a); font-size:15px; font-weight:700}
        #pipeline .pl-s{fill:var(--muted,#5a6b73); font-size:12.5px}
        #pipeline .pl-in{fill:var(--bg-alt,#f6f9f9); stroke:var(--border,#d4dcdf); stroke-width:1}
        #pipeline .pl-it{fill:var(--text,#12222a); font-size:12.5px; font-weight:700}
        #pipeline .pl-is{fill:var(--muted,#5a6b73); font-size:11.5px}
        #pipeline .pl-arrow{stroke:var(--faint,#87969d); stroke-width:1.6; fill:none}
        #pipeline .pl-head{fill:var(--faint,#87969d)}
      </style>
      <defs><marker id="pl-head" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
        <path class="pl-head" d="M0,0 L10,5 L0,10 z"/></marker></defs>
      ${steps.slice(0, -1).map((_, i) => arrow(i)).join("")}
      ${steps.map(box).join("")}
      ${inputs.map(input).join("")}
    </svg>`;
  }

  async function init() {
    const m = await DataLoader.loadJSON("data/meta.json");
    const rels = await DataLoader.loadJSON("data/releases.json").catch(() => []);
    const src = Object.entries(m.cn_sources || {}).filter(([, v]) => v > 0);
    const cal = m.calibration || { analytic_shuffled_pct: 6.05, perm_shuffled_calls: 0 };
    const relTable = !rels.length ? "" : `
      <h3 id="releases">Releases</h3>
      <p>Each release is a full rebuild and rescore, and the dashboard always shows the newest. Counts are
      specific super-enhancers at permutation FDR ≤ 0.10.</p>
      <div class="card scroll" style="margin:0 0 12px"><table class="tbl">
        <thead><tr><th>release</th><th>date</th><th class="num">cell lines</th><th class="num">samples</th>
          <th class="num">SE loci</th><th class="num">lineage calls</th><th class="num">disease calls</th>
          <th>copy number</th></tr></thead>
        <tbody>${rels.slice().reverse().map(r => `<tr>
          <td><b>${U.esc(r.version)}</b></td><td class="mono" style="white-space:nowrap">${U.esc(r.date)}</td>
          <td class="num mono">${fmt(r.n_lines)}</td><td class="num mono">${fmt(r.n_samples)}</td>
          <td class="num mono">${fmt(r.n_ses)}</td><td class="num mono">${fmt(r.n_lineage_calls)}</td>
          <td class="num mono">${fmt(r.n_disease_calls)}</td>
          <td>${Object.entries(r.cn_sources || {}).filter(([, v]) => v > 0).map(([k, v]) => `${v} ${U.esc(k)}`).join(", ")}</td>
        </tr>`).join("")}</tbody></table></div>
      <ul>${rels.slice().reverse().map(r => `<li><b>${U.esc(r.version)}</b> (${U.esc(r.date)})${r.title ? `: ${U.esc(r.title)}` : ""}.
        ${U.esc(r.notes || "")}</li>`).join("")}</ul>`;

    U.el("about-body").innerHTML = `
      <h3>What this is</h3>
      <p>A reference atlas of <b>lineage-specific super-enhancers</b> across the cancer cell-line panel.
      Super-enhancers (SEs) are large clusters of active enhancers, marked by broad <b>H3K27ac</b>, that
      drive the genes defining a cell's identity. SE-CaCTS asks, for every SE, <em>how specific is it to one
      cancer lineage or disease?</em>, and reads out the SEs that most distinguish each group.</p>

      ${relTable}

      <h3>Data sources</h3>
      <ul>
        <li><b>H3K27ac ChIP-seq:</b> ChIP-Atlas (hg38); ${fmt(m.n_pull)} experiments pulled,
          ${fmt(m.n_samples)} passing the ≥ 2,000-peak QC gate.</li>
        <li><b>Cell-line annotation and copy number:</b> DepMap 2026q1 <code>Model.csv</code> (Oncotree
          lineage, disease and subtype), <code>OmicsCNGeneWGS.csv</code> (WGS gene-level copy number, the
          preferred source), <code>OmicsCNGeneMC_WES.csv</code> (WES), and the protein-coding expression
          matrix (used for the concordance layer). Cell Model Passports WES pureCN (2025) copy number for
          lines DepMap WGS does not cover. Lines outside DepMap take Oncotree labels via Cellosaurus and NCIt.</li>
        <li><b>Gene coordinates:</b> Ensembl GRCh38.106.</li>
        <li><b>Engines:</b> <code>cnrose</code> (SE calling), <code>pyCaCTS</code> (JSD specificity and the
          permutation null).</li>
      </ul>

      <h3>Pipeline</h3>
      <div class="card pipeline-card">${pipelineSVG(m)}</div>

      <h3>Pipeline, step by step</h3>
      <ul>
        <li><b>SE calling with <code>cnrose</code>.</b> Each H3K27ac experiment is called for super-enhancers
          with a bigWig-native, copy-number-aware reimplementation of ROSE, validated bit-for-bit against
          ROSE2. It stitches nearby enhancer peaks and applies the tangent-cutoff rule, but reads signal
          straight from the coverage track and can divide out copy number.</li>
        <li><b>Union atlas.</b> Per-sample SE calls are merged into one catalogue (≥ 25% reciprocal
          overlap) of <b>${fmt(m.n_ses)} SE loci</b>, and signal is quantified for every locus in
          every sample, normalized across studies (S3norm) behind a QC gate, and collapsed from
          ${fmt(m.n_samples)} experiments to <b>${fmt(m.n_lines)} cancer cell lines</b>.</li>
        <li><b>Specificity (CaCTS JSD).</b> For each Oncotree group, the per-group mean signal of every SE is
          scored by Jensen-Shannon divergence against a perfectly group-specific profile (via
          <code>pyCaCTS</code>). Lower score means more group-specific.</li>
        <li><b>Significance</b> comes from a label-permutation FDR (next section).</li>
      </ul>

      <h3>The specificity call: a permutation FDR</h3>
      <p>For each group we ask which SEs are significantly more concentrated there than elsewhere. The null is
      <b>measured, not assumed</b>: shuffle which cell line carries which group label, recompute the JSD,
      repeat 1,000 times, and take a Benjamini-Hochberg FDR against that empirical null. This replaced a
      normal-approximation null that <b>failed calibration outright</b>: run on shuffled labels, where
      nothing real exists to find, it called ${cal.analytic_shuffled_pct}% of tests "specific", while the
      permutation null made ${fmt(cal.perm_shuffled_calls)} calls.</p>
      <div class="callout">Copy number is corrected at scoring time, from measured copy number for every line
      (${src.map(([k, v]) => `${k} for ${v}`).join(", ")}). Correction is not merely a penalty on amplified
      signal. At the group level it mainly <b>rescues</b> real, copy-neutral specificity that amplicon
      variance was masking in the permutation null (see the CN ablation tab).</div>

      <h3>Resolution: read this first</h3>
      <p>The panel supports specificity <b>calls at the lineage and primary-disease levels</b>. At the
      subtype level the atlas shows <b>rankings only</b>: ${m.n_subtypes_single} of ${m.n_subtypes} subtypes
      contain a single cell line and ${m.n_subtypes_le4} contain four or fewer, and because the permutation
      preserves group size, a random handful of lines is as "specific" as the real grouping. Rankings (which SE
      is most concentrated in a group) stay meaningful there; significance calls do not.</p>
      <p><b>Cell lines</b> are scored by what their independent studies agree on (a consensus: the value at least
      75% of the studies reach), and each line is compared four ways: with <b>all lines</b>, and with the other
      lines of its <b>subtype</b>, <b>primary disease</b> and <b>lineage</b> (the null shuffles labels only among
      those relatives, so a super-enhancer shared across the group is not called). A comparison is tested only
      when the line has at least two independent studies and the group holds at least four lines; otherwise the
      list is a ranking. Shuffled labels give no calls in any of the four comparisons. "Specific" means specific
      H3K27ac signal: the <b>Called</b> column says whether the line's own experiments call a super-enhancer there.</p>

      <h3>SE to gene links, and why one gene can appear many times</h3>
      <p>Each SE is labelled with its nearest protein-coding gene. That is a <b>locational label, not a scored
      regulatory association</b>. Two consequences to keep in mind when reading the atlas:</p>
      <ul>
        <li><b>Distinct rows near one gene usually tile a single SE domain.</b> A gene can appear several
          times in a group's top SEs. Those are different loci (check the length and coordinates), but they
          often belong to one super-enhancer region: a large stitched domain plus the smaller peaks nested
          inside it (the catalogue keeps them separate when their reciprocal overlap is under 25%). Count SE
          <em>domains</em>, not rows.</li>
        <li><b>The nearest protein-coding gene can be a bystander.</b> The strongest lymphoid SEs sit at the
          immunoglobulin loci (IGH on chr14q32, IGL on chr22); their nearest annotated protein-coding genes
          (e.g. TMEM121) are not the functional target, because the immunoglobulin genes themselves are not in
          the protein-coding annotation set.</li>
      </ul>
      <p>The <b>Concordance</b> tab is the aggregate cross-check on these links: SE-proximal genes are
      themselves group-specific in expression far above chance, and the effect decays with distance.</p>

      <h3>What is NOT claimed</h3>
      <ul>
        <li>Any specific-SE <b>count</b> from the analytic (normal-approximation) null.</li>
        <li>Subtype- or cell-line-level specificity <b>calls</b>.</li>
        <li>That a given SE regulates its nearest gene (proximity annotation only).</li>
        <li>That each row is an independent regulatory element (nested and tiling loci, see above).</li>
      </ul>

      <h3>Credit</h3>
      <p>Method: <b>CaCTS</b> (Reddy <em>et al.</em>, <em>Sci. Adv.</em> 2021;
      <a href="https://github.com/lawrenson-lab/CaCTS" target="_blank" rel="noopener">lawrenson-lab/CaCTS</a>),
      here adapted from gene expression to the super-enhancer layer. SE calling reimplements ROSE2. Sister
      project: <b>pyCaCTS</b> (the master-transcription-factor atlas on the same JSD engine). Data: ChIP-Atlas
      H3K27ac, DepMap, Cell Model Passports, Ensembl.</p>`;
  }
  return { init };
})();
