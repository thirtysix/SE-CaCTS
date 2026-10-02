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
    const inferred = (m.cn_sources || {})["inferred from ChIP input"] || 0;
    const measured = src.filter(([k]) => k !== "inferred from ChIP input");
    const nMeasured = measured.reduce((t, [, v]) => t + v, 0);
    const card = (title, body, cls = "", id = "") =>
      `<section class="card about-card ${cls}"${id ? ` id="${id}"` : ""}><div class="card-h"><h3>${title}</h3></div><div class="card-b">${body}</div></section>`;
    const relBody = !rels.length ? "" : `
      <p>Each release is a full rebuild and rescore; the dashboard shows the newest. Counts are specific
      super-enhancers at permutation FDR ≤ 0.10.</p>
      <div class="scroll" style="margin:0 0 12px"><table class="tbl">
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

    U.el("about-body").innerHTML = `<div class="about-grid">
    ${card("What this is", `
      <p>A reference atlas of <b>lineage-specific super-enhancers</b> across ${fmt(m.n_lines)} cancer cell lines.
      Super-enhancers (SEs) are large clusters of active enhancers, marked by broad <b>H3K27ac</b>, that drive the
      genes defining a cell's identity. SE-CaCTS asks, for every SE, <em>how specific is it to one cancer lineage,
      disease or cell line?</em>, and reads out the SEs that most distinguish each group.</p>
      <ul>
        <li><b>SE atlas:</b> specific SEs per lineage and primary disease (calls), per subtype (rankings), and per
          cell line, compared with all lines and with the other lines of its lineage, disease and subtype.</li>
        <li><b>Genomic View (IGV):</b> any line in a genome browser (${U.link("igv", "igv.js")}), with those
          comparisons as tracks, its copy number and its H3K27ac coverage.</li>
        <li><b>SE finder</b>, <b>CN ablation</b> and <b>Concordance</b>: a per-gene lookup, what copy-number
          correction changes, and the cross-check against expression.</li>
      </ul>`, "span2")}

    ${card("Pipeline", `<div class="pipeline-card">${pipelineSVG(m)}</div>
      <ul>
        <li><b>Baseline experiments only.</b> Each experiment's ${U.link("geo", "GEO")} record was read and classified as untreated,
          control arm, perturbed, or not H3K27ac (by a language model, which agreed with a reference labelling of 461
          experiments on 98% of keep-or-drop decisions). Drug-treated, knocked-down and otherwise perturbed
          experiments, input controls and other marks are removed.</li>
        <li><b>SE calling with ${U.link("cnrose", "<code>cnrose</code>")}</b>, a bigWig-native, copy-number-aware
          reimplementation of ${U.link("rose", "ROSE")}, validated bit-for-bit against ${U.link("rose2", "ROSE2")}.</li>
        <li><b>Union atlas.</b> Per-sample SE calls are merged into one catalogue (≥ 25% reciprocal overlap) of
          <b>${fmt(m.n_ses)} SE loci</b>; signal is quantified for every locus in every sample, normalized across
          studies (${U.link("s3norm", "S3norm")}) behind a QC gate, and collapsed from ${fmt(m.n_samples)} experiments to
          <b>${fmt(m.n_lines)} cell lines</b>.</li>
        <li><b>Specificity (CaCTS JSD).</b> For each group, every SE's signal profile is scored by Jensen-Shannon
          divergence against a perfectly group-specific profile (via ${U.link("pycacts", "<code>pyCaCTS</code>")});
          lower is more specific.</li>
        <li><b>Significance</b> from a label-permutation FDR (see Statistics).</li>
      </ul>`, "span2")}

    ${card("Data", `<ul>
        <li><b>H3K27ac ChIP-seq:</b> ${U.link("chipatlas", "ChIP-Atlas")} (hg38), plus ${U.link("sra", "SRA")}
          experiments it lacks run through its own v1 pipeline (identical outputs on validation samples). ${fmt(m.n_pull)} untreated or control experiments
          pulled, ${fmt(m.n_samples)} passing the ≥ 2,000-peak QC gate.</li>
        <li><b>Cell-line labels:</b> ${U.link("depmap", "DepMap")} 2026q1 <code>Model.csv</code>
          (${U.link("oncotree", "OncoTree")} lineage, disease, subtype); lines outside DepMap via
          ${U.link("cellosaurus", "Cellosaurus")} and ${U.link("ncit", "NCIt")}.</li>
        <li><b>Copy number:</b> ${measured.map(([k, v]) => `${k} for ${v}`).join(", ")} lines (${fmt(nMeasured)} measured;
          sources: ${U.link("depmap", "DepMap")}, ${U.link("depmap_ref", "Tsherniak <i>et al.</i> 2017")};
          ${U.link("cmp", "Cell Model Passports")}, ${U.link("cmp_ref", "van der Meer <i>et al.</i> 2019")};
          ${U.link("ccle", "CCLE")}, ${U.link("ccle_ref", "Ghandi <i>et al.</i> 2019")});
          ${inferred ? `inferred from the experiments' own input controls for ${fmt(inferred)} lines no measured source
          covers, a weaker correction (the <b>Measured copy number only</b> analysis leaves them out).` : ""}</li>
        <li><b>Expression:</b> ${U.link("depmap", "DepMap")} protein-coding expression, for the concordance layer.</li>
        <li><b>Gene coordinates:</b> ${U.link("ensembl", "Ensembl GRCh38.106")}. <b>Engines:</b>
          ${U.link("cnrose", "<code>cnrose</code>")} (SE calling), ${U.link("pycacts", "<code>pyCaCTS</code>")} (JSD
          specificity and the permutation null). <b>Genome browser:</b> ${U.link("igv", "igv.js")}.</li>
      </ul>`)}

    ${card("Statistics", `
      <p><b>Permutation FDR.</b> The null is <b>measured, not assumed</b>: shuffle which cell line carries which
      group label (group sizes kept), recompute the score, repeat 1,000 times, and take a Benjamini-Hochberg FDR
      against that null. On shuffled labels, where nothing real exists, it makes ${fmt(cal.perm_shuffled_calls)}
      calls; the normal-approximation null it replaced called ${cal.analytic_shuffled_pct}% of tests.</p>
      <p><b>Copy number</b> is divided out at scoring time. At the group level this mostly <b>rescues</b> real,
      copy-neutral specificity that amplicon variance was hiding in the null (see CN ablation).</p>
      <p><b>What counts as a specific super-enhancer.</b> Every locus that any experiment calls a super-enhancer is
      scored in every line, so a group can have the most H3K27ac at a locus that none of its own experiments calls a
      super-enhancer. Since release v3.0.1 a call needs both: the group's signal is specific (permutation FDR ≤ 0.10)
      <em>and</em> at least one experiment of the group (lineage, disease, subtype, or the line itself) calls a
      super-enhancer overlapping the locus. The rule removed 18% of lineage calls and 74% of per-line calls; FDRs are
      those computed over every locus.</p>
      <p><b>Levels.</b> Calls at <b>lineage</b> and <b>primary disease</b>. <b>Subtypes</b> are rankings only:
      ${m.n_subtypes_single} of ${m.n_subtypes} hold a single line and ${m.n_subtypes_le4} hold four or fewer.
      <b>Cell lines</b> are scored by what their independent studies agree on (the value at least 75% of them
      reach) and compared four ways: with all lines, and with the other lines of the same lineage, disease and
      subtype (shuffling labels only among those relatives, so an SE shared across the group is not called). A
      comparison is tested only with at least two independent studies and at least four lines in the group;
      shuffled labels give no calls in any of the four.</p>
      <p><b>Analysis selector.</b> The lineage and disease calls can be switched to the atlas scored without the
      lines whose copy number is inferred, or without copy-number correction at all.</p>`)}

    ${card("Reading the results", `<ul>
        <li><b>One table per cell line.</b> Rows are the SEs specific vs all lines, ranked by one score; each
          <b>FDR</b> column says whether the SE also stands out against that comparison (<b>&gt; 0.10</b> tested, not
          called; <b>–</b> not tested). <b>Called vs</b> lists the comparisons that call it: a lineage programme
          reads A·L without D.</li>
        <li><b>Called</b> (x/N experiments) says whether the line's own experiments call an SE at the locus:
          "specific" means specific H3K27ac signal, which is not always an SE in that line.</li>
        <li><b>Genes within 100 kb, not a target call.</b> Each SE lists every protein-coding gene whose body lies
          within 100 kb of it, nearest first, with the distance (<b>+N</b> shows the rest; the filter and the SE finder
          match any of them). Where none lies within 100 kb, the nearest gene is shown and marked. Proximity is not a
          scored link: the target may be further away (long-range enhancers act over a megabase) or not the nearest,
          and the strongest lymphoid SEs sit at the immunoglobulin loci, whose neighbouring protein-coding genes (e.g.
          TMEM121) are bystanders. <span class="conc-badge" style="margin:0">⇌</span> marks a gene that is itself
          specific to the group in DepMap expression; the Concordance tab is the aggregate check.</li>
        <li><b>Rows near one gene often tile one SE domain</b> (↳ tiles #N): count domains, not rows.</li>
        <li><span class="flag-chip" style="margin:0">⚠</span> marks known artifact classes: loci on chrY (presence
          follows the line's sex) and copy number below 0.3 (deep deletions), removed in the next release.</li>
      </ul>`)}

    ${card("What is NOT claimed", `<ul>
        <li>Any count from the analytic (normal-approximation) null.</li>
        <li>Subtype-level calls: their number depends on how the multiple-testing correction is pooled.</li>
        <li>That a line-level specific locus is a super-enhancer in that line (see <b>Called</b>).</li>
        <li>That an SE regulates any listed gene (nearest or not), or that each row is an independent element.</li>
        <li>Loci flagged <span class="flag-chip" style="margin:0">⚠</span>.</li>
      </ul>`)}

    ${relBody ? card("Releases", relBody, "span2", "releases") : ""}

    ${card("Credit", `<p>Method: <b>CaCTS</b> (${U.link("cacts", "Reddy <em>et al.</em>, <em>Sci. Adv.</em> 2021")};
      code ${U.link("cactscode", "lawrenson-lab/CaCTS")}), adapted here from gene expression to the super-enhancer
      layer. SE calling reimplements ${U.link("rose", "ROSE")} / ${U.link("rose2", "ROSE2")}; cross-study normalization
      by ${U.link("s3norm", "S3norm")}; genome browser ${U.link("igv", "igv.js")}. Sister project:
      ${U.link("pycacts", "<b>pyCaCTS</b>")} (the master-transcription-factor atlas on the same JSD engine).
      Data: ${U.link("chipatlas", "ChIP-Atlas")}, ${U.link("geo", "GEO")} / ${U.link("sra", "SRA")},
      ${U.link("depmap", "DepMap")}, ${U.link("cmp", "Cell Model Passports")}, ${U.link("ccle", "CCLE")},
      ${U.link("cellosaurus", "Cellosaurus")}, ${U.link("ncit", "NCIt")}, ${U.link("oncotree", "OncoTree")},
      ${U.link("ensembl", "Ensembl")}. Source: ${U.link("repo", "github.com/thirtysix/SE-CaCTS")}.</p>
      <p style="margin:10px 0 0"><b>Cell-line data references</b></p>
      <ul>
        <li>${U.link("depmap", "DepMap")} (release 26Q1: cell-line labels, WGS and WES copy number, expression).
          Tsherniak A <i>et al.</i> Defining a cancer dependency map. <i>Cell</i> 2017;170:564–576.
          ${U.link("depmap_ref", "doi:10.1016/j.cell.2017.06.010")}</li>
        <li>${U.link("cmp", "Cell Model Passports")} (WES copy number). van der Meer D <i>et al.</i> Cell Model
          Passports: a hub for clinical, genetic and functional datasets of preclinical cancer models.
          <i>Nucleic Acids Res</i> 2019;47:D923–D929. ${U.link("cmp_ref", "doi:10.1093/nar/gky872")}</li>
        <li>${U.link("ccle", "CCLE")} (SNP6 copy-number segments, via ${U.link("cbioportal", "cBioPortal")}). Ghandi M
          <i>et al.</i> Next-generation characterization of the Cancer Cell Line Encyclopedia. <i>Nature</i>
          2019;569:503–508. ${U.link("ccle_ref", "doi:10.1038/s41586-019-1186-3")}</li>
      </ul>`, "span2")}
    </div>`;
  }
  return { init };
})();
