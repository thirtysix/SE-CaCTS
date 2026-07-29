# Step A — novelty check (result)

*Deep-research run, 2026-07-16. 3 search angles → 12 primary sources fetched → 48 claims extracted →
25 verified by 3-vote adversarial voting (24 confirmed 3–0, 1 refuted 1–2, 0 unverified). This file is the
durable evidence for ROADMAP Step A; `PRIOR_ART.md` carries the folded-in synthesis.*

## Verdict

**The claim "to our knowledge, the first copy-number-aware per-super-enhancer lineage-specificity score"
is still defensible as of 2026.** No published resource computes a *copy-number-corrected*, *per-SE*,
*lineage/cancer-type specificity* score from H3K27ac across a compendium. The two elements SE-CaCTS combines
(a per-SE JSD specificity score + copy-number correction of the H3K27ac signal) exist only **separately and
in different domains** — never jointly, and never even individually applied to super-enhancers.

This *resolves the one open item flagged as potentially project-killing* in the original PRIOR_ART.md
("SEdb 3.0 only partially read"). SEdb 3.0 is now fully checked and does not preempt the method.

## Q1 — SEdb 3.0 (2026): fully checked, no specificity metric, no CN correction

- **SEdb 3.0**, *Nucleic Acids Res* 2026, 54(D1):D322, DOI 10.1093/nar/gkaf1294 (PMC12807715; published
  2026-01-06). The largest SE compendium to date: **3,478,186 SEs from 5,387 H3K27ac ChIP-seq samples**
  across four species (human/mouse/Arabidopsis/maize; hg38/mm10/TAIR10/B73).
- Pipeline (verified against the NAR article **and** PMC full text): matched-input screening → ChIPQC →
  Bowtie2 v2.5.4 → MACS2 v2.2.9.1 (input as background) → **ROSE** SE calling. Verbatim: *"SE calling was
  then performed using the ROSE algorithm."*
- **No compendium-wide per-SE specificity statistic** (no tau, JSD, Shannon entropy, or prevalence) and
  **no copy-number / amplicon correction** — a full-text search for copy number / CNV / amplicon /
  amplification / ploidy returns **zero** hits. The two *new* 3.0 tools are **SE BLAST alignment** and
  **SE-driven core-TF enrichment**; neither is a specificity score. "Differential-overlapping-SE" is a
  pairwise comparison and "tissue-category" is a browsing filter, not a computed statistic.
- (No conflation risk: the Shannon-entropy database is **SEA / Super-Enhancer Archive**, a *different*
  resource — not SEdb.)

## Q2 — 2024–2026 CN-aware / lineage-specific-SE work: none removes CN to isolate specificity

- **Su, Peters, Soltys & Chan**, *BMC Genomics* 2025, 26:306, DOI 10.1186/s12864-025-11442-y (PMC11951689,
  PMID 40155863; bioRxiv 2024.04.11.588815). Confirms the CN-normalization recipe exactly: CNVkit WGS mode
  (`--method wgs --target-avg-size 50,000`) → for copy-**gain** peaks (log2 CNR > 0) divide counts by CNR,
  for copy-**loss** peaks multiply by CNR → "averaged signal per gene copy." **But** applied to general
  **ATAC-seq and G-quadruplex (BG4) ChIP-seq** peaks — *never super-enhancers, never H3K27ac* — and only two
  paired disease/control comparisons (Bloom-syndrome fibroblasts; Down-syndrome trisomy-21 lymphoblasts),
  **not a cancer/lineage compendium**, with **no specificity statistic**. It *removes* the CN–signal coupling
  (the direction SE-CaCTS wants). *(Note: baseline said "Su/Chen" — corresponding author is Chan; "Chen" was
  a typo. Recipe attribution is otherwise exact.)*
- **Zhang et al.**, *Nat Genet* 2016, 48:176–182, DOI 10.1038/ng.3470 (PMC4857881, PMID 26656844). The
  canonical amplicon-SE paper — deliberately **embraces** the CN↔H3K27ac coupling: intersects GISTIC focal
  amplifications (12 tumor types) with tissue-matched ROSE super-enhancers to *find* amplified driver SEs
  (KLF5, USP12, PARD6B, MYC-LASE, MYC-ECSE). No per-SE specificity statistic; no CN correction of H3K27ac.
  The **inverse** of SE-CaCTS.

## Q3 — JSD per-SE across a cancer compendium: not done

- **CaCTS** (Reddy et al., *Sci Adv* 2021, 7:eabf6123, DOI 10.1126/sciadv.abf6123; PMC8612691) is confirmed
  as the JSD template — *"an entropy-based measure of Jensen-Shannon divergence (JSD),"* score = −log10 JSD
  via the R `jsd` package — **but on TCGA TF RNA-seq**, not SEs; SEs appear only as post-hoc ROSE2/GSEA
  validation. It **embraces** copy number (candidate MTFs coincide with CN gain, P = 8.8×10⁻⁴), never
  correcting for it. So JSD has **never** been computed per-super-enhancer.
- **cSEAdb** (*PLOS Comput Biol* 2024, DOI 10.1371/journal.pcbi.1011873; PMC10883583) — closest cancer
  analog, but scores SE constituents by **presence/prevalence** (active-prevalence < ~0.21 = cell-specific;
  threshold set at the variation-of-information inflection point), **not JSD/tau/entropy**; normalization is
  DESeq2 RLE depth + genome-wide coverage scaling + zero-imputation — **no per-locus CN correction**.
  ⚠️ A WebSearch summarizer's gloss that the scaling step "accounts for copy number variations" is a
  **hallucination absent from the paper** — flagged and discarded during verification.
- **SEgene** (*npj Syst Biol Appl* 2025, DOI 10.1038/s41540-025-00533-x; PMC12089303) — peak-to-gene links
  (Pearson of H3K27ac CPM vs RNA-seq TPM within ±1 Mbp) and ranks SEs by **commonness** (the *inverse* of
  specificity); no per-SE JSD/tau/entropy, no CN. *(The one refuted claim in the run — "SEgene has no
  specificity metric of any kind" — was killed 1–2 because commonness is technically a frequency metric;
  this does not change the verdict.)*
- **Trends Genet** 2022 review (DOI 10.1016/j.tig.2022.11.001, PMID 35803787) quantifies SE relationships
  with sample-level **Jaccard co-occurrence** networks (MCL clustering), not a per-SE specificity statistic;
  corroborated by a *Brief Bioinform* 2025 survey (bbaf596): "No mentions of Jensen-Shannon divergence, tau,
  Shannon entropy, or specificity statistics."

## Q4 — amplicon-SE papers do not implicitly CN-correct

- **Kudo/Yamamoto et al.** (*Mol Syst Biol* 2025, DOI 10.1038/s44320-025-00098-1; PMC12130324; bioRxiv
  2024.11.11.622904) nominate MYB-SE in GI adenocarcinoma by **qualitative** cross-tissue H3K27ac
  inspection and **explicitly decline**: *"we cannot determine whether the H3K27ac signal of MYB-SE is
  associated with its amplification"* (WGS and H3K27ac from unmatched sources — the same matched-data
  problem SE-CaCTS must solve).
- **Watanabe et al.** (2024, PMC11462929, PMID 38923173) map the MYC BENC-CNA amplicon in B-ALL via strong
  H3K27ac/H3K4me1 but argue the **amplification creates** the SE signal — CN-correcting would erase their
  thesis. Another embrace-the-coupling case.

## Closest analogs — how each falls short

| analog | has per-SE specificity? | metric | CN-corrected? | compendium | gap vs SE-CaCTS |
|---|---|---|---|---|---|
| **SEdb 3.0** (NAR 2026) | No | — (calls only) | No | 5,387 H3K27ac samples | no specificity statistic at all; no CN |
| **cSEAdb** (PLOS CB 2024) | Yes (per-constituent) | presence/prevalence (~0.21) | No | NCI-60 / 28 types | non-JSD categorical metric; no CN; narrow panel |
| **CaCTS** (Sci Adv 2021) | Yes (per-TF) | −log10 JSD + null FDR | No (embraces CN) | TCGA 34 types | on TFs not SEs; no CN |
| **Su et al.** (BMC Gen 2025) | No | — | **Yes** (÷CNR) | 2 paired samples | CN recipe but on ATAC/G4 not SEs; no specificity |
| **Zhang et al.** (Nat Genet 2016) | No | — (intersection) | uses CN (opposite) | TCGA×H3K27ac, 12 types | *finds* amplified SEs; inverse goal |
| **SEgene** (npj SBA 2025) | commonness (inverse) | Pearson SE→gene | No | multi-sample | ranks common not specific; no CN |

## Caveats & residual open questions (carried into ROADMAP)

- **"First / to our knowledge" is bounded by what the harness indexed.** An unindexed preprint or niche
  resource could exist; none surfaced across SE databases, CN-normalization methods, amplicon-SE papers, and
  the JSD-specificity literature.
- **Live SEdb 3.0 website vs paper.** The verdict rests on the 2026 paper + PMC full text; a running-site
  feature not in the publication can't be fully excluded — worth a quick direct check of the live tool
  before writing "first" in a manuscript.
- **ecDNA / pediatric amplicon-SE subfield** (neuroblastoma MYCN, medulloblastoma) was **not** exhaustively
  covered — a possible edge case for implicit CN-correction.
- **Method soundness (new, non-prior-art):** assigning a single copy-number ratio to a ROSE-stitched,
  multi-constituent SE that spans heterogeneous CN segments is non-trivial — a real design question for
  Phase 3, not a novelty threat.
- **Validation data:** need matched WGS/CN + H3K27ac at cell-line scale (e.g. DepMap CN aligned to an
  H3K27ac compendium) to *empirically* show CN correction separates amplified loci (MYC) from genuine
  lineage specificity. (The ChIP-Atlas × local DepMap join is exactly this.)

## References

- SEdb 3.0 — *Nucleic Acids Res* 2026, 54(D1):D322 — https://academic.oup.com/nar/article/54/D1/D322/8373947 (DOI 10.1093/nar/gkaf1294; PMC12807715)
- cSEAdb — *PLOS Comput Biol* 2024 — https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011873 (PMC10883583)
- CaCTS (Reddy et al.) — *Sci Adv* 2021, 7:eabf6123 — https://www.science.org/doi/10.1126/sciadv.abf6123 (PMC8612691; github.com/lawrenson-lab/CaCTS)
- Su et al. — *BMC Genomics* 2025, 26:306 — https://link.springer.com/article/10.1186/s12864-025-11442-y (PMC11951689; PMID 40155863)
- Zhang et al. — *Nat Genet* 2016, 48:176 — https://www.nature.com/articles/ng.3470 (PMC4857881; PMID 26656844)
- SEgene — *npj Syst Biol Appl* 2025 — https://www.nature.com/articles/s41540-025-00533-x (PMC12089303; PMID 40389443)
- Kudo/Yamamoto et al. — *Mol Syst Biol* 2025 — https://pmc.ncbi.nlm.nih.gov/articles/PMC12130324/ (DOI 10.1038/s44320-025-00098-1)
- Watanabe et al. — 2024 — https://pmc.ncbi.nlm.nih.gov/articles/PMC11462929/ (PMID 38923173)
- Trends Genet review — 2022 — https://www.cell.com/trends/genetics/fulltext/S0168-9525(22)00148-2 (DOI 10.1016/j.tig.2022.11.001; PMID 35803787)
- Brief Bioinform survey (corroborating) — 2025 — bbaf596
