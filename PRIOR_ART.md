# Prior art — is there already a "CaCTS for super-enhancers"?

Synthesis of a deep multi-source literature search (2026-07-16; 21 primary sources fetched, 25 claims verified
3–0 by adversarial voting, 0 refuted). Two of the most on-point papers (cSEAdb, Zhang 2016) were surfaced by the
search but under-weighted in its auto-synthesis; they were pulled from the run journal and are given full weight
here.

> **Step A update (2026-07-16, follow-up run) — novelty confirmed, SEdb 3.0 now fully checked.** A second
> deep-research pass (12 primary sources, 24 claims verified 3–0) closed the one project-killing open item: **SEdb
> 3.0 (NAR 2026) is fully read and does not preempt the method** — still ROSE-on-input-normalized H3K27ac, no
> specificity statistic, no CN correction. No 2024–26 work computes a copy-number-corrected per-SE
> lineage-specificity score. **Verdict holds.** Full evidence, citations, and the closest-analog table are in
> **`STEP_A_NOVELTY_CHECK.md`**; the landscape table and caveats below are updated to match.

## Bottom line

**Partly done — but not the version proposed here.** A per-super-enhancer specificity score from H3K27ac across a
compendium exists in a few flavors, and a cancer-lineage one exists; but a **copy-number-corrected** per-SE
specificity score — the one that stops an amplicon (MYC/CCAT1) from masquerading as lineage specificity — has
**not been published**. That CN-aware dimension is the genuinely novel core; breadth and TF/SE integration are
additional differentiators.

## Landscape, closest analog first

| method (year) | unit | specificity statistic | compendium | CN-aware? |
|---|---|---|---|---|
| **cSEAdb** — "epigenomic fingerprint of human cancers … at the constituent level" (PLOS Comput Biol 2024) | constituent enhancer within SEs | **presence/prevalence** (2-component mixture → active/inactive per constituent, flag if prevalence < ~0.21, threshold by minimizing variation-of-information) | **NCI-60: 60 lines / 28 cancer types** | **No** |
| **Ryu et al.** (BMC Bioinformatics 2019, 20(Suppl 3):127) | ROSE SE domain | **tau** (Yanai) on background-subtracted log2 H3K27ac RPM | 30 *normal/general* cell-tissue types | No |
| **SEA v4.0** (NAR 2025) | SE | **length-normalized Shannon entropy** (≈log2(n) = common; ≈0 = specific) | atlas-scale, mixed | No |
| **Zhang et al.** (Nature Genetics 2016, 48(2):176; Meyerson lab) | amplified SE | *none* — CN × H3K27ac **intersection** to find amplified driver SEs | TCGA CN × tissue H3K27ac, 12 tumor types | uses CN (opposite goal) |
| **SEdb 3.0** (NAR 2026) — *fully checked in Step A* | SE | ROSE calls only; new 3.0 tools = SE-BLAST + core-TF enrichment; **no compendium-wide specificity statistic** | **3,478,186 SEs / 5,387 H3K27ac samples**, 4 species | No (zero CN/amplicon terms in the paper) |
| **SEgene** (npj Syst Biol Appl 2025) | SE | peak-to-gene Pearson (H3K27ac CPM vs RNA-seq TPM, ±1 Mbp); ranks **commonness** (inverse of specificity) | multi-sample | No |
| **Su et al.** (BMC Genomics 2025) — the CN recipe | general ATAC / G4-ChIP peaks (**not SE, not H3K27ac**) | *none* — divide counts by CNVkit copy-number ratio ("signal per gene copy") | 2 paired disease/control samples | **Yes** — but never joined to SE specificity |
| **dbSUPER, dbCoRC** | SE | pairwise overlap / presence-absence only; **no compendium-wide statistic** | large; per-tissue | No |
| **CaCTS** (Reddy et al., Sci Adv 2021) — the template | TF (not SE) | **−log10 JSD** of expression vs one-hot ideal; empirical-null 5% FDR (1,000 label permutations) | TCGA pan-cancer RNA-seq (9,691 samples, 34 types / 140 subtypes) | No; **embraces** CN (MTFs coincide with CN gain, P=8.8e-4); SEs used only as downstream ROSE2 validation |

## The two papers that most shape the verdict

**cSEAdb (PLOS Comput Biol 2024)** — the closest *cancer* analog. Per-element lineage-specificity of SEs across
the NCI-60. But its metric is **presence/prevalence** from a two-component mixture model (active/inactive per
constituent), *not* a continuous divergence index — no JSD/tau/entropy/Gini — and it applies **no copy-number
correction**. So "score each SE's cancer-type specificity across a compendium" is precedented, but with a
categorical metric, a narrow 28-type panel, and no CN.

**Zhang et al. 2016 (Nat Genet)** — the closest *CN + SE* work, and instructive because it does the **inverse**
of SE-CaCTS. It combined somatic copy-number with lineage H3K27ac across 12 tumor types to *find focally amplified
super-enhancers* driving cancer genes — KLF5, USP12, PARD6B, and two MYC 3′ super-enhancers (MYC-LASE in lung
adeno, MYC-ECSE in endometrial). It **embraces** the CN↔H3K27ac coupling to nominate amplified drivers. SE-CaCTS
wants to **remove** that coupling so a locus isn't crowned "specific" merely because it is amplified. That the
OVCAR3 MYC/CCAT1 SE is both amplified and lineage-associated is exactly the confound Zhang exploited — a clean
motivation for why the correction is needed and non-trivial.

## Supporting facts (each verified 3–0)

- **CaCTS** uses `CaCTS_{i,j} = −log10 JSD(x_i, u_j)` on TCGA RNA-seq TF expression; SEs appear only as an
  orthogonal validation (ROSE2 across ~20 tumor types); no per-SE score; no CN in the score.
- **General vs individualized metrics matter.** tau/Shannon entropy (Ryu, SEA) are *general* (one number per SE =
  "how cell-restricted"), so they do **not** say *which* lineage. CaCTS is *individualized* (per gene × cancer
  type). The individualized/one-vs-rest form is what SE-CaCTS needs and is the less-trodden path on SEs.
- **The JSD machinery is off-the-shelf:** `tspex` (Camargo et al. 2020) implements Jensen-Shannon specificity
  (JSS), tau, Gini, entropy, etc., with both "general" and per-element "individualized" outputs — framed for
  expression but mechanically signal-agnostic.
- **The tissue-specificity benchmark** (Kryuchkova-Mostacci & Robinson-Rechavi 2017, Brief Bioinform 18(2):205)
  compared nine expression metrics and crowned **tau**; JSD was not even in the benchmark — so "tau" is the
  safest *defended* choice and "JSD on SEs" is genuinely new territory.
- **CN normalization of epigenomic signal is a solved, reusable recipe** (Su/Chen et al. 2025, BMC Genomics):
  estimate CN (CNVkit, from WGS or ChIP/ATAC input in 50 kb bins), then divide read counts by the copy-number
  ratio in gained regions (scale up losses). Established generally; **never joined to SE specificity**. Their
  worked ChIP case was G-quadruplex, not H3K27ac, and they never mention MYC — the MYC illustration is our
  accurate extension.
- **Adjacent atlases don't fill the gap:** EpiMap (Boix et al. 2021, Nature) builds tissue-specific enhancer
  *modules* but not SEs and no portable per-element specificity score; CRCmapper/coltron/dbCoRC map
  per-cell-type core regulatory circuitry with no cross-compendium specificity statistic.

## Novelty verdict

- "Specificity statistic on individual SEs from H3K27ac across a compendium" → **precedented** (Ryu tau; SEA
  entropy; cSEAdb prevalence).
- "The CaCTS JSD form ported to SEs" → **not found** (everyone uses tau/entropy/prevalence).
- "**Copy-number-corrected** SE specificity score" → **not published.** Zhang 2016 does the inverse; Su/Chen 2025
  gives the recipe but not for SEs.

⇒ Frame as **novel-in-combination**, genuinely **novel in the CN-aware dimension**, with breadth (~295 cancer
cell types vs cSEAdb's 28) and integration with the master-TF/CRC layer as the additional novelties. Do *not*
claim wholly unprecedented — cite Ryu / cSEAdb / SEA as the lineage.

## Caveats & open prior-art questions

- Several conclusions are **claims-from-absence** ("X does not correct for copy number") — each verified against
  full method text, but absence can't be proven exhaustively. "First / to our knowledge" is bounded by what the
  search harness could index (an unindexed preprint could exist; none surfaced).
- ~~**SEdb 3.0 (2026)** only partially read.~~ **Resolved (Step A):** SEdb 3.0 fully checked against the NAR 2026
  paper + PMC full text — still ROSE-on-input-normalized H3K27ac, no specificity/CN metric. See
  `STEP_A_NOVELTY_CHECK.md`.
- ~~Has anyone applied **JSD specifically** to per-SE H3K27ac across a cancer-lineage compendium?~~ **Answered: no.**
- ~~Do any amplicon-SE papers **implicitly** CN-correct when nominating lineage-specific SEs?~~ **Answered: no** —
  Kudo/MYB-SE 2025 explicitly declines; Watanabe/MYC-B-ALL 2024 embraces the coupling.
- **Still open (worth a final check before manuscript "first" language):** (a) the *live* SEdb 3.0 website could
  ship an undocumented specificity/CN feature not in the paper; (b) the **ecDNA / pediatric amplicon-SE** subfield
  (neuroblastoma MYCN, medulloblastoma) was not exhaustively covered.

## References (with links)

- Reddy et al., *Sci Adv* 2021 (CaCTS) — https://www.science.org/doi/10.1126/sciadv.abf6123
- Ryu et al., *BMC Bioinformatics* 2019 — https://pmc.ncbi.nlm.nih.gov/articles/PMC6439976
- cSEAdb, *PLOS Comput Biol* 2024 — https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011873
- SEA v4.0, *NAR* 2025 — https://pmc.ncbi.nlm.nih.gov/articles/PMC12807652/
- SEdb 2.0, *NAR* 2023 (D280) — https://academic.oup.com/nar/article/51/D1/D280/6786195
- dbSUPER (Khan & Zhang), *NAR* 2016 — https://pmc.ncbi.nlm.nih.gov/articles/PMC4702767
- Zhang et al., *Nature Genetics* 2016 — https://www.nature.com/articles/ng.3470
- Su/Chen et al., *BMC Genomics* 2025 (CN normalization) — https://link.springer.com/article/10.1186/s12864-025-11442-y
- tspex (Camargo et al. 2020) — https://apcamargo.github.io/tspex/
- Kryuchkova-Mostacci & Robinson-Rechavi, *Brief Bioinform* 2017 — https://pmc.ncbi.nlm.nih.gov/articles/PMC5444245
- EpiMap (Boix et al.), *Nature* 2021 — https://pubmed.ncbi.nlm.nih.gov/33536621/

*Added in Step A (2026-07-16 follow-up; see `STEP_A_NOVELTY_CHECK.md`):*
- SEdb 3.0, *Nucleic Acids Res* 2026, 54(D1):D322 — https://academic.oup.com/nar/article/54/D1/D322/8373947 (PMC12807715)
- SEgene, *npj Syst Biol Appl* 2025 — https://www.nature.com/articles/s41540-025-00533-x (PMC12089303)
- Kudo/Yamamoto et al. (MYB-SE), *Mol Syst Biol* 2025 — https://pmc.ncbi.nlm.nih.gov/articles/PMC12130324/
- Watanabe et al. (MYC BENC-CNA, B-ALL) 2024 — https://pmc.ncbi.nlm.nih.gov/articles/PMC11462929/ (PMID 38923173)
- Trends in Genetics SE-subtype review 2022 — https://www.cell.com/trends/genetics/fulltext/S0168-9525(22)00148-2 (PMID 35803787)
