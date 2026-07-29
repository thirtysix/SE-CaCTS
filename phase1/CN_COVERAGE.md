# CN coverage — the lines we can investigate, by copy-number provenance

**What this answers:** for every cell line in scope, do we have copy-number (CN) data, and from
where — DepMap WGS (already in the atlas), another measured database (a lookup), signal-based
**inference**, or nothing? CN matters because the SE caller is CN-aware: scoring-time correction
(`cnrose.cn.correct`, DESIGN §3.2) needs a per-line CN profile, so CN coverage bounds how much of the
panel can be corrected rather than run agnostic-only.

**Provenance / dates.** Everything here is reproducible from local files. The funnel and per-line
manifest counts come from `phase1/data/phase1_manifest.tsv`; the four-tier CN split is produced by
**`phase1/scripts/12_cn_source_audit.py` → `phase1/data/cn_source.tsv`** (one row per line, with the tier
and the underlying DR flags). It reproduces the **method** of the one-off 2026-07-18 audit
(`cnrose/DESIGN.md §4`), but now runs on the **identity-CORRECTED** manifest (see
[`CELL_LINE_IDENTITY.md`](CELL_LINE_IDENTITY.md) — a sentinel-sink fix and RRID-canonicalization in
`scripts/11`). §4's numbers were 254/119/25 on the buggy pre-fix identities; the corrected split is
**260/113/25** — several homonym/subline CVCLs were fixed to their canonical human accessions, whose
Cellosaurus `DR` blocks carry more measured-CN cross-refs (so lines moved C→B). **This supersedes §4.**
Line counts are distinct CVCL unless noted.

---

## 1. The funnel (ChIP-Atlas H3K27ac → analysis-ready)

| Stage | Experiments | Lines (CVCL) |
|---|---:|---:|
| All H3K27ac hg38 in ChIP-Atlas | 11,827 | 740 |
| → resolved to a CVCL identity (sink-guarded) | 5,373 | — |
| QC-pass | 10,775 | 722 |
| QC-pass **+ DepMap-joined** | 3,583 | 463 |
| QC-pass **+ DepMap WGS CN** ← **the pull (built)** | **2,917** | **324** |

The atlas as built is **2,916** usable experiments — SRX20868733 (retired, 404) was dropped. (Before the
identity fix the CVCL-resolved count read 9,579 and the line count 721; ~4,200 of those were a single
"NA"-sentinel sink and are now correctly unmapped — the DepMap-anchored rows were unaffected.)

---

## 2. The 722 lines, by CN provenance

"All the lines we would investigate" = the full **722** QC-pass CVCL lines (current pull + expansion).
Every line is in exactly one tier; the four tiers partition the 722 exactly.

| Tier | Lines | % of 722 | CN status | Phase-3 action |
|---|---:|---:|---|---|
| **A — DepMap WGS** (current pull) | 324 | 45% | measured, **in the atlas now** | none — already corrected-capable |
| **B — other measured DB** | 260 | 36% | measured elsewhere (DepMap-merged WES/SNP, Progenetix, CMP, COSMIC) | **lookup** — new `CNProvider` backend, no inference |
| **C — input-inferred** | 113 | 16% | none direct; has a matched input | **inference** from signal (the actual modeling task) |
| **D — orphan** | 25 | 3% | none | agnostic catalog only (`cn_source=none`) |
| **Total** | **722** | 100% | | |

- **CN available (A+B+C) = 697 / 722 = 96.5%.** Only **25 lines (3.5%)** are true orphans.
- Expansion arithmetic: current pull **324** + net-new **398** = **722**; net-new **398** = B(260) + C(113) + D(25).
- Experiments: pull **2,917** + net-new **~5,900** ≈ **8,800**.

**Per-resource reach among the 398 net-new lines** (lines are multiply-covered, so these overlap and
do **not** sum to 260): Progenetix 251 · CMP 159 · DepMap-merged WES/SNP 146 · COSMIC 97 · GDSC 100.
Resolution order (cheapest/best first, `DESIGN.md §4`): **DepMapMerged** (one download of
`OmicsCNGene.csv`) → **Progenetix** (CVCL-native, biggest reach) → **CellModelPassports**
(uniquely *absolute integer* CN) → **CosmicCLP** (licence friction; top-up only) → **InputInferred**
(the 113) → orphan. (Reaches rose vs §4's 245/153/140/95/98 because the identity fix corrected homonym
CVCLs to human accessions with richer DR blocks.)

---

## 3. What actually needs *inference* (the number that's easy to overstate)

The headline "~94% of net-new lines have CN" (= (260+113)/398) is a **coverage** figure, not the
inference workload. Of the whole 722:

- **324** already have DepMap WGS CN (in the atlas) — nothing to do;
- **260** are a **database lookup** — a new backend, no signal modeling;
- **113** (16% of 722) require **genuine signal-based CN inference** — this is the hard part;
- **25** stay CN-agnostic.

So the real inference target is **~113 lines, not ~680.** Wording elsewhere ("CN inference for
the ~94% of lines lacking DepMap CN") conflated coverage with workload — corrected here.

### Two candidate inference substrates (open design decision)

The 113 can get CN from two different data sources, and the docs don't yet commit:

1. **Matched-input inference** (`DESIGN.md §4`): CN from each line's **matched Input-control bigWig**
   via HMCan / Control-FREEC / CNVkit. Needs a real input, and needs **gating on input read depth**
   ("has an input" ≠ "has a *good* input"; coarse on focal amplicons).
2. **H3K27ac-background inference** (`ROADMAP.md` Phase 3): LILY-style CN from the **H3K27ac signal
   itself** (the 1 kb genome archive), no matched input required — the durable purpose of that archive.

Either way, **validate inferred CN against the 324 Tier-A lines that have DepMap WGS CN** before
trusting it. This choice (and the validation) is the Phase-3 design task that is still open.

---

## 4. Reproducible detail from the manifest (Phase-3 prep)

Computed here from `phase1_manifest.tsv`; unlike the tier split, these are re-derivable any time.

- **Read depth** (H3K27ac experiments; a *proxy* — the input-depth gate in Tier C needs the
  *input-control* depth, which is not in this manifest): pull median **32.6 M** (p10 14.8 M, p90 63.6 M);
  net-new median **30.4 M** (p10 15.8 M, p90 60.4 M); floor **5.0 M** (QC min). Net-new is not
  systematically shallower than the pull.
- **Lineage gap — now RESOLVED (was 259 lineage-less).** The 259 net-new lines lacked an Oncotree lineage
  only because the join never read Cellosaurus' own disease field, not because the data was missing.
  `scripts/14_cellosaurus_lineage.py` recovers it from Cellosaurus `DI` (NCIt disease) via an empirical
  NCIt→Oncotree crosswalk, and gates the set to human cancer lines. Net effect on the full scope:

  | | before | after |
  |---|---:|---:|
  | lines with a lineage | 462 | **631** (460 DepMap + 171 recovered) |
  | lineage-less | 259 | 0 (every *included* line has one) |
  | excluded (non-cancer / non-human / DepMap-Non-Cancerous) | — | 91 |

  The input-inferred tier that looked "95% lineage-less" was mostly **non-cancer immortalised-normal**
  lines (HEK293, HaCaT, BEAS-2B, …); once those are excluded per scope, the included bucket C is 53 and all
  have a lineage. Full per-line detail: [`CELL_LINE_IDENTITY.md`](CELL_LINE_IDENTITY.md) +
  `data/lineage_resolved.tsv`.

---

## How the tiers are decided (the audit predicate)

`12_cn_source_audit.py` classifies each QC-pass CVCL, all from local files, no network:

- **Tier A** — the line is in the pull (`phase1_manifest.tsv` `has_cn==1`: its DepMap `ModelID` is in
  `OmicsCNGeneWGS.csv`).
- **Tier B (measured_other)** — not A, but the Cellosaurus entry carries a `DR` cross-ref to a
  measured-CN resource: `DepMap` (→ merged `OmicsCNGene.csv`, WES/SNP — recovers non-WGS DepMap lines),
  `cancercelllines` (Progenetix), `Cell_Model_Passport` (CMP), or `Cosmic-CLP` (COSMIC Cell Line
  Project). **Note:** the predicate uses `Cosmic-CLP`, *not* the 18k generic `Cosmic;` mutation-sample
  rows (nearly every line has those). GDSC is reported for reach only — its id == the Cosmic-CLP id, so
  it adds ~no lines beyond COSMIC.
- **Tier C (input_inferred)** — not A/B, but the line has a matched hg38 **Input-control** experiment in
  `experimentList.tab` (matched by CVCL via the same Cellosaurus name map, or by raw cell string).
- **Tier D (none)** — none of the above; agnostic catalog only.

Two soft spots to know: a `DR` cross-ref means the *resource lists the line*, taken as a proxy for "CN
obtainable there" (true for the CN databases Progenetix/CMP; for DepMap it means the merged file — which
we don't hold locally — would carry it). And matched-input is **line-level** ("has an input somewhere"),
not same-study/same-sample — the stricter matched-input quality gate (input read depth) is a Phase-3
step, not done here.

## Reproduce

```bash
cd phase1/scripts
python3 12_cn_source_audit.py     # -> ../data/cn_source.tsv + the reconciliation report vs §4 on stderr
```

The report prints the 722-line partition, the net-new measured/input/orphan split, and per-resource
reach, each next to the pre-fix §4 figure (the shift is expected — see the provenance note at the top).
`cn_source.tsv` (checked in) is the per-line authority; `cnrose/DESIGN.md §4` is the design write-up of
the same tiers on the pre-fix identities.
