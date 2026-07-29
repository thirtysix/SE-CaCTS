# SE-CaCTS, explained simply

*A jargon-free overview — what this project is and why. Evergreen: no volatile numbers (those live in
`README.md` / `ROADMAP.md` / `DATA_SOURCES.md`).*

## The one-sentence version

We're building a way to find the DNA "control switches" that define each type of cancer — and to do it
honestly, so a switch only counts if it's genuinely special to that cancer, not just because the cancer made
extra copies of it.

## The pieces, in plain terms

- **Super-enhancers** are large stretches of DNA that act like master switches: they turn on the genes that
  make a cell what it *is*. Different cancers keep different switches "on."
- **How we see a switch is "on":** a chemical tag called **H3K27ac** sits on active DNA. Lots of public
  experiments have measured this tag across hundreds of cancer cell lines, so we can see which switches are
  active where.
- **The existing idea we're borrowing (CaCTS):** a published method that ranks *genes* by how *uniquely* active
  they are in one cancer type versus all others — the more one-cancer-specific, the more likely it's a "master"
  regulator. We apply that same specificity idea to the **switches (super-enhancers)** instead of the genes.
- **The honest twist (the new part):** cancers often duplicate chunks of their DNA (**amplification**). More
  copies of a switch → more tag signal → it *looks* more active and more "special," even if per-copy it's
  ordinary. We **correct for the number of DNA copies**, so a switch is only called cancer-specific if it's more
  active than its extra copies alone would explain. No one has published this copy-number-corrected version.
- **How we know the answers aren't noise:** any method like this needs a way to say "this result is stronger
  than chance." The honest way to check *that* is to scramble the labels — deliberately tell the method the
  wrong cancer type for every cell line, so there is genuinely nothing real left to find, and see how much it
  still "discovers." A trustworthy test finds essentially nothing. **Our first one didn't** — on scrambled
  labels it kept confidently reporting findings, which meant its confidence numbers were meaningless. We
  replaced it with one built from those scrambles directly, and re-checked: on scrambled data it now finds
  nothing, while the real biology still comes through. Every number we report comes from the replacement.

## Inputs → output

- **In:** the big public collection of switch-activity measurements across cancer cell lines; plus, for the same
  cell lines, their DNA copy numbers, gene activity, and which genes they can't survive without (DepMap).
- **Out:** for each cancer type, a ranked, trustworthy list of its defining switches — split into two honest
  buckets: switches that are *genuinely, specifically wired* for that cancer, versus ones that only look special
  because they're amplified.

## Why it matters

The defining switches of a cancer are promising drug targets and tell us how that cancer is wired. Doing it
across many cancer types at once, with the copy-number honesty built in, is what makes this both broadly useful
and genuinely new.

## How we're going about it (high level)

1. Check nobody's already done it (they haven't) and that the data is trustworthy enough (it is, with the right
   normalization) — both checked.
2. Assemble the data: match every switch-activity experiment to its cell line and that line's copy-number data.
3. Call the switches, measure their activity everywhere, make experiments from different labs comparable.
4. Correct for copy number, then score each switch's cancer-type specificity.
5. Validate against known biology, and connect the specific switches to their master genes.

Steps 1–5 are **done**. The method rediscovers the textbook switches for several cancer types without ever
being told what they were, and independently rediscovers well-known amplifications as amplifications rather
than as biology — which is the check that matters most, because that confusion is exactly what the method
exists to prevent.

## One honest limitation, stated up front

The results are trustworthy at the level of broad **cancer lineages** and **primary diseases**, but not at
finer subdivisions. The reason is simply how many cell lines we have per group: many fine-grained subtypes are
represented by a single cell line, and a group of one cannot be distinguished from a lucky draw of one. For
those finer levels the project shows a **ranked list** but deliberately refuses to attach a confidence claim
to it. Fixing this needs more cell lines per subtype — not more experiments on the lines we already have.

## Seeing it

The results are browsable as a website — every cancer type's specific switches, the copy-number comparison,
and the validation — at **<https://thirtysix.github.io/SE-CaCTS/>**. The code and data are open.
