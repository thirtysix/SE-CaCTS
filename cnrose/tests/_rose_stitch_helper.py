"""Ground-truth stitch oracle — runs under the se_rose env (which has the rose2 package).

Reads a 3-col constituents BED and prints ROSE2's own stitched regions (chrom<TAB>start<TAB>end),
via rose2.utils.LocusCollection.stitchCollection(window, 'both'). This is the reference cnrose.stitch
must reproduce. Invoked by ../tests/validate_vs_rose2.py with the se_rose python.

    ~/miniconda3/envs/se_rose/bin/python _rose_stitch_helper.py constituents.bed 12500
"""
import os
import sys

# rose2/__init__.py prints a banner to stdout on import — suppress it so only region lines reach stdout.
_saved = os.dup(1)
_null = os.open(os.devnull, os.O_WRONLY)
os.dup2(_null, 1)
try:
    from rose2 import utils
finally:
    sys.stdout.flush()
    os.dup2(_saved, 1)
    os.close(_null)
    os.close(_saved)

bed, window = sys.argv[1], int(sys.argv[2])
loci = []
with open(bed) as fh:
    for line in fh:
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.rstrip("\n").split("\t")
        loci.append(utils.Locus(f[0], int(f[1]), int(f[2]), ".", f"{f[0]}:{f[1]}-{f[2]}"))

lc = utils.LocusCollection(loci, 50)
stitched = lc.stitchCollection(window, "both")
regs = sorted((lo.chr(), lo.start(), lo.end()) for lo in stitched.getLoci())
for c, s, e in regs:
    sys.stdout.write(f"{c}\t{s}\t{e}\n")
