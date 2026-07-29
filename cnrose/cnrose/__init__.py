"""cnrose — a bigWig-native, external-CN, ROSE-style super-enhancer caller.

Agnostic core (this milestone): io.quantify + stitch + callsuper, ported faithfully from ROSE2
(see ../DESIGN.md §7 for the validation contract). CN provider layer lands next (DESIGN.md §8 step 2).
"""
__version__ = "0.1.0-dev"

from .stitch import stitch
from .callsuper import calculate_cutoff, call_super
from .io import quantify

__all__ = ["stitch", "calculate_cutoff", "call_super", "quantify"]
