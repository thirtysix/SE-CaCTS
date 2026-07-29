"""cnrose.cn — copy-number providers + the separable correction transform (DESIGN.md §3, §4, §6)."""
from .base import CNTrack, CNProvider, ChainProvider, correct

__all__ = ["CNTrack", "CNProvider", "ChainProvider", "correct"]
