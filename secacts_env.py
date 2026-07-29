#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Machine-local configuration for SE-CaCTS, read from `.env` (gitignored) or the environment.

Nothing in the tracked tree may hardcode a local filesystem path or an HPC allocation id;
they live in `.env` only. `cp sample.env .env` and edit. See `sample.env` for every key.

Usage (the repo root must be on sys.path — every caller already computes it from __file__):

    sys.path.insert(0, SECACTS)
    from secacts_env import DATAROOT, CSC_PROJECT, CACHE_DIR

`DATAROOT` is resolved lazily and raises a clear, actionable error if unset, so importing
this module never fails on a machine that only needs, say, `CACHE_DIR`.

Real environment variables take precedence over `.env`, so a one-off override works:

    SECACTS_DATAROOT=/data/mirror python phase2/score_pilot.py ...

Deliberately dependency-free (no python-dotenv): this is imported by scripts that run on a
bare cluster venv carrying only pyBigWig/numpy/scipy/pandas.
"""
from __future__ import annotations

import os

ROOT = os.path.dirname(os.path.abspath(__file__))
#: Overridable for tests and for staged cluster copies, matching secacts_env.sh.
ENV_FILE = os.environ.get("SECACTS_ENV_FILE") or os.path.join(ROOT, ".env")


def _parse_env_file(path: str) -> dict:
    """Minimal KEY=VALUE reader: ignores blanks/#comments, strips one layer of quotes,
    tolerates `export KEY=...`, and expands $VARS against the current environment."""
    out: dict = {}
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            if line.startswith("export "):
                line = line[len("export "):].lstrip()
            key, _, val = line.partition("=")
            key, val = key.strip(), val.strip()
            if len(val) >= 2 and val[0] == val[-1] and val[0] in "\"'":
                val = val[1:-1]
            out[key] = os.path.expandvars(os.path.expanduser(val))
    return out


_FILE_ENV = _parse_env_file(ENV_FILE)


def get(key: str, default=None):
    """Real environment first, then `.env`, then the default."""
    v = os.environ.get(key)
    if v not in (None, ""):
        return v
    v = _FILE_ENV.get(key)
    return v if v not in (None, "") else default


def require(key: str, why: str = "") -> str:
    v = get(key)
    if not v:
        msg = [f"[secacts_env] {key} is not set."]
        if why:
            msg.append(f"  {why}")
        msg += [
            f"  Fix: cp {os.path.join(ROOT, 'sample.env')} {ENV_FILE} and set {key},",
            f"       or export {key}=... for this run.",
        ]
        raise SystemExit("\n".join(msg))
    return v


_DATAROOT_WHY = (
    "It is the parent directory holding DepMap/, chip-atlas/, cellosaurus/, 0.human_genome/ "
    "and the sibling pyCaCTS checkout (see sample.env for the exact layout)."
)


def __getattr__(name: str):
    """PEP 562 module-level lazy attributes.

    `from secacts_env import DATAROOT` resolves here and returns a plain `str`, raising a
    clear SystemExit only for callers that actually ask for it — so a script needing only
    CACHE_DIR still imports cleanly on a machine with no `.env`. Using a module __getattr__
    rather than a lazy str subclass keeps DATAROOT an ordinary string, so f-strings,
    os.path.join and pathlib all behave exactly as they look.
    """
    if name == "DATAROOT":
        return require("SECACTS_DATAROOT", _DATAROOT_WHY)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


#: SLURM allocation for the Phase-2 HPC pull. Optional — everything downstream runs locally.
CSC_PROJECT = get("SECACTS_CSC_PROJECT")

#: Parsed-once caches (gene coordinates, ...). Repo-local and gitignored by default.
CACHE_DIR = get("SECACTS_CACHE_DIR") or os.path.join(ROOT, ".cache")


def cache_path(name: str) -> str:
    """Absolute path for a named cache artifact, creating the directory on demand."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, name)


def pycacts_path() -> str:
    """The sibling pyCaCTS checkout — the JSD scoring engine this project reuses.

    Calls require() rather than referencing DATAROOT: a module-level __getattr__ (PEP 562)
    only fires for *external* attribute access, so a bare `DATAROOT` inside this module
    would raise NameError.
    """
    return os.path.join(require("SECACTS_DATAROOT", _DATAROOT_WHY), "002.AI_projects", "pyCaCTS")


if __name__ == "__main__":                              # `python secacts_env.py` = show config
    print(f"repo root       {ROOT}")
    print(f".env            {ENV_FILE} ({'found' if os.path.exists(ENV_FILE) else 'MISSING'})")
    print(f"DATAROOT        {get('SECACTS_DATAROOT') or '** UNSET — cp sample.env .env **'}")
    print(f"CSC_PROJECT     {CSC_PROJECT or '(unset; only needed for the HPC pull)'}")
    print(f"CACHE_DIR       {CACHE_DIR}")
