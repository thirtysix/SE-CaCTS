# Running the Step-B pilot on CSC Roihu (`$SECACTS_CSC_PROJECT`)

Roihu supersedes Puhti (compute off **31 Jul 2026**) and Mahti (**31 Aug 2026**). CPU-only pipeline →
ports cleanly (no ARM/GPU recompile concerns here). Quotas are tighter: **15 GiB** projappl, **250 GiB**
scratch — the pilot needs a few GB, fine.

## One-time setup (yours to do)

1. **Add Roihu as a service** to `$SECACTS_CSC_PROJECT` in MyCSC (existing billing units carry over; no new
   project needed).
2. **SSH alias** — add to `~/.ssh/config` (you already have `puhti`/`mahti`). ⚠️ Roihu's login nodes are
   **`roihu-cpu.csc.fi`** (86.50.172.16) and **`roihu-gpu.csc.fi`** (86.50.172.21) — **NOT** `roihu.csc.fi`
   (that round-robin timed out from our network). Force IPv4 (laptop has no IPv6 route) + offer only the CSC key:
   ```
   Host roihu roihu-cpu.csc.fi
       HostName roihu-cpu.csc.fi
       User barker
       IdentityFile ~/.ssh/id_ed25519
       IdentitiesOnly yes
       AddressFamily inet
       StrictHostKeyChecking accept-new
   # plus a matching `Host roihu-gpu` block -> roihu-gpu.csc.fi (GPU nodes are ARM GH200)
   ```
   Verify (prefix with `!`):  `! ssh roihu 'hostname; ls -d /scratch/$SECACTS_CSC_PROJECT'`
   ✅ **The real gotcha — Roihu uses SSH _certificates_, not raw keys.** The Puhti key is rejected
   (`Permission denied (publickey)`) until you sign it. **Every 24 h**, to open a *new* connection:
   - MyCSC → Profile → SSH public keys → ⋮ → **Sign and download SSH certificate**  *(or CLI
     `python3 ~/.ssh/csc_cert.py -u barker ~/.ssh/id_ed25519.pub`)*, then save it as **`~/.ssh/id_ed25519-cert.pub`**
     (OpenSSH auto-loads it next to the key — no `CertificateFile` line needed).
   - Also **accept the Roihu ToS once** (MyCSC → Projects → Services → Roihu). Check validity:
     `ssh-keygen -L -f ~/.ssh/id_ed25519-cert.pub | grep Valid`.
   - **Submitted jobs keep running after the cert expires**; you only re-sign to log back in (SSH/rsync/scp).
     An already-open connection also survives expiry (cert checked only at connect).
3. **Toolchain (checked 2026-07-18 — access works):** **deepTools is NOT a Roihu module → containerize with
   tykky** (`00_build_env.sh` → `/projappl/$SECACTS_CSC_PROJECT/secacts_tykky`; `pilot.slurm` already falls back to
   it). `/scratch/$SECACTS_CSC_PROJECT` + `/projappl/$SECACTS_CSC_PROJECT` exist. **CPU partitions:** `small` (default,
   3-day — used by `pilot.slurm`), `medium`/`large` (1.5-day), **`longrun` (10-day — for the big pull)**,
   `interactive`, `test`, `hugemem`. Login node `roihu-cpu-login1`.

## Run (from the laptop, once the alias works)

```bash
conda activate atac_hdac
cd pilot/scripts
bash 03_download_bigwigs.sh --go        # download the ~13 bigWigs LOCALLY (login-node discipline: don't pull multi-GB on the login node)
cd ../roihu
bash stage_and_submit.sh                # rsync scripts+metadata+bigWigs to scratch, build env, sbatch
# monitor / pull results as printed by the script
```

## Layout on Roihu

```
/scratch/$SECACTS_CSC_PROJECT/se-cacts/pilot/   scripts/ data/{selection.tsv,bigwigs/} results/
/projappl/$SECACTS_CSC_PROJECT/secacts_tykky/   deepTools + numpy/pandas/scipy/matplotlib (if no module)
```

## Notes / discipline

- **Login node = only internet.** For the *pilot* we download the ~13 bigWigs on the laptop and rsync up
  (respects the "no multi-GB downloads on the login node" rule). For the **Phase-1** full pull (~1,789
  experiments, 100s of GB) neither laptop-download nor a naive login-node loop scales — that needs a
  proper staged/chunked transfer plan (a separate design step; ChIP-Atlas also offers bulk endpoints).
- **Reality check:** the pilot compute (multiBigwigSummary over 13 bigWigs + a numpy PCA) is minutes and a
  few GB RAM — it runs fine locally. Roihu here is mainly a **rehearsal** of the Roihu flow before the
  heavy Phase-1 work. The big HPC payoff is Phase 1–4, not this pilot.
- `pilot.slurm` sources `zz-csc-env.sh` *before* `set -u` (the CSC env script isn't `-u`-clean) and runs
  `python -u` so progress streams to the `.out` live.
