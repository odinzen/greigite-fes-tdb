#!/usr/bin/env python3
"""Regenerate every published manuscript figure in one command.

Runs each manuscript/make_Figure_*.py (Fig 1-6, S1-S6) with the current Python
interpreter; each script writes its PNG(s) to artifacts/figures/. Build the
thermodynamic databases first (see the README "Quick start"), since the figure
scripts read them from artifacts/tdb/.

    python manuscript/make_all_figures.py
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

scripts = sorted(HERE.glob("make_Figure_*.py"))
if not scripts:
    sys.exit("no make_Figure_*.py scripts found next to this driver")

print(f"Regenerating {len(scripts)} manuscript figure scripts -> artifacts/figures/\n")
failed = []
for s in scripts:
    print(f"  {s.name} ...", flush=True)
    r = subprocess.run([sys.executable, str(s)], capture_output=True, text=True)
    if r.returncode != 0:
        failed.append(s.name)
        print(f"    FAILED (exit {r.returncode}):\n{r.stderr.strip()[-2000:]}\n")

if failed:
    sys.exit(f"\n{len(failed)} script(s) failed: {', '.join(failed)}")
print(f"\nAll {len(scripts)} figure scripts succeeded.")
