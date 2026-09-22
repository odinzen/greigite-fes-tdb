#!/usr/bin/env python3
"""Regenerate every manuscript figure (Fig. 1-4, S1-S7) in one command.

Runs each manuscript/make_Figure_*.py plus the porewater (Fig. 3) and kinetics (Fig. 4)
scripts with the current Python interpreter; each writes its PNG(s) to artifacts/figures/.
Build the databases and run engine/validate_fes_engine.py first (see the README "Quick
start"): the Fe-S figures read artifacts/tdb/ and artifacts/fes_engine_boundaries.json.
Fig. 3 needs phreeqpython; its PHREEQC grid (artifacts/aqueous/ehph_fields.npz) is computed
on first run. Fig. S8 (powder XRD) was produced outside this repository.

    python manuscript/make_all_figures.py
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent

scripts = sorted(HERE.glob("make_Figure_*.py"))
if not scripts:
    sys.exit("no make_Figure_*.py scripts found next to this driver")
if not (REPO / "artifacts" / "aqueous" / "ehph_fields.npz").exists():
    scripts.append(REPO / "aqueous" / "build_ehph_diagram.py")
scripts += [REPO / "aqueous" / "make_ehph_figure.py", REPO / "kinetics" / "make_kinetics_figure.py"]

print(f"Regenerating {len(scripts)} figure scripts -> artifacts/figures/\n")
failed = []
for s in scripts:
    print(f"  {s.relative_to(REPO)} ...", flush=True)
    r = subprocess.run([sys.executable, str(s)], capture_output=True, text=True)
    if r.returncode != 0:
        failed.append(s.name)
        print(f"    FAILED (exit {r.returncode}):\n{r.stderr.strip()[-2000:]}\n")

if failed:
    sys.exit(f"\n{len(failed)} script(s) failed: {', '.join(failed)}")
print(f"\nAll {len(scripts)} figure scripts succeeded.")
