"""Porewater Eh-pH stability diagram for the Fe-S(-C-O-H) system with stable greigite.

Predominance method: over a pH x pe grid, at fixed low total dissolved Fe and diagenetic total S
and C, compute the saturation index (SI) of every candidate Fe mineral at the imposed (pH, pe) and
take the most-supersaturated solid (SI per Fe atom) as the predominant phase. Mackinawite is
excluded (metastable precursor). Greigite log_K is the CODATA-referenced value from
derive_greigite_logk.py; the grid is run for nominal / +1sigma / -1sigma to bound the greigite
field.

Needs phreeqpython (pip install phreeqpython) + numpy. The reference database is phreeqpython's
bundled phreeqc.dat; the Fe-S(-O) phases it lacks are grafted on in memory. Writes the field grids
to artifacts/aqueous/ehph_fields.npz (gitignored); make_ehph_figure.py renders it.
"""
import re, os, pathlib, numpy as np
import phreeqpython
from phreeqpython import PhreeqPython

REPO = pathlib.Path(__file__).resolve().parents[1]
ART = REPO / "artifacts" / "aqueous"
ART.mkdir(parents=True, exist_ok=True)
DB_SRC = pathlib.Path(os.path.dirname(phreeqpython.__file__)) / "database" / "phreeqc.dat"

GREIGITE_LOGK = {"nominal": -68.95, "hi": -65.12, "lo": -72.79}  # from LOGK_DERIVATION.md
PYRRHOTITE_LOGK = -6.03

# candidate solids for the predominance map (Mackinawite excluded = metastable precursor)
CANDIDATES = ["Greigite", "Pyrite", "Pyrrhotite", "Magnetite", "Hematite",
              "Goethite", "Siderite", "Fe(OH)3(a)"]
# Fe atoms per dissolution reaction: predominance compares SI/n_Fe (lowest imposed a_Fe wins),
# NOT raw SI, which is biased toward multi-Fe phases (hematite/magnetite).
N_FE = {"Greigite": 3, "Pyrite": 1, "Pyrrhotite": 1, "Magnetite": 3,
        "Hematite": 2, "Goethite": 1, "Siderite": 1, "Fe(OH)3(a)": 1}


def _analytic_logk(text, T=298.15):
    import math
    nums = re.findall(r"-?\d+\.?\d*[eE]?[-+]?\d+|-?\d+\.?\d*", text)
    A = [float(c) for c in nums[:5]] + [0.0] * (5 - len(nums[:5]))
    return A[0] + A[1] * T + A[2] / T + A[3] * math.log10(T) + A[4] / T ** 2


def _inject_logk(txt):
    """Some IPhreeqc builds reject a PHASES entry that carries only `-analytic` and no explicit
    `log_k`. Buffer each entry; if it has a reaction and an -analytic line but no log_k, insert
    log_k(25C) from the analytic expression so the database loads."""
    def is_header(l):
        return bool(l) and not l[0].isspace() and l.strip() and "=" not in l and not l.strip().startswith("-")
    out, entry = [], []

    def flush():
        if entry:
            has_logk = any(re.search(r"log_k", x, re.I) for x in entry)
            has_rxn = any("=" in x for x in entry)
            ana = next((x for x in entry if re.match(r"\s*-analytic", x)), None)
            if has_rxn and not has_logk and ana:
                for x in entry:
                    out.append(x)
                    if "=" in x and not x.strip().startswith("-"):
                        out.append(f"        log_k   {_analytic_logk(ana):.4f}")
            else:
                out.extend(entry)
    for l in txt.splitlines():
        if is_header(l):
            flush(); entry = [l]
        else:
            entry.append(l)
    flush()
    return "\n".join(out)


def make_db(greigite_logk):
    """phreeqc.dat (loads clean; pyrite log_k -18.479 = the reference we validated against) plus
    the Fe-S(-O) phases it lacks: greigite (this work), pyrrhotite/troilite, magnetite. The new
    phases are inserted INSIDE the existing PHASES block (content after the file's trailing END is
    ignored by the loader), right after the `PHASES` header."""
    txt = DB_SRC.read_text(errors="ignore")
    inject = f"""Greigite
        Fe3S4 + 4H+ = 2Fe+3 + Fe+2 + 4HS-
        log_k   {greigite_logk}
Pyrrhotite
        FeS + H+ = Fe+2 + HS-
        log_k   {PYRRHOTITE_LOGK}
Magnetite
        Fe3O4 + 8H+ = 2Fe+3 + Fe+2 + 4H2O
        log_k   3.737
"""
    m = re.search(r"^PHASES\s*\n", txt, re.M)
    if not m:
        raise RuntimeError("no PHASES section found in database")
    txt = txt[:m.end()] + inject + txt[m.end():]
    out = ART / "_db_tmp.dat"
    out.write_text(txt)
    return out


def predominant(pp, pH, pe, feT=1e-6, sT=1e-3, cT=1e-3):
    """Return (name, SI) of the most-supersaturated candidate solid, or ('aqueous', si)."""
    try:
        # total S with NO fixed valence so the imposed pe sets the sulfate/sulfide split
        sol = pp.add_solution({"units": "mol/kgw", "pH": pH, "pe": pe,
                               "Fe": feT, "S": sT, "C(4)": cT, "Na": 2 * sT, "Cl": 1e-9})
    except Exception:
        return ("error", np.nan)
    best, best_norm, best_si = "aqueous", -99.0, np.nan
    for m in CANDIDATES:
        try:
            si = sol.si(m)
        except Exception:
            continue
        if si is None:
            continue
        norm = si / N_FE[m]          # SI per Fe atom -> lowest imposed a_Fe wins
        if norm > best_norm:
            best, best_norm, best_si = m, norm, si
    pp.remove_solutions([sol.number])
    # a solid predominates only where it is at/above saturation (SI/n_Fe >= 0)
    return (best, best_si) if best_norm >= 0 else ("aqueous", best_norm)


if __name__ == "__main__":
    pH_grid = np.arange(2.0, 12.001, 0.1)
    pe_grid = np.arange(-12.0, 17.001, 0.2)
    labels = {name: i for i, name in enumerate(["aqueous"] + CANDIDATES + ["error"])}

    fields = {}
    for variant, glk in GREIGITE_LOGK.items():
        db = make_db(glk)
        pp = PhreeqPython(database=str(db))
        grid = np.zeros((len(pe_grid), len(pH_grid)), dtype=int)
        for j, pe in enumerate(pe_grid):
            for i, pH in enumerate(pH_grid):
                name, _ = predominant(pp, pH, pe)
                grid[j, i] = labels[name]
        fields[variant] = grid
        del pp
        print(f"[{variant:7}] greigite log_K={glk}: field counts:",
              {n: int((grid == labels[n]).sum()) for n in ["Greigite", "Pyrite", "Pyrrhotite",
               "Hematite", "Goethite", "Magnetite", "Siderite", "aqueous"]})

    # ---- sanity checks against known Fe-S geochemistry ----
    print("\nSanity checks (nominal):")
    pp = PhreeqPython(database=str(make_db(GREIGITE_LOGK["nominal"])))
    for (pH, pe, expect) in [(7.0, -3.0, "sulfidic reducing -> pyrite/greigite/pyrrhotite"),
                             (7.0, 8.0, "oxidizing -> Fe(III) oxide"),
                             (4.0, 10.0, "acidic oxidizing -> aqueous/oxide"),
                             (8.0, -6.0, "very reducing sulfidic -> pyrrhotite/greigite")]:
        name, si = predominant(pp, pH, pe)
        print(f"  pH={pH:4} pe={pe:5} -> {name:11} (SI={si:6.2f})   [{expect}]")

    np.savez(ART / "ehph_fields.npz", pH=pH_grid, pe=pe_grid,
             labels=np.array(list(labels.keys())), **fields)
    print("\nsaved", ART / "ehph_fields.npz")
