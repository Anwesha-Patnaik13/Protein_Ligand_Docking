import os
import csv
from rdkit import Chem
from rdkit.Chem import rdMolAlign, AllChem

THRESHOLD = 2.0
N_POSES = 20

# ---- Check 1: required files exist ----
required = [
    "protein/receptor.pdbqt",
    "ligand/ligand.pdbqt",
    "ligand/crystal_A.pdb",
    "ligand/imatinib_original.sdf",
    "results/docked_v2.pdbqt",
] + [f"results/v2_pose{i}.pdb" for i in range(1, N_POSES + 1)]

missing = [f for f in required if not os.path.exists(f)]
print("CHECK 1 - Required files:", "PASS" if not missing else f"FAIL, missing {missing}")
if missing:
    raise SystemExit(1)

# ---- Reference and template ----
ref_raw = Chem.MolFromPDBFile("ligand/crystal_A.pdb", removeHs=True)
suppl = Chem.SDMolSupplier("ligand/imatinib_original.sdf", removeHs=True)
template = Chem.GetMolFrags(suppl[0], asMols=True)[0]
ref = AllChem.AssignBondOrdersFromTemplate(template, ref_raw)

# ---- Affinities from Vina output ----
affinities = [float(l.split()[3]) for l in open("results/docked_v2.pdbqt")
              if l.startswith("REMARK VINA RESULT")]
print("CHECK 2 - Number of poses matches:",
      "PASS" if len(affinities) == N_POSES else f"FAIL ({len(affinities)} found)")

# ---- RMSD for every pose ----
rows = []
for i, aff in enumerate(affinities, start=1):
    raw = Chem.MolFromPDBFile(f"results/v2_pose{i}.pdb", removeHs=True)
    pose = AllChem.AssignBondOrdersFromTemplate(template, raw)
    rmsd = rdMolAlign.CalcRMS(pose, ref)
    rows.append((i, aff, round(rmsd, 3), "YES" if rmsd < THRESHOLD else "NO"))

# ---- Check 3: values are sensible ----
ok = all(0 <= r[2] < 50 for r in rows)
print("CHECK 3 - RMSD values are finite and reasonable:", "PASS" if ok else "FAIL")

# ---- Save table ----
with open("results/validation_results.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Pose", "Affinity_kcal_mol", "RMSD_A", "Within_2A"])
    w.writerows(rows)

# ---- Summary ----
top = rows[0]
best = min(rows, key=lambda r: r[2])
n_good = sum(1 for r in rows if r[3] == "YES")
print()
print(f"Top-ranked pose: pose {top[0]}, {top[1]} kcal/mol, RMSD {top[2]} A")
print(f"Closest pose to crystal: pose {best[0]}, {best[1]} kcal/mol, RMSD {best[2]} A")
print(f"Poses within {THRESHOLD} A: {n_good} of {len(rows)}")
print("CHECK 4 - Top-ranked pose reproduces crystal pose:",
      "PASS" if top[2] < THRESHOLD else "FAIL")
print("CHECK 5 - Any pose reproduces crystal pose:",
      "PASS" if best[2] < THRESHOLD else "FAIL")
print("Table saved to results/validation_results.csv")