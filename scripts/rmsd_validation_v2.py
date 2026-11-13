import csv
from rdkit import Chem
from rdkit.Chem import rdMolAlign, AllChem

ref_raw = Chem.MolFromPDBFile("ligand/crystal_A.pdb", removeHs=True)
suppl = Chem.SDMolSupplier("ligand/imatinib_original.sdf", removeHs=True)
template = Chem.GetMolFrags(suppl[0], asMols=True)[0]
ref = AllChem.AssignBondOrdersFromTemplate(template, ref_raw)

affinities = []
for line in open("results/docked_v2.pdbqt"):
    if line.startswith("REMARK VINA RESULT"):
        affinities.append(float(line.split()[3]))

rows = []
print("Pose  Affinity  RMSD (A)")
for i, aff in enumerate(affinities, start=1):
    raw = Chem.MolFromPDBFile(f"results/v2_pose{i}.pdb", removeHs=True)
    pose = AllChem.AssignBondOrdersFromTemplate(template, raw)
    rmsd = rdMolAlign.CalcRMS(pose, ref)
    rows.append((i, aff, round(rmsd, 3)))
    print(f"{i:>4}  {aff:>8}  {rmsd:8.3f}")

with open("results/rmsd_v2.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Pose", "Affinity_kcal_mol", "RMSD_A"])
    w.writerows(rows)

best = min(rows, key=lambda r: r[2])
print()
print(f"Top-ranked pose RMSD: {rows[0][2]} A")
print(f"Closest pose: pose {best[0]} with RMSD {best[2]} A")
print("Docking SUCCESS" if rows[0][2] < 2.0 else "Top pose above 2.0 A threshold")