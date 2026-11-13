from rdkit import Chem
from rdkit.Chem import rdMolAlign, AllChem

# Reference: crystal imatinib (chain A), heavy atoms only
ref_raw = Chem.MolFromPDBFile("ligand/crystal_A.pdb", removeHs=True)

# Template with correct bond orders, from the original SDF (first fragment only)
suppl = Chem.SDMolSupplier("ligand/imatinib_original.sdf", removeHs=True)
frags = Chem.GetMolFrags(suppl[0], asMols=True)
template = frags[0]

ref = AllChem.AssignBondOrdersFromTemplate(template, ref_raw)

print("Pose  Affinity-rank  RMSD (A)")
results = []
for i in range(1, 10):
    pose_raw = Chem.MolFromPDBFile(f"results/pose{i}.pdb", removeHs=True)
    pose = AllChem.AssignBondOrdersFromTemplate(template, pose_raw)
    # No alignment: compare the docked pose directly in the protein frame
    rmsd = rdMolAlign.CalcRMS(pose, ref)
    results.append((i, rmsd))
    print(f"{i:>4}  {i:>13}  {rmsd:8.3f}")

best = min(results, key=lambda x: x[1])
print()
print(f"Top-ranked pose RMSD: {results[0][1]:.3f} A")
print(f"Closest pose: pose {best[0]} with RMSD {best[1]:.3f} A")
print("Docking SUCCESS" if results[0][1] < 2.0 else "Top pose above 2.0 A threshold")