import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results/validation_results.csv")

fig, ax = plt.subplots(figsize=(8, 5))
colors = ["seagreen" if r < 2.0 else "firebrick" for r in df["RMSD_A"]]
ax.scatter(df["RMSD_A"], df["Affinity_kcal_mol"], c=colors, s=70, edgecolor="black")

ax.axvline(2.0, color="gray", linestyle="--", label="2.0 A threshold")
for _, row in df.iterrows():
    if row["Pose"] in (1, 4):
        ax.annotate(f"pose {int(row['Pose'])}",
                    (row["RMSD_A"], row["Affinity_kcal_mol"]),
                    textcoords="offset points", xytext=(8, 6))

ax.set_xlabel("RMSD to crystal ligand (A)")
ax.set_ylabel("Vina affinity (kcal/mol)")
ax.set_title("Imatinib redocking into Abl kinase (1IEP)")
ax.legend()
plt.tight_layout()
plt.savefig("figures/affinity_vs_rmsd.png", dpi=300)
print("Saved figures/affinity_vs_rmsd.png")