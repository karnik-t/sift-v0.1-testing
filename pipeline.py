import sys
import subprocess
import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
model = joblib.load("egfr_model.pkl")

def smiles_to_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return np.array(fp_gen.GetFingerprint(mol))

def largest_fragment(smiles):
    return max(smiles.split("."), key=len)

def dock_smiles(smiles, name):
    smiles = largest_fragment(smiles)
    with open(f"{name}.smi", "w") as f:
        f.write(smiles + "\n")
    subprocess.run(["obabel", f"{name}.smi", "-O", f"{name}.pdbqt", "--gen3d"],
                    check=True, capture_output=True)
    result = subprocess.run([
        "./vina", "--receptor", "receptor.pdbqt", "--ligand", f"{name}.pdbqt",
        "--center_x", "22.014", "--center_y", "0.253", "--center_z", "52.794",
        "--size_x", "20", "--size_y", "20", "--size_z", "20",
        "--out", f"{name}_docked.pdbqt"
    ], capture_output=True, text=True)
    for line in result.stdout.splitlines():
        if line.strip().startswith("1"):
            return float(line.split()[1])
    return None

def run_pipeline(csv_path, top_n):
    df = pd.read_csv(csv_path)
    df["fingerprint"] = df["smiles"].apply(smiles_to_fp)
    df = df.dropna(subset=["fingerprint"])
    X = np.stack(df["fingerprint"].values)
    df["predicted_prob_active"] = model.predict_proba(X)[:, 1]
    top = df.sort_values("predicted_prob_active", ascending=False).head(top_n)

    print(f"Screened {len(df)} candidates. Docking top {top_n}...\n")
    results = []
    for _, row in top.iterrows():
        safe_name = "".join(c if c.isalnum() else "_" for c in str(row["name"]))[:30]
        score = dock_smiles(row["smiles"], safe_name)
        results.append((row["name"], row["predicted_prob_active"], score))
        print(f"{row['name']:30s}  screening={row['predicted_prob_active']:.3f}  docking={score}")

    return results

if __name__ == "__main__":
    csv_path = sys.argv[1]
    top_n = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    run_pipeline(csv_path, top_n)
