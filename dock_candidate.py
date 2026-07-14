import sys
import subprocess

def dock_smiles(smiles, name="candidate"):
    # Write the SMILES to a temporary file
    with open(f"{name}.smi", "w") as f:
        f.write(smiles + "\n")

    # Generate a 3D structure and convert to PDBQT in one step
    subprocess.run([
        "obabel", f"{name}.smi", "-O", f"{name}.pdbqt", "--gen3d"
    ], check=True)

    # Dock against the EGFR pocket (same coordinates as the erlotinib validation)
    subprocess.run([
        "./vina",
        "--receptor", "receptor.pdbqt",
        "--ligand", f"{name}.pdbqt",
        "--center_x", "22.014", "--center_y", "0.253", "--center_z", "52.794",
        "--size_x", "20", "--size_y", "20", "--size_z", "20",
        "--out", f"{name}_docked.pdbqt"
    ])

if __name__ == "__main__":
    smiles = sys.argv[1]
    name = sys.argv[2] if len(sys.argv) > 2 else "candidate"
    dock_smiles(smiles, name)
