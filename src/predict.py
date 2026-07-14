import sys
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
model = joblib.load("egfr_model.pkl")

def predict(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = np.array(fp_gen.GetFingerprint(mol)).reshape(1, -1)
    return model.predict_proba(fp)[0][1]

if __name__ == "__main__":
    smiles = sys.argv[1]
    prob = predict(smiles)
    if prob is None:
        print("Invalid SMILES string")
    else:
        print(f"Predicted probability of EGFR activity: {prob:.3f}")

