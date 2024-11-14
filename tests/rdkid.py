from rdkit.Chem import Draw
from rdkit import Chem

# Create a list of molecules
smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

# Draw the molecules in a grid
img = Draw.MolsToGridImage(molecules, molsPerRow=2, subImgSize=(200, 200), returnPNG=False)
img.show()