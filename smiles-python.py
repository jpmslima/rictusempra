from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

def draw_molecule_from_smiles(smiles_string):
    """
    Interprets a SMILES string and returns a Pillow Image of the chemical structure.

    Args:
        smiles_string (str): The SMILES string to interpret.

    Returns:
        PIL.Image.Image: A Pillow Image object of the chemical structure.
                         Returns None if the SMILES string is invalid.
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is not None:
            img = MolToImage(mol)
            return img
        else:
            return None
    except Exception as e:
        print(f"Error processing SMILES: {e}")
        return None

if __name__ == '__main__':
    # Example usage
    test_smiles = 'CCO' # Ethanol
    image = draw_molecule_from_smiles(test_smiles)
    if image:
        image.save('ethanol.png')
        print(f"Image of {test_smiles} saved as ethanol.png")
    else:
        print(f"Could not draw molecule for {test_smiles}")
