import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolToImage
from streamlit_molstar import st_molstar

# --- Functions ---

def draw_molecule_from_smiles(smiles_string):
    """Generates a 2D image of a molecule from a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        return MolToImage(mol, size=(500, 500)) if mol else None
    except Exception:
        return None

def get_3d_mol_block(smiles_string):
    """Generates a 3D MolBlock from a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None: return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=1)
        AllChem.UFFOptimizeMolee(mol)
        return Chem.MolToMolBlock(mol)
    except Exception:
        return None

# --- Streamlit App Layout ---

st.set_page_config(page_title="SMILES to Structure Viewer 🧪", layout="centered")

st.title("SMILES to Structure Viewer")
st.markdown("Enter a SMILES string below to visualize its 2D and 3D chemical structures.")

smiles_input = st.text_input("SMILES string:", 'O=C(O)c1ccccc1C(=O)O') # Phthalic acid

if smiles_input:
    tab1, tab2 = st.tabs(["2D Structure", "3D Structure"])

    with tab1:
        st.subheader("2D Visualization")
        img = draw_molecule_from_smiles(smiles_input)
        if img:
            st.image(img, caption=f"2D Structure of {smiles_input}", use_container_width=True)
        else:
            st.error("Invalid SMILES string. Please check your input.")

    with tab2:
        st.subheader("Interactive 3D Visualization")
        mol_block = get_3d_mol_block(smiles_input)
        if mol_block:
            # CORRECTED: Pass the data using the 'files' argument
            st_molstar(files=[{'data': mol_block, 'extension': 'mol'}], key="molstar_viewer", height=500)
        else:
            st.error("Could not generate a 3D model for the given SMILES string.")
else:
    st.warning("Please enter a SMILES string to get started.")