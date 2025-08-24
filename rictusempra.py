import streamlit as st
import subprocess
import os
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
from openbabel import pybel
from streamlit_molstar import st_molstar

# --- Helper Functions ---

def smiles_to_2d_image(smiles: str):
    """Generates a 2D image from a SMILES string using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return MolToImage(mol, size=(400, 400))
    except Exception:
        return None

def smiles_to_mol2_file(smiles: str, filename: str) -> str | None:
    """
    Generates a 3D structure, saves it as a MOL2 file, and returns the path.
    """
    try:
        mol = pybel.readstring("smi", smiles)
        mol.make3D()
        mol.write("mol2", filename, overwrite=True)
        return filename
    except Exception:
        return None

def run_dimorphite(smiles: str, min_ph: float, max_ph: float) -> str | None:
    """Runs Dimorphite-DL and returns the protonated SMILES string."""
    try:
        output_filename = 'protonated.smi'
        command = [
            'dimorphite_dl',
            '--ph_min', str(min_ph),
            '--ph_max', str(max_ph),
            '--output_file', output_filename,
            smiles
        ]
        subprocess.run(command, capture_output=True, text=True, check=True)

        if os.path.exists(output_filename):
            with open(output_filename, 'r') as f:
                line = f.readline()
                if line:
                    protonated_smiles = line.split()[0]
                    os.remove(output_filename)
                    return protonated_smiles
    except subprocess.CalledProcessError as e:
        st.error(f"Dimorphite-DL Error: {e.stderr}")
    except Exception as e:
        st.error(f"An error occurred: {e}")
    
    if os.path.exists(output_filename):
        os.remove(output_filename)
        
    return None

def read_file_content(filepath: str) -> str:
    """Reads and returns the content of a file."""
    with open(filepath, 'r') as f:
        return f.read()

# --- Streamlit App ---

st.set_page_config(page_title="Molecule Protonation Tool", layout="wide")
st.title("Rictusempra - Interactive Molecule Protonation Tool 🧪")

# --- Sidebar for Controls ---

with st.sidebar:
    st.image("Rictusempra.png", use_container_width=True)
    st.sidebar.markdown('*Rictusempra* is a web-based cheminformatics tool for interactively visualizing small molecules and calculating their most likely protonation state at a given physiological pH. \n'
                  'It provides a simple interface to generate 2D and 3D molecular structures and prepare them for further computational chemistry tasks like molecular docking or simulation.')
    st.sidebar.markdown('Please see the [documentation](https://github.com/jpmslima/rictusempra) for more information.')
    st.sidebar.markdown('Developed by the [EvoMol-Lab](github.com/evomol-lab).\n'
                    '[BioME](bioinfo.imd.ufrn.br), UFRN, Brazil.')
    st.header("Controls")
    
    smiles_input = st.text_input("Enter SMILES string:", "c1ccccc1C(=O)O")
    
    with st.form("protonation_form"):
        st.write("Set pH for Protonation")
        min_ph = st.number_input("Minimum pH", value=7.2, min_value=0.0, max_value=14.0, step=0.1)
        max_ph = st.number_input("Maximum pH", value=7.6, min_value=0.0, max_value=14.0, step=0.1)
        
        submitted = st.form_submit_button("Calculate Protonation State")

# --- Main Panel for Results ---

if smiles_input:
    # --- Initial Structure Display ---
    st.header("Initial Structure", divider="rainbow")
    
    initial_mol2_path = smiles_to_mol2_file(smiles_input, "initial.mol2")
    img_initial = smiles_to_2d_image(smiles_input)

    if initial_mol2_path and img_initial:
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("2D Structure")
            st.image(img_initial, use_container_width=True)
        with col2:
            st.subheader("3D Structure")
            st_molstar(initial_mol2_path, key="molstar_initial")
            
            mol2_content = read_file_content(initial_mol2_path)
            st.download_button(
                label="Download .mol2 File",
                data=mol2_content,
                file_name="initial_structure.mol2",
                mime="chemical/x-mol2"
            )
    else:
        st.error("Invalid SMILES string. Please check your input.")

    # --- Protonated Structure Display (if calculation was run) ---
    if submitted:
        st.header("Protonated Structure", divider="rainbow")
        with st.spinner("Running Dimorphite-DL..."):
            protonated_smiles = run_dimorphite(smiles_input, min_ph, max_ph)
        
        if protonated_smiles:
            st.success("Protonation calculation complete!")
            st.code(protonated_smiles, language="smiles")
            
            protonated_mol2_path = smiles_to_mol2_file(protonated_smiles, "protonated.mol2")
            img_protonated = smiles_to_2d_image(protonated_smiles)
            
            if protonated_mol2_path and img_protonated:
                col3, col4 = st.columns(2)
                with col3:
                    st.subheader("Protonated 2D Structure")
                    st.image(img_protonated, use_container_width=True)
                with col4:
                    st.subheader("Protonated 3D Structure")
                    st_molstar(protonated_mol2_path, key="molstar_protonated")
                    
                    protonated_mol2_content = read_file_content(protonated_mol2_path)
                    st.download_button(
                        label="Download Protonated .mol2 File",
                        data=protonated_mol2_content,
                        file_name="protonated_structure.mol2",
                        mime="chemical/x-mol2"
                    )
        else:
            st.error("Could not determine the protonated structure.")