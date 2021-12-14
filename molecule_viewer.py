# Import dependencies
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG


sideb=st.sidebar

# Set header title
st.title('Molecule Viewer')

# Read dataset (CSV)
uploaded_file = sideb.file_uploader("Choose a CSV file that contains a SMILES column")

if uploaded_file is not None:
    df = pd.read_csv(uploaded_file)
    smiles_lst = df['SMILES'].tolist()

    selected_SMILES = sideb.selectbox('Select SMILES to visualize', smiles_lst)
    mol = Chem.MolFromSmiles(selected_SMILES)
    d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    img = d2d.GetDrawingText()

    with open('tmp/mol_image.png', 'wb') as png_file:
        png_file.write(img)

    st.image('tmp/mol_image.png')
else:
    query = sideb.text_input("Please enter your SMILES here:")
    mol = Chem.MolFromSmiles(query)
    d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    img = d2d.GetDrawingText()

    with open('tmp/mol_image.png', 'wb') as png_file:
        png_file.write(img)

    st.image('tmp/mol_image.png')