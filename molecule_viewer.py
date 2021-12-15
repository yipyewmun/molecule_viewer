# Import dependencies
from pandas.io.parsers import read_csv
import streamlit as st
import streamlit.components.v1 as components
from streamlit import caching
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from rdkit import DataStructs


@st.cache(allow_output_mutation=True)
def get_data():
    return []

sideb=st.sidebar

# Set header title
st.title('Molecule Viewer')

# Compare molecules
compare_molecules = sideb.checkbox('Compare Molecule Structures?')
if compare_molecules:
    with open('tmp/list.txt', 'w') as f:

        uploaded_file = sideb.file_uploader("Choose a CSV file that contains the SMILES_1 & SMILES_2 columns")
        if uploaded_file is not None:
            df = read_csv(uploaded_file)
            smiles_1 = df['SMILES_1'].tolist()
            smiles_2 = df['SMILES_2'].tolist()
            
            # selected_smiles_1 = sideb.selectbox('Select SMILES to visualize', smiles_1)
            selected_smiles_1 = sideb.selectbox('Select SMILES to visualize', range(len(smiles_1)), format_func=lambda x: smiles_1[x])
            selected_smiles_2 = smiles_2[selected_smiles_1]

            mol_1 = Chem.MolFromSmiles(smiles_1[selected_smiles_1])
            d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
            d2d.DrawMolecule(mol_1)
            d2d.FinishDrawing()
            img_1 = d2d.GetDrawingText()
            with open('tmp/mol_image_1.png', 'wb') as png_file:
                png_file.write(img_1)

            mol_2 = Chem.MolFromSmiles(selected_smiles_2)
            d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
            d2d.DrawMolecule(mol_2)
            d2d.FinishDrawing()
            img_2 = d2d.GetDrawingText()
            with open('tmp/mol_image_2.png', 'wb') as png_file:
                png_file.write(img_2)

            tanimoto_similarity = DataStructs.FingerprintSimilarity(Chem.RDKFingerprint(mol_1), Chem.RDKFingerprint(mol_2))
            st.write('Tanimoto Similarity: ' + str(tanimoto_similarity))
            st.image('tmp/mol_image_1.png')
            st.image('tmp/mol_image_2.png')

            add_to_list = sideb.button('Add to list?')
            if add_to_list:
                get_data().append((smiles_1[selected_smiles_1], selected_smiles_2))

            clear_list = sideb.button('Clear List')
            if clear_list:
                download_list = sideb.download_button('Download list when completed', str([]))
                if download_list:
                    st.download_button('Download CSV', str([]))
            else:
                download_list = sideb.download_button('Download list when completed', str(get_data()))
                if download_list:
                    st.download_button('Download CSV', str(get_data()))

        else:
            query_1 = sideb.text_input("Please enter your 1st SMILES here:")
            mol_1 = Chem.MolFromSmiles(query_1)
            d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
            d2d.DrawMolecule(mol_1)
            d2d.FinishDrawing()
            img_1 = d2d.GetDrawingText()

            with open('tmp/mol_image_1.png', 'wb') as png_file:
                png_file.write(img_1)
            
            query_2 = sideb.text_input("Please enter your 2nd SMILES here:")
            mol_2 = Chem.MolFromSmiles(query_2)
            d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
            d2d.DrawMolecule(mol_2)
            d2d.FinishDrawing()
            img_2 = d2d.GetDrawingText()

            with open('tmp/mol_image_2.png', 'wb') as png_file:
                png_file.write(img_2)

            tanimoto_similarity = DataStructs.FingerprintSimilarity(Chem.RDKFingerprint(mol_1), Chem.RDKFingerprint(mol_2))
            st.write('Tanimoto Similarity: ' + str(tanimoto_similarity))
            st.image('tmp/mol_image_1.png')
            st.image('tmp/mol_image_2.png')
else:

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