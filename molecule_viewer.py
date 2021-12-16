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

@st.cache(allow_output_mutation=True)
def get_idx():
    return [0]

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
            selected_smiles_1 = sideb.selectbox('Select Molecule Pair to visualize (' + str(len(smiles_1)) + ' pairs in total)' , range(len(smiles_1)), format_func=lambda x: x+1)
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
            st.write(str(selected_smiles_1+1) + '. Tanimoto Similarity: ' + str(tanimoto_similarity))

            col1, col2 = st.columns([1,1])
            with col1:
                st.image('tmp/mol_image_1.png')
            with col2:
                st.image('tmp/mol_image_2.png')

            add_to_list = sideb.button('Add to list?')
            if add_to_list:
                if (smiles_1[selected_smiles_1], selected_smiles_2) in get_data():
                    pass
                else:
                    get_data().append((smiles_1[selected_smiles_1], selected_smiles_2))

            clear_list = sideb.button('Clear List')
            if clear_list:
                get_data().clear()

            if len(get_data()) > 0:
                download_list = sideb.download_button('Download list when completed ' + '(' + str(len(get_data())) + ' items)', str(get_data()))

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

        # selected_SMILES = sideb.selectbox('Select SMILES to visualize', smiles_lst)
        # mol = Chem.MolFromSmiles('')
        # d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
        # d2d.DrawMolecule(mol)
        # d2d.FinishDrawing()
        # img = d2d.GetDrawingText()

        col1, col2 = st.columns([1,1])
        with col1:
            previous_structure = st.button('Previous')
        with col2:
            next_structure = st.button('Next')
        
        if previous_structure:
            # Get previous index
            prev_idx = get_idx()[-1] - 1
            if prev_idx < 0:
                pass
            elif prev_idx == 0:
                mol = Chem.MolFromSmiles(smiles_lst[prev_idx])
                d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
                d2d.DrawMolecule(mol)
                d2d.FinishDrawing()
                img = d2d.GetDrawingText()

                print(prev_idx)
                get_idx().append(prev_idx)
                print(get_idx())
            else:
                mol = Chem.MolFromSmiles(smiles_lst[prev_idx])
                d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
                d2d.DrawMolecule(mol)
                d2d.FinishDrawing()
                img = d2d.GetDrawingText()

                # print(prev_idx)
                # get_idx().append(prev_idx)
                # print(get_idx())

        if next_structure:
            # Get next index
            next_idx = get_idx()[-1] + 1

            mol = Chem.MolFromSmiles(smiles_lst[next_idx])
            d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
            d2d.DrawMolecule(mol)
            d2d.FinishDrawing()
            img = d2d.GetDrawingText()

            # print(next_idx)
            # get_idx().append(next_idx)
            # print(get_idx())
        
        with open('tmp/mol_image.png', 'wb') as png_file:
            png_file.write(img)

        st.image('tmp/mol_image.png')
    else:
        get_idx().clear()
        get_idx().append(0)

        query = sideb.text_input("Please enter your SMILES here:")
        mol = Chem.MolFromSmiles(query)
        d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        img = d2d.GetDrawingText()

        with open('tmp/mol_image.png', 'wb') as png_file:
            png_file.write(img)

        st.image('tmp/mol_image.png')