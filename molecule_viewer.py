# Import dependencies
from pandas.io.parsers import read_csv
import streamlit as st
import streamlit.components.v1 as components
from streamlit import caching
from bokeh.models.widgets import Div
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

def generate_image(smiles):
    mol = Chem.MolFromSmiles(smiles)
    d2d = rdMolDraw2D.MolDraw2DCairo(350,300)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    img = d2d.GetDrawingText()

    return img

def customButton_searchPubChem(button_text, link):
    if st.button('Search Molecule in PubChem'):
        js = "window.open('" + link + "')"  # New tab or window
        html = '<img src onerror="{}">'.format(js)
        div = Div(text=html)
        st.bokeh_chart(div)


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
            
            selected_smiles_1 = sideb.selectbox('Select Molecule Pair to visualize (' + str(len(smiles_1)) + ' pairs in total)' , range(len(smiles_1)), format_func=lambda x: x+1)
            selected_smiles_2 = smiles_2[selected_smiles_1]

            img_1 = generate_image(smiles_1[selected_smiles_1])
            with open('tmp/mol_image_1.png', 'wb') as png_file:
                png_file.write(img_1)

            img_2 = generate_image(selected_smiles_2)
            with open('tmp/mol_image_2.png', 'wb') as png_file:
                png_file.write(img_2)

            tanimoto_similarity = DataStructs.FingerprintSimilarity(Chem.RDKFingerprint(Chem.MolFromSmiles(smiles_1[selected_smiles_1])), Chem.RDKFingerprint(Chem.MolFromSmiles(selected_smiles_2)))
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
            img_1 = generate_image(query_1)
            with open('tmp/mol_image_1.png', 'wb') as png_file:
                png_file.write(img_1)
            
            query_2 = sideb.text_input("Please enter your 2nd SMILES here:")
            img_2 = generate_image(query_2)
            with open('tmp/mol_image_2.png', 'wb') as png_file:
                png_file.write(img_2)

            tanimoto_similarity = DataStructs.FingerprintSimilarity(Chem.RDKFingerprint(Chem.MolFromSmiles(query_1)), Chem.RDKFingerprint(Chem.MolFromSmiles(query_2)))
            st.write('Tanimoto Similarity: ' + str(tanimoto_similarity))
            st.image('tmp/mol_image_1.png')
            st.image('tmp/mol_image_2.png')
else:

    # Read dataset (CSV)
    uploaded_file = sideb.file_uploader("Choose a CSV file that contains a SMILES column")

    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        smiles_lst = df['SMILES'].tolist()

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
                if smiles_lst[prev_idx] == "-":
                    st.image('tmp/No or Missing Structure.png')
                else:
                    img = generate_image(smiles_lst[prev_idx])
                    with open('tmp/mol_image.png', 'wb') as png_file:
                        png_file.write(img)

                    st.image('tmp/mol_image.png')
                    customButton_searchPubChem('Search Molecule in PubChem', 'https://pubchem.ncbi.nlm.nih.gov/#query=' + smiles_lst[prev_idx])

                get_idx().append(prev_idx)

            else:
                if smiles_lst[prev_idx] == "-":
                    st.image('tmp/No or Missing Structure.png')
                else:
                    img = generate_image(smiles_lst[prev_idx])
                    with open('tmp/mol_image.png', 'wb') as png_file:
                        png_file.write(img)

                    st.image('tmp/mol_image.png')
                    customButton_searchPubChem('Search Molecule in PubChem', 'https://pubchem.ncbi.nlm.nih.gov/#query=' + smiles_lst[prev_idx])

                get_idx().append(prev_idx)

        elif next_structure:
            # Get next index
            next_idx = get_idx()[-1] + 1
            if next_idx == len(smiles_lst):
                pass
            else:
                if smiles_lst[next_idx] == "-":
                    st.image('tmp/No or Missing Structure.png')
                else:
                    img = generate_image(smiles_lst[next_idx])
                    with open('tmp/mol_image.png', 'wb') as png_file:
                        png_file.write(img)

                    st.image('tmp/mol_image.png')
                    customButton_searchPubChem('Search Molecule in PubChem', 'https://pubchem.ncbi.nlm.nih.gov/#query=' + smiles_lst[next_idx])

                get_idx().append(next_idx)

        else:
            if smiles_lst[0] == "-":
                st.image('tmp/No or Missing Structure.png')
            else:
                img = generate_image(smiles_lst[0])
                with open('tmp/mol_image.png', 'wb') as png_file:
                    png_file.write(img)

                st.image('tmp/mol_image.png')
                customButton_searchPubChem('Search Molecule in PubChem', 'https://pubchem.ncbi.nlm.nih.gov/#query=' + smiles_lst[0])
        
    else:
        get_idx().clear()
        get_idx().append(0)

        query = sideb.text_input("Please enter your SMILES here:")
        img = generate_image(query)
        with open('tmp/mol_image.png', 'wb') as png_file:
            png_file.write(img)

        st.image('tmp/mol_image.png')