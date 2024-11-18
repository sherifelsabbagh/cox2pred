import streamlit as st
import pandas as pd
import subprocess
import os
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

st.set_page_config(
    page_title="Virtual Screening")


st.header("Virtual Screening")
st.warning("""
[Please be sure that file format such this](https://drive.google.com/file/d/1Qa-0JKVa9ccuSMd4Km6Z5nP-RLeT1Gp0/view?usp=sharing)
""")

def PUbchemfp_desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints  -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')
    
def descriptors(smiles):
    mols = [Chem.MolFromSmiles(s) for s in smiles] 
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    desc_names = calc.GetDescriptorNames()
    
    Mol_descriptors =[]
    for mol in mols:

        descriptors = calc.CalcDescriptors(mol)
        Mol_descriptors.append(descriptors)
    return Mol_descriptors,desc_names 
    

# Model
def the_model(input_data):
    load_model = joblib.load(rf_model.pkl')
    # Make prediction
    prediction = load_model.predict(input_data)
    prediction_probability=load_model.predict_proba(input_data)
    
    x=pd.DataFrame(prediction_probability,columns=["Inactive probability","Active_probability"])
    st.header('Prediction Result')
    
    prediction_output = pd.Series(prediction, name='Activity')
    #proba_output=pd.Series(prediction_probability,name="prediction_proba")
    
    molecule_name = pd.Series(reading_data[1], name='Molecule CHEMBL id/Molecule Name ')
    
    Result= pd.concat([molecule_name, prediction_output,x], axis=1)
    
    result = []
    for x in Result["Activity"]:
        if x==1:
            result.append("Active")
        if x==0:
            result.append("Inactive")
    Result["Activity"]=result
    st.write(Result)
    prediction_csv = Result.to_csv(index=False)
    st.download_button(label="Download prediction result",data=prediction_csv,file_name="My_result.csv")



uplouded_file=st.file_uploader("Please upload your input file", type=['txt'])


if st.button('Predict'):
    reading_data = pd.read_table(uplouded_file, sep='\t', header=None)
    reading_data.to_csv('molecule.smi', sep = '\t', index = False,header=None)
    st.subheader('The input data')
    st.write(reading_data)



    PUbchemfp_desc_calc()

 
    st.subheader('Calculated Pubchem_Fingerprint descriptors')
    pubfp_calc = pd.read_csv("descriptors_output.csv")
    pubfp_calc.drop('Name', axis=1, inplace=True)
    st.write(pubfp_calc)
    st.write(pubfp_calc.shape)
    # Read descriptor feature used in training
    st.header('Descriptors Subset ')
    feature_list = list(pd.read_csv('Pubchem_features.csv').columns)
    desc_subset = pubfp_calc[feature_list]
    st.write(desc_subset)
    st.write(desc_subset.shape)
    #st.download_button(label="Download PubchemFingerprinter subset Descriptor",data=desc_subset,file_name="Descriptor_subset.csv")
    first_column = reading_data.iloc[:, 0] 
    Mol_descriptors,desc_names =descriptors(first_column)
    df_with_200_descriptors = pd.DataFrame(Mol_descriptors,columns=desc_names)
    df=df_with_200_descriptors[["MolWt","MolLogP","NumHAcceptors","NumHDonors"]]
    st.subheader("Lipinski Rule of 5 Descriptors")
    st.write(df)
    
    
    the_model(desc_subset)
else:
    st.warning('Limit 250 compounds per file')
    
    
