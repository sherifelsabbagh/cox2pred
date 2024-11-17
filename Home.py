import streamlit as st

st.set_page_config(
    page_title="Home")

st.title('LRRK2 Activity prediction App')
st.info('The LRRK2 Activity Prediction App can be used to predict whether a molecule is active or inactive for LRRK2 target protein .')

st.subheader("App function:")
st.write("The sidebar of App Contain four function page")
st.write("The first page can be used for prediction of only one molecule. It should be provided in SMILES format.")
st.write("The second page is used for prediction of multiple molecules provided in SMILES format .")
st.write("The third page has information about model performance.")
st.write("The last page has information about LRRK2 protein.")

st.warning('Limit 250 compounds per file')
