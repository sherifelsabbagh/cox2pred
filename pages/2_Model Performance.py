import streamlit as st
from PIL import Image

# Page configuration
st.set_page_config(
  page_title='Model performance',
  initial_sidebar_state='expanded')

st.title('Model Peformance')

st.write("We chosed PubChem fingerprints for building the model using RandomForest Classifier. SMOTE technique was applied to overcome imbalance problem.")


with st.expander('Performance', expanded=True):
    st.write("Sensitivity (SN) : 90%")
    st.write("Specificity (SP) : 64%")
    st.write("Matthews’s correlation coefficient (MCC) : 0.52")
    st.write("Accuracy (Q) : 84%")
    Auc_i = Image.open('AUC.png')
    st.image(Auc_i, use_column_width=True)
