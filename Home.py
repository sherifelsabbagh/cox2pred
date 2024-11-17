import streamlit as st

st.set_page_config(
    page_title="Home")

st.image("https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcSByQqTdmeNAG0Fhb7TAVN2X8BM9BOX6g0A0g&s")

st.image("https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcT9i6jUsvad581B5oqj4CQUhVUkeV6zdTaHbw&s")

st.title('COX2 Activity Prediction App')

st.info('This App can be used to predict whether a molecule would be active or inactive against COX2 target protein.')

st.subheader("App Usage:")
st.write("The sidebar of the app contains three pages:")
st.write("The first page can be used for prediction of multiple molecules provided in SMILES format (Virtual screening).")
st.write("The second page contains information regarding the built model performance.")
st.write("You can find information about COX2 protein in the third page.")

st.info("© Faculty of Pharmacy, Galala university, 2024")




