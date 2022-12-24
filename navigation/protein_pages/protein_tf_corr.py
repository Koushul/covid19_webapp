import streamlit as st
import pandas as pd
from utils.dataloaders import Dataset
from texts.references import CELL_NAMES


def display_protein_tf_corr_data(option: str, cell_type: str, ds: Dataset, disabled = False):
    if not disabled: 
        st.markdown(f'### 🔗 Correlation between {option} expression and inferred TF activities in {CELL_NAMES.get(cell_type, cell_type)}')
        
    try:
        csv, img = ds.get_protein_tf_data(cell_type, option)            
        if not disabled:
            st.image(img)
            exp = st.expander('Raw Data')            
            exp.dataframe(pd.read_csv(csv))
        return csv, img
    
    
    except IndexError:
        if not disabled: 
            st.warning('No high correlation matches found.')
        return None, None