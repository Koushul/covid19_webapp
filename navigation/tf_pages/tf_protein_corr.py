import streamlit as st
import pandas as pd
from utils.dataloaders import Dataset
from texts.references import CELL_NAMES

def display_tf_protein_corr_data(option: str, cell_type: str, ds: Dataset, disabled = False):
    if not disabled: st.markdown(f'### 🔗 Correlation between {option} and surface protein expression in {CELL_NAMES.get(cell_type, cell_type)}')
    
    try:
        csv, img = ds.get_tf_protein_data(cell_type, option)
        if not disabled: 
            st.image(img)
            exp = st.expander('Raw Data')
            exp.dataframe(pd.read_csv(csv))
        return csv, img
        
    except IndexError:
        st.warning('No high correlation matches found.')
        return None, None
        
        