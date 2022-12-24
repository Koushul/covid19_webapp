from typing import Tuple, Union
import streamlit as st
from utils.api import load_tftargets
from texts.references import CELL_NAMES
from utils.components import create_st_button
from utils.dataloaders import Dataset
import pandas as pd

def display_tf_target_data(option: str, cell_type: str, ds: Dataset, disabled = False) -> Tuple[str, pd.DataFrame]:
    if not disabled:
        st.markdown(f'### 🎯 {option}-Target genes correlation in {CELL_NAMES.get(cell_type, cell_type)}') 
        
    with st.spinner('Loading heatmap...'):
        try:
            img = ds.get_tf_target_gene_img(cell_type, option)
            if not disabled:
                st.image(img)
        except:
            if not disabled:
                st.warning('No high correlation matches found.')
            
    if not disabled: st.markdown('---')
    targets = load_tftargets(tf=option)
    if not disabled: st.subheader(f'{option} target genes based on hTFtarget (N={len(targets)})')
    
    left, right = st.columns(2)
    
    
    if len(targets) > 0 and not disabled:
        with left:
            target_name = st.selectbox('Targets', targets.target_name.values)
            target_ensg = targets[targets.target_name==target_name].target_id.values[0]
        
        with left:
            _s = '!'
            create_st_button(
                link_text=f'View epigenomic states for {target_name}', 
                link_url=f'http://bioinfo.life.hust.edu.cn/hTFtarget#{_s}/epigenomic_states?gene={target_ensg}&tf={option}&scale=2')
        
        st.table(targets)
        
        
    return img, targets