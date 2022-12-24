import streamlit as st
from texts.descriptions import DatasetDesc
from texts.references import CELL_NAMES
from utils.components import resize, ridge_plot

def display_protein_expressions(option: str, ds: DatasetDesc, cell_type: str, disabled=False):
    try:
        if not disabled: 
            st.markdown(f'### 📈 {option} expression in {CELL_NAMES.get(cell_type, cell_type)}')  
        
        df_plt = ds.get_protein_df(cell_type, option)
        
        groups = list(set(df_plt['patient_group'].values))

        if 'healthy' in groups:
            good = 'healthy'
            bad = ', '.join(list(set(groups).difference([good])))
            bad = f'impacted ({bad})'
        elif 'Healthy' in groups:
            good = 'Healthy'
            bad = ', '.join(list(set(groups).difference([good])))
            bad = f'Impacted ({bad})'
            
        else:
            good = groups[0]
            bad = groups[1]
        
        
        pval, adj_pval = ds.get_protein_pvals(cell_type, option)
        if pval < 0.05:
            st.success(f'{good} vs {bad}: p-value = {pval:.3f}')
        else:
            st.error(f'{good} vs {bad}: p-value = {pval:.3f}')
        
        
        c, d, _ = st.columns([2, 6, 2])
        if disabled:
            img = ridge_plot(df_plt, f'{option} expression', column=None)
        else:
            img = ridge_plot(df_plt, f'{option} expression', column=c)
        img = resize(img)
        if not disabled:
            d.image(img)
            
        return img
    
    except (FileNotFoundError, KeyError):
        st.warning('No high correlation matches found.')
        return None
    