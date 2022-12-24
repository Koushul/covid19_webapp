import streamlit as st

from utils.components import ridge_plot
import hydralit_components as hc

import time
import pandas as pd

def viral_entry_page(datacore):
    """
    protein-TFs_higCorr_viral/Figure
    mRNA-TFs_higCorr_viral/Figure
    generate density plots for mRNA(viral entry) from mRNAs_all_…viral.csv and get pvalue from mRNA_wilcox…_viral.csv
    generate density plots for protein(viral entry) from proteinss_all_…viral.csv and get pvalue from protein_wilcox…_viral.csv
        some of the data file may not exist sine no viral protein or genes for that dataset or celltype
    
    """
    a, b = st.columns(2)
    
    datasets = datacore.name2ds.keys()
    

    a.markdown('##### &nbsp;')            
    sbs = b.checkbox('Side-by-Side Plots', value=True)    
        
    p_ds = a.selectbox(f'Dataset ({len(datasets)})', datasets, 0)
    cell_types = datacore.name2ds[p_ds].celltypes
    cell_type = b.selectbox(f'Cell Type ({len(cell_types)})', cell_types, 0)
    
    s_proteins = list(datacore.name2ds[p_ds].get_viral_proteins(cell_type))
    s_mrna = list(datacore.name2ds[p_ds].get_viral_mrnas(cell_type))
    

    
    show_protein_tf = b.checkbox('Viral-entry related factors-TF correlations')
    show_mrna_tf = b.checkbox('Viral-entry related factors mRNA-TF correlations')
    
    show_protein_density  = a.checkbox('Viral-entry related factors protein distributions')
    show_mrna_density = a.checkbox('Viral-entry related factors mRNA expression distributions')
    
    
    
    st.markdown('---')  

    if sbs:
        if show_protein_tf or show_mrna_tf:
            st.caption('Click on image corner to enlarge')
        
        columns = st.columns(len(s_proteins))
        columns_2 = st.columns(len(s_mrna))
        
        
        for c, p in zip(columns, s_proteins):
            c.markdown('##### ' + p)
          
        if show_protein_tf:
            # st.markdown('#### Viral protein-TF correlations')
            
            for cols, option_protein in zip(columns, s_proteins):
                with cols:
                    with hc.HyLoader(f'🎨 Painting {option_protein}', hc.Loaders.standard_loaders,index=[3]):
                        csv, img = datacore.name2ds[p_ds].get_viral_protein_data(cell_type, option_protein)
                        st.image(img)
                        
        st.markdown('---')
            
            
        if show_protein_density:
            for cols, protein in zip(columns, s_proteins):
                with cols:
                    with hc.HyLoader(f'🎨 Painting {protein}', hc.Loaders.standard_loaders,index=[3]):
                        df = datacore.name2ds[p_ds].get_viral_protein_df(cell_type, protein)
                        img = ridge_plot(df, xlabel=f'{protein} expression', column=None)
                        xx, yy = img.size
                        dd = 2
                        xx = int(xx/dd)
                        yy = int(yy/dd)
                        st.image(img.resize((xx, yy)))
        st.markdown('---')
        
        if show_mrna_tf:
            # st.markdown('#### Viral protein-TF correlations')
            
            for c, p in zip(columns_2, s_mrna):
                c.markdown('##### ' + p)
            for cols, option_protein in zip(columns_2, s_mrna):
                with cols:
                    with hc.HyLoader(f'🎨 Painting {option_protein}', hc.Loaders.standard_loaders,index=[3]):
                        csv, img = datacore.name2ds[p_ds].get_viral_mrna_data(cell_type, option_protein)
                        st.image(img)
                        
        st.markdown('---')
            

            
        if show_mrna_density:
            for cols, option_protein in zip(columns_2, s_mrna):
                with cols:
                    with hc.HyLoader(f'🎨 Painting {option_protein}', hc.Loaders.standard_loaders,index=[3]):
                        df = datacore.name2ds[p_ds].get_viral_mrna_df(cell_type, option_protein)
                        img = ridge_plot(df, xlabel=f'{option_protein} expression', column=None)
                        xx, yy = img.size
                        dd = 2
                        xx = int(xx/dd)
                        yy = int(yy/dd)
                        st.image(img.resize((xx, yy)))
            
            
    
    else:
        if show_protein_tf:
            for option_protein in s_proteins:
                st.markdown('##### '+option_protein)
                
                with hc.HyLoader(f'🎨 Painting {option_protein}', hc.Loaders.standard_loaders,index=[3]):
                    csv, img = datacore.name2ds[p_ds].get_viral_protein_data(cell_type, option_protein)
                    st.image(img)
                    
        st.markdown('---')
        
        if show_protein_density:
            for protein in s_proteins:
                with hc.HyLoader(f'🎨 Painting {protein}', hc.Loaders.standard_loaders,index=[3]):
                    df = datacore.name2ds[p_ds].get_viral_protein_df(cell_type, protein)
                    img = ridge_plot(df, xlabel=f'{protein} expression', column=None)
                    xx, yy = img.size
                    dd = 2
                    xx = int(xx/dd)
                    yy = int(yy/dd)
                    st.image(img.resize((xx, yy)))
        
        if show_mrna_tf:
            
            for option_protein in s_mrna:
                st.markdown('##### ' + option_protein)

                with hc.HyLoader(f'🎨 Painting {option_protein}', hc.Loaders.standard_loaders,index=[3]):
                    csv, img = datacore.name2ds[p_ds].get_viral_mrna_data(cell_type, option_protein)
                    st.image(img)
                    
        st.markdown('---')
        
        if show_mrna_density:
            for mrna in s_mrna:
                with hc.HyLoader(f'🎨 Painting {mrna}', hc.Loaders.standard_loaders,index=[3]):
                    
                    df = datacore.name2ds[p_ds].get_viral_mrna_df(cell_type, mrna)
                    img = ridge_plot(df, xlabel=f'{mrna} expression', column=None)
                    xx, yy = img.size
                    dd = 2
                    xx = int(xx/dd)
                    yy = int(yy/dd)
                    st.image(img.resize((xx, yy)))
                    
