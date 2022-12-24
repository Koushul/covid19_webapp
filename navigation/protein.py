import random
import pandas as pd
import streamlit as st
import hydralit_components as hc
from navigation.protein_pages.expression import display_protein_expressions
from navigation.protein_pages.litcovid import display_litcovid_data
from navigation.protein_pages.protein_tf_corr import display_protein_tf_corr_data

from utils.api import litcov_search
from utils.components import display_ppi_data
from io import BytesIO
from PIL import Image

def generate_color():
    random_number = random.randint(0,16777215)
    hex_number = str(hex(random_number))
    hex_number ='#'+ hex_number[2:]
    return hex_number
   


protein_bar = {'txc_inactive': 'black','menu_background':'white','txc_active':'white','option_active':'blue'}

def download_button(desc: str, data, fname: str, label: str = 'Download', format: str ='text/csv'):
    left_col, right_col = st.columns([2, 1])
    right_col.download_button(
        label='📥 '+label,
        data=data,
        file_name=fname,
        mime=format,
    )
        
    left_col.write(desc)

def img_path2buffer(img_path: str):
    img = Image.open(img_path)
    buf = BytesIO()
    img.save(buf, format='PNG')
    png = buf.getvalue()
    return png



def protein_page(datacore):
    st.info('Search surface protein across COVID-19 patients and health individuals. For a single surface protein, SPaRTAN COVID-19db supports exploring its expression for each cell type, correlated transcription factors, protein-protein interactions and relevant literature.')
    
    _, _, edge = st.columns(3)
    
    checkbox = edge.checkbox('Show high correlations only')    
    a, b, c = st.columns(3)
    st.session_state['datasets'] = ["GSE155673"]
    datasets = datacore.name2ds.keys()
    p_ds = a.selectbox(f'Dataset ({len(datasets)})', datasets, 0)
    cell_types = datacore.name2ds[p_ds].celltypes
    cell_type = b.selectbox(f'Cell Type ({len(cell_types)})', cell_types, 0)
    
    if not checkbox:
        s_proteins = list(datacore.name2ds[p_ds].total_surface_proteins)
    else:
        s_proteins = list(datacore.name2ds[p_ds].surface_proteins(cell_type))

    selected_protein = c.selectbox(
        f"Surface protein of interest ({len(s_proteins)})",
        s_proteins)
            
    try:    
        litcovdf = litcov_search(selected_protein)
    except:
        litcovdf = pd.DataFrame([])
        
    subtabs = ['Protein expression', 'Protein-TF correlation', 'Drugs', f'Literature references ({len(litcovdf)})', 'PPI', 'Download']

    option_data2 = [{'icon': icon, 'label':label} for 
            label, icon in zip(
                subtabs, 
                ['📈', '🔗', '💊', '📚', ' 🔄','📩']
            )
    ]
    
    chosen_subtab = hc.option_bar(
        option_definition=option_data2,
        title='',
        key='Sub2',
        override_theme=protein_bar,
        horizontal_orientation=True)
        
        
    if chosen_subtab == 'Protein expression':
        display_protein_expressions(selected_protein, datacore.name2ds[p_ds], cell_type)
        
    elif chosen_subtab == 'Protein-TF correlation':
        display_protein_tf_corr_data(selected_protein, cell_type, datacore.name2ds[p_ds])
        
    elif chosen_subtab == 'Drugs':
        drugs_df = pd.read_csv('./data/drugs/repurposing_drugs_20180907.txt', sep='\t')
        drugs_df = drugs_df.dropna(subset=['target'], axis=0)
        st.caption('Examples: CD44, CD38...')
        # for p in s_proteins:
        #     _drugs_df = drugs_df[drugs_df.target.apply(lambda x: p in x.split('|'))]
        #     if len(_drugs_df) > 0:
        #         st.title(p)
        #         # st.table(_drugs_df)
        #         st.write(_drugs_df.to_dict())
        
                    
        _drugs_df = drugs_df[drugs_df.target.apply(lambda x: selected_protein in x.split('|'))]
        if len(_drugs_df) > 0:
            # st.title(selected_protein)
            st.table(_drugs_df.drop('target', axis=1))     
            
        else:
            st.warning(f'No drugs found for {selected_protein}.')       

            
    elif 'Literature references' in chosen_subtab:
        display_litcovid_data(litcovdf, selected_protein)
    
    elif chosen_subtab == 'PPI':
        display_ppi_data(selected_protein)
        
        
    elif chosen_subtab == 'Download':
        info = st.empty()
        info.info('⌛ Generating files')
        pbar = st.progress(10)
        
        tf_activities = display_protein_expressions(
            selected_protein, datacore.name2ds[p_ds], cell_type, disabled=True) 
        pbar.progress(25)
        pbar.progress(30)
        
        csv_path, img_path = display_protein_tf_corr_data(
            selected_protein, cell_type, datacore.name2ds[p_ds], disabled = True)
        pbar.progress(40)
        
        html_ppi, wiki_csv =  display_ppi_data(selected_protein, disabled = True)
        
        if tf_activities is not None:
            buf = BytesIO()
            tf_activities.save(buf, format='PNG')
            png = buf.getvalue()
            
            pbar.progress(45)
            
            download_button(
                    desc = f'{selected_protein} expression density plots\n',
                    label = 'Download as PNG', 
                    data = png, 
                    fname = f'{selected_protein}_tf_activities.png',
                    format="application/octet-stream"
            )

        pbar.progress(55)
        
        pbar.progress(65)
                
        download_button(
                desc = 'Literature references\n', 
                label = 'Download as CSV', 
                data = litcovdf.to_csv(), 
                fname = f'{selected_protein}_litcovid_references.csv'
        )
        
        pbar.progress(75)
        if img_path:
            download_button(
                    desc = f'{selected_protein} protein-tf correlation heatmap\n', 
                    label = 'Download as PNG', 
                    data = img_path2buffer(img_path), 
                    fname = f'{selected_protein}_protein-TF_corr.png'
            )
            
            
            pbar.progress(85)
            
            download_button(
                    desc = f'{selected_protein} tf-protein correlation csv\n', 
                    label = 'Download as CSV', 
                    data = pd.read_csv(csv_path).to_csv(), 
                    fname = f'{selected_protein}_TF-protein_corr.csv'
        )
        
        if html_ppi:
            download_button(
                    desc = f'{selected_protein} StringDB network\n', 
                    label = 'Download as PNG', 
                    data = html_ppi, 
                    fname = f'{selected_protein}_ppi.html'
            )
            
        

        pbar.progress(100)
        info.success('Downloads ready!')
        
        

            
            
        

