import streamlit as st
from navigation.tf_pages.activities import display_tf_activities
from navigation.tf_pages.litcovid import display_litcovid_data
from navigation.tf_pages.mrna import display_tf_target_data
from navigation.tf_pages.tf_protein_corr import display_tf_protein_corr_data
from texts.descriptions import Desc
from utils.api import litcov_search, load_drug_info

from utils.components import display_ppi_data
import hydralit_components as hc
import pandas as pd
from io import BytesIO
from PIL import Image

tf_bar = {'txc_inactive': 'black','menu_background':'white','txc_active':'white','option_active':'blueviolet'}

drugs = load_drug_info()


def download_button(desc: str, data, fname: str, label: str = 'Download', format: str ='text/csv'):
    left_col, right_col = st.columns([2, 1])
    right_col.download_button(
        label=label,
        data=data,
        file_name=fname,
        mime=format,
    )
        
    left_col.write(desc)

@st.cache   
def img_path2buffer(img_path: str):
    img = Image.open(img_path)
    buf = BytesIO()
    img.save(buf, format='PNG')
    png = buf.getvalue()
    return png

    

    
def tf_page(datacore):
    st.info(Desc.tf)
    _, _, edge = st.columns(3)
    checkbox = edge.checkbox('Show high correlations only')
    aa, ba, ca = st.columns(3)
    st.session_state['datasets'] = ["GSE155673"]
    
    datasets = datacore.name2ds.keys()
    tf_ds = aa.selectbox(f'Dataset ({len(datasets)})', datasets, 1)
    cell_types = datacore.name2ds[tf_ds].celltypes
    cell_type = ba.selectbox(f'Cell Type ({len(cell_types)})', cell_types, 0)

    if not checkbox:
        tfs = datacore.name2ds[tf_ds].total_transcription_factors            
    else:
        tfs = list(datacore.name2ds[tf_ds].transcription_factors(cell_type))
    
    selected_tf = ca.selectbox(
        f'Transcription factor of interest ({len(tfs)})', tfs)   
    
    litcovdf = litcov_search(selected_tf)
    targets = drugs[drugs.Target.apply(lambda x: selected_tf in x)]

    subtabs = [
        'TF activities', 
        'TF-Target gene', 
        'TF-protein correlation', 
        f'Drugs',
        f'Literature references ({len(litcovdf)})', 
        'PPI', 
        'Download'
    ]

    option_data2 = [
        {'icon': icon, 'label':label}
            for label, icon in zip(
                            subtabs, 
                            ['📈', '🎯', '🔗', '💊', '📚','🔄', '📩']
                        )
    ]
    
    chosen_subtab = hc.option_bar(
        option_definition=option_data2,
        title='',
        key='Sub1',
        override_theme=tf_bar,
        horizontal_orientation=True)
    
    
    
    
    
    if chosen_subtab == 'TF activities':  
        display_tf_activities(selected_tf, datacore.name2ds[tf_ds], cell_type)   
        
    elif chosen_subtab == 'TF-Target gene':
        display_tf_target_data(selected_tf, cell_type, ds=datacore.name2ds[tf_ds])
    
    elif chosen_subtab == 'TF-protein correlation':
        display_tf_protein_corr_data(selected_tf, cell_type, datacore.name2ds[tf_ds])
            
    elif 'Drugs' in chosen_subtab:    
        
        # for ds in datacore.datasets:
        #     for ct in ds.celltypes:
        #         for x in ds.transcription_factors(ct):
        #             pval, _ = ds.get_tf_pvals(ct, x)
        #             if pval < 0.05:
        #                 t = drugs[drugs.Target.apply(lambda z: x in z)]
        #                 if len(t) > 0 and len(litcovdf) > 2:
        #                     st.write(f'{x} | {pval} | {ds.name} | {ct}')
                        
                        
            
        if len(targets) > 0:
            st.table(targets.drop(['master regulators', 'nTarget', 'ChEMBL'], axis=1))
        else:
            st.warning(f'No drugs found for {selected_tf}.')
        
    elif 'Literature references' in chosen_subtab:
        display_litcovid_data(litcovdf, selected_tf)
        
    elif chosen_subtab == 'PPI':
        display_ppi_data(selected_tf)
        
    elif chosen_subtab == 'Download':
        info = st.empty()
        info.info('⌛ Generating files')
        
        pbar = st.progress(10)
        
        tf_activities = display_tf_activities(
            selected_tf, datacore.name2ds[tf_ds], cell_type, disabled=True) 
        pbar.progress(25)
        
        tf_target_img, tf_target_df = display_tf_target_data(
            selected_tf, cell_type, ds=datacore.name2ds[tf_ds], disabled = True)
        pbar.progress(30)
        
        csv_path, img_path = display_tf_protein_corr_data(
            selected_tf, cell_type, datacore.name2ds[tf_ds], disabled = True)
        pbar.progress(40)
        
        html_ppi, wiki_csv =  display_ppi_data(selected_tf, disabled = True)
        
        if tf_activities is not None:
            buf = BytesIO()
            tf_activities.save(buf, format='PNG')
            png = buf.getvalue()
            
            pbar.progress(45)
            
            download_button(
                    desc = 'TF activities density plots\n',
                    label = 'Download as PNG', 
                    data = png, 
                    fname = f'{selected_tf}_tf_activities.png',
                    format="application/octet-stream"
            )
        
        download_button(
                desc = 'TF-Target gene heatmap\n', 
                label = 'Download as PNG', 
                data = img_path2buffer(tf_target_img), 
                fname = tf_target_img
        )
        
        pbar.progress(55)
        
        
        download_button(
                desc = f'{selected_tf} target genes from hTFtarget\n', 
                label = 'Download as CSV', 
                data = tf_target_df.to_csv(), 
                fname = f'{selected_tf}_htftargets.csv'
        )
        
        pbar.progress(65)
                
        download_button(
                desc = 'Literature references\n', 
                label = 'Download as CSV', 
                data = litcovdf.to_csv(), 
                fname = f'{selected_tf}_litcovid_references.csv'
        )
        
        pbar.progress(75)
        
        download_button(
                desc = f'{selected_tf} tf-protein correlation heatmap\n', 
                label = 'Download as PNG', 
                data = img_path2buffer(img_path), 
                fname = f'{selected_tf}_TF-protein_corr.png'
        )
        
        
        pbar.progress(85)
        
        download_button(
                desc = f'{selected_tf} tf-protein correlation csv\n', 
                label = 'Download as CSV', 
                data = pd.read_csv(csv_path).to_csv(), 
                fname = f'{selected_tf}_TF-protein_corr.csv'
        )
        
        if html_ppi:
            download_button(
                    desc = f'{selected_tf} StringDB network\n', 
                    label = 'Download as PNG', 
                    data = html_ppi, 
                    fname = f'{selected_tf}_ppi.html'
            )
            
        

        pbar.progress(100)
        info.success('Downloads ready!')
        
    
