from navigation.datasets import dataset_page
from navigation.doc import doc_page
from navigation.home import home_page
from navigation.multi_dataset_search import multi_dataset_page
from navigation.protein import protein_page
from navigation.stats import stats_page
from navigation.tf import tf_page
from navigation.viral_entry import viral_entry_page
from utils.dataloaders import Dataset
from utils.components import footer_style
from utils.dataloaders import make_datacore

from PIL import Image
import streamlit as st

import hydralit_components as hc


st.set_page_config(
        page_title='COVID-19db of Immune Cell States',
        page_icon="./assets/logos/spartan.png",
        # layout="wide",
        initial_sidebar_state="expanded",
)



Image.MAX_IMAGE_PIXELS = None

HOME = "Home"
DATASET = "Datasets"
TF = "Transcription factors"
PROTEINS = "Surface proteins"
VIRAL_ENTRY = "Viral-entry related factors"
MULTI = "Multi datasets"
DOC = "Docs"
STATS = "Usage"


datacore =  make_datacore()()
ds = datacore.metadata




hide_menu_style = """
        <style>
        #MainMenu {visibility: hidden;}
        </style>
        """
# st.markdown(hide_menu_style, unsafe_allow_html=True)
st.markdown(footer_style, unsafe_allow_html=True) ## Footer

tabs = [
    HOME,
    DATASET,
    PROTEINS,
    TF,
    VIRAL_ENTRY,
    DOC, 
    STATS
]

option_data = [
   {'icon': "🏠", 'label':HOME},
   {'icon':"💾",'label':DATASET},
   {'icon': "🔎", 'label':PROTEINS},
   {'icon': "🔎", 'label':TF},
   {'icon': "🔍", 'label':VIRAL_ENTRY},
   {'icon': "💼 💼", 'label':"Multi datasets"},
   {'icon': "📋", 'label':DOC},
   {'icon': "📈", 'label':STATS},
]

over_theme = {'txc_inactive': 'black','menu_background':'white','txc_active':'white','option_active':'red'}
font_fmt = {'font-class':'h3','font-size':'50%'}

chosen_tab = hc.option_bar(
    option_definition=option_data,
    title='',
    key='PrimaryOptionx',
    override_theme=over_theme,
    horizontal_orientation=True)

st.markdown('----')
    
#######--------------------#######
#######--------------------#######

if chosen_tab == HOME:
    home_page()
    
elif chosen_tab == DATASET:
    dataset_page(ds, datacore)
    
elif chosen_tab == TF:
    tf_page(datacore)
    
elif chosen_tab == PROTEINS:
    protein_page(datacore)
        
elif chosen_tab == VIRAL_ENTRY:
    viral_entry_page(datacore)
        
elif chosen_tab == MULTI:
    multi_dataset_page(datacore)
          
elif chosen_tab == DOC:
    doc_page()
    
elif chosen_tab == STATS:
    stats_page()
        
