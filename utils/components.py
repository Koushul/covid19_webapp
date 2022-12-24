from io import TextIOWrapper
import re
from typing import Tuple
import uuid
import numpy as np
import stringdb
from PIL import Image
import extra_streamlit_components as stx
import streamlit as st
from pyvis.network import Network
import streamlit.components.v1 as components
from matplotlib import colors
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


footer_style = f"""
    <style>
        MainMenu {{visibility: hidden;}}
        footer {{visibility: hidden;}}
        footer:after {{
            content:'Created by Koushul using Streamlit'; 
            visibility: visible;
            display: block;
            position: relative;
            #background-color: red;
            padding: 5px;
            top: 2px;
        }}
    </style>
"""



def get_ppi_edge_list(gene, neighbors, size):
    genes = [gene]
    string_ids = stringdb.get_string_ids(genes)
    enrichment_df = stringdb.get_interaction_partners(string_ids.queryItem, limit=size)
    extras = list(enrichment_df['preferredName_B'].values[:neighbors])
    string_ids = stringdb.get_string_ids(extras+genes)
    enrichment_df = stringdb.get_interaction_partners(string_ids.queryItem, limit=size)

    return string_ids, enrichment_df[['preferredName_A',	'preferredName_B', 'score']].values


def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, color='black', fontsize=35,
            ha="left", va="center", transform=ax.transAxes)
    ax.tick_params(axis='x', labelsize=30)
    
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})

## https://matplotlib.org/stable/gallery/color/named_colors.html
palette_mapper = {
    'healthy': colors.to_hex('seagreen'),
    'mild': colors.to_hex('salmon'),
    'moderate': colors.to_hex('lightskyblue'),
    'critical': colors.to_hex('red'),
    'severe': colors.to_hex('magenta'),
    'stable': colors.to_hex('orange'),
    'progressive': colors.to_hex('blueviolet'),
}

severity_order = [
    'healthy',
    'mild',
    'moderate',
    'severe',
    'critical'
]

def st_url(text, link, tags='######'):
    st.write(tags+f" [{text}]({link})")

def totags(tags):
    t = '<div>'
    for tag, color in zip(tags, ['red', 'blue', 'green', 'orange']):
        t += f"<span class='highlight {color}'>{tag}</span>\t"
    t += '</div>'
    return t

def make_palette(column, groups):
    return [colors.to_rgb(column.color_picker(g, palette_mapper.get(g.lower(), colors.to_hex('black')))) for g in groups]    


def fig2img(fig):
    """Convert a Matplotlib figure to a PIL Image and return it"""
    import io
    buf = io.BytesIO()
    fig.savefig(buf)
    buf.seek(0)
    img = Image.open(buf)
    return img

def ridge_plot(df_plt: pd.DataFrame, xlabel: str, column: st.container, group: str='patient_group'):
    groups_ = np.sort(np.unique(df_plt[group].values))
    groups_ = np.array([i.lower() for i in groups_])
    
    # ## dirty sorting
    # groups = []
    # for g in severity_order:
    #     if g in groups_:
    #         groups.append(g)
    
    groups = groups_
                
    if column is not None:
        palette = make_palette(column=column, groups=groups)
    else:
        palette = sns.color_palette('pastel')
        
    g = sns.FacetGrid(df_plt, 
        palette = palette,
        row=group, hue=group, aspect=5)
        
    g.map_dataframe(sns.kdeplot, x="Pair_correlations", fill=True, alpha=1)
    g.map_dataframe(sns.kdeplot, x="Pair_correlations", color='black')
    g.map(label, group)
    g.fig.subplots_adjust(hspace=-.25)
    g.set_titles("")
    g.set(yticks=[], xlabel=xlabel)
    g.despine( left=True)
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    plt.xlabel(xlabel, fontsize=36, labelpad=25)
    plt.gcf().subplots_adjust(bottom=0.25)
    # plt.tight_layout()
    img = fig2img(g.fig)
    
    return img

def create_st_button(link_text, link_url, hover_color="#4EC5F1", st_col=None):    
    button_uuid = uuid.uuid1()
    button_id = re.sub("\d+", "", str(button_uuid))
    
    button_css = f"""
        <style>
            #{button_id} {{
                background-color: none;
                color: rgb(38, 39, 48);
                padding: 0.25em 0.38em;
                position: relative;
                text-decoration: none;
                border-radius: 4px;
                border-width: 2px;
                border-style: solid;
                border-color: rgb(13, 242, 201);
                border-image: initial;
            }}
            #{button_id}:hover {{
                border-color: {hover_color};
                color: {hover_color};
            }}
            #{button_id}:active {{
                box-shadow: none;
                background-color: {hover_color};
                color: white;
                }}
        </style> """

    html_str = f'<a href="{link_url}" target="_blank" id="{button_id}";>{link_text}</a><br></br>'

    if st_col is None:
        st.markdown(button_css + html_str, unsafe_allow_html=True)
    else:
        st_col.markdown(button_css + html_str, unsafe_allow_html=True)


def make_tab(name, desc=''):
    return stx.TabBarItemData(id=name, title=name, description=desc)


def display_ppi_data(option: str, disabled: bool = False) -> Tuple[TextIOWrapper, pd.DataFrame]:
    if not disabled:
        st.markdown(f'#### Protein-Protein Interactions for {option} from StringDB')
        st.caption('Threshold of significance to include an interaction: 400')
        st.caption('Learn more: https://string-db.org/help/api/')
        left_col, right_col = st.columns(2)
        size = left_col.slider('Interactors', 2, 100, 10)
        neighbors = right_col.slider('Neighbors', 1, 10, 3)
    else:
        size = 10
        neighbors = 3
        

    net = Network()
        
    try:
        caption, pairs = get_ppi_edge_list(option, neighbors, size)
        preferred_name = caption[caption['queryItem']==option].preferredName.values[0]

        for i, j, v in pairs:
            net.add_node(j, label=j)
            
        for n in net.nodes:
            if n['label'] == preferred_name:
                n['color'] = 'red'
                 
        for pair in pairs:
            net.add_edge(pair[0], pair[1], value=pair[2]/10)
    except:
        if not disabled: st.warning('No PPI Information available from stringdb')
        return None, None
    # net.show_buttons(filter_=['physics'])
    
    net.show('./utils/ppi.html')
        
    HtmlFile = open("./utils/ppi.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    if not disabled: components.html(source_code, height = 650)
    
    
    if not disabled: 
        st.markdown(f'#### Pathways & drugs associated with {option} derived from Wiki-CORONA')
    
    try:
        df = pd.read_html(f'http://severus.dbmi.pitt.edu/corona/index.php/search?q={option}')[0]
        df.columns = ['Interactant Symbol', 'Name', 'Associated Pathways', 'Binding Drugs', 'Associated Diseases']
        if not disabled: st.dataframe(df)
    except:
        if not disabled: st.warning('No interactors found from Wiki-CORONA')
        df = None
        
        
    return HtmlFile, df
    
def download_button(desc, data, left_col, right_col, label='Download as CSV'):
    right_col.download_button(
        label=label,
        data=data,
        file_name='projD.csv',
        mime='text/csv',
        key=uuid.uuid1()
    )
        
    left_col.write(desc)
    
    
def local_css(file_name='./assets/style.css'):
    with open(file_name) as f:
        st.markdown('<style>{}</style>'.format(f.read()), unsafe_allow_html=True)
        
        
def resize(img: Image, factor: float = 2) -> Image:
    x, y = img.size
    x = int(x/factor)
    y = int(y/factor)
    img = img.resize((x, y))
    return img