import re
import uuid

import extra_streamlit_components as stx
import streamlit as st
from navigation.protein import generate_color, get_ppi_edge_list
from pyvis.network import Network
import streamlit.components.v1 as components
from matplotlib import colors


hide_streamlit_style = f"""
    <style>
        MainMenu {{visibility: hidden;}}
        footer {{visibility: hidden;}}
        footer:after {{
            content:'Created by Koushul using Streamlit and BioRender'; 
            visibility: visible;
            display: block;
            position: relative;
            #background-color: red;
            padding: 5px;
            top: 2px;
        }}
    </style>
"""

from matplotlib import rcParams
import seaborn as sns
import matplotlib.pyplot as plt
from PIL import Image

# rcParams['figure.figsize'] = 11.7,8.27

def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, color='black', fontsize=15,
            ha="left", va="center", transform=ax.transAxes)
    
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})
# palette = sns.color_palette("Set2", 12)

def fig2img(fig):
    """Convert a Matplotlib figure to a PIL Image and return it"""
    import io
    buf = io.BytesIO()
    fig.savefig(buf)
    buf.seek(0)
    img = Image.open(buf)
    return img

def ridge_plot(df_plt, xlabel=''):
    green = colors.to_hex('green')
    salmon = colors.to_hex('salmon')
    lightskyblue = colors.to_hex('lightskyblue')
    orchid = colors.to_hex('orchid')

    g = sns.FacetGrid(df_plt, palette=
                    [
                        colors.to_rgb(st.color_picker('healthy', green)), 
                        colors.to_rgb(st.color_picker('mild', salmon)),
                        colors.to_rgb(st.color_picker('moderate', lightskyblue)),
                        colors.to_rgb(st.color_picker('severe', orchid))
                    ], 
        row="patient_group", hue="patient_group", aspect=5)
    
    g.map_dataframe(sns.kdeplot, x="Pair_correlations", fill=True, alpha=1)
    g.map_dataframe(sns.kdeplot, x="Pair_correlations", color='black')
    g.map(label, "patient_group")
    g.fig.subplots_adjust(hspace=-.5)
    g.set_titles("")
    g.set(yticks=[], xlabel=xlabel)
    g.despine( left=True)
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)
    # g.savefig('ridge.png', dpi=180)
    
    img = fig2img(g.fig)
    
    return img

def create_st_button(link_text, link_url, hover_color="#e78ac3", st_col=None):    
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
                border-color: rgb(230, 234, 241);
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


def ppi(option):
    st.header(f'Protein-Protein Interactions for {option} from stringdb')
    st.caption('Threshold of significance to include an interaction: 400')
    st.caption('Learn more: https://string-db.org/help/api/')
    left_col, right_col = st.columns(2)
    size = left_col.slider('Interactors', 2, 100, 10)
    neighbors = right_col.slider('Neighbors', 1, 10, 3)

    net = Network()
    try:
        caption, pairs = get_ppi_edge_list(option, neighbors, size)
        for i, j, v in pairs:
            net.add_node(j, label=j)
            
        for n in net.nodes:
            if n['label'] == option:
                n['color'] = 'red'
                 
        for pair in pairs:
            net.add_edge(pair[0], pair[1], value=pair[2]/10)
    except:
        st.warning('No PPI Information available from stringdb')
        return
    # net.show_buttons(filter_=['physics'])

    net.show('./utils/ppi.html')
    HtmlFile = open("./utils/ppi.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    components.html(source_code, height = 700, width=1000)
    
    # st.download_button(
    #     label='Download Graph',
    #     data=source_code,
    #     file_name = 'ppi.html'
    # )
    
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