import streamlit as st

from st_aggrid import AgGrid, GridUpdateMode
from st_aggrid.grid_options_builder import GridOptionsBuilder

from utils.components import local_css
from texts.descriptions import DatasetDesc, Desc

from PIL import Image


@st.cache
def totags(tags):
    t = '<div>'
    for tag, color in zip(tags, ['red', 'blue', 'green', 'orange', 'purple', 'yellow']):
        t += f"<span class='highlight {color}'>{tag}</span>\t"
    t += '</div>'
    return t


from st_aggrid import JsCode


def dataset_page(ds, datacore):
        
    st.header(f'Dataset Browser')
    st.info(Desc.ds)
    st.caption('Select one or more datasets')
    ds2 = ds[['name', 'patients', 'cell_count', 'surface_proteins', 'transcription_factors', 'pmid']]
    
    # ds2.pmid = ds2.pmid.apply(lambda x: f'https://pubmed.ncbi.nlm.nih.gov/{x}/')
    
    gd = GridOptionsBuilder.from_dataframe(ds2)
    
    if 'datasets' in st.session_state:
        pre_selected_rows = [list(ds.name).index(i) for i in st.session_state['datasets']]
    else:
        pre_selected_rows = []
                
                
    gd.configure_selection(selection_mode='multiple', use_checkbox=True, pre_selected_rows=pre_selected_rows)
    
    
    cell_renderer =  JsCode("""function(params) {return `<a href='https://pubmed.ncbi.nlm.nih.gov/${params.value}/' target="_blank">${params.value}</a>`}""")
    

    gd.configure_column("pmid", cellRenderer=cell_renderer)
    
    grid_table = AgGrid(ds2, 
        gridOptions=gd.build(),
        update_mode=GridUpdateMode.SELECTION_CHANGED, 
        fit_columns_on_grid_load=True,
        height=175,
        width='100%',
        theme="streamlit",
        allow_unsafe_jscode=True,
    )

        
    selected_row = grid_table["selected_rows"]
    
    st.session_state['datasets'] = [x['name'] for x in selected_row]
    
    for idx, row_json in enumerate(selected_row):
        
        st.markdown('#### ' + row_json['name'] + ': ' + ds[ds.name==row_json['name']].title.values[0])
        st.markdown(DatasetDesc[row_json['name']])
        st.markdown(ds[ds.name==row_json['name']].url.values[0])
        
        exp = st.expander('View UMAP')

        try:
            umap = Image.open(datacore.name2ds[row_json['name']].umap)
            exp.image(umap.resize((600, 400)))
        except:
            exp.write('No UMAP available.')

        local_css()
        
        st.markdown(
            totags([
                row_json['cell_count'] + ' cells', 
                str(row_json['patients']) + ' patients',
                str(row_json['surface_proteins']) + ' surface proteins',
                str(row_json['transcription_factors']) + ' transcription factors',
                
                ],
            ), 
            unsafe_allow_html=True)
        
        st.markdown('---')