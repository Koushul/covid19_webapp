import streamlit as st

def multi_dataset_page(datacore):
    st.info("""
        Discover transcription factors that are consistently and significantly different between healthy vs impacted patients
            
    """)
    
    analysis_type = st.selectbox('Search for', ['Transcription Factors', 'Surface Proteins'])
                                 
    thresh = st.slider('p-value threshold', min_value=0.0, max_value=0.5, value=0.05, step=0.01)
    
    names = list(datacore.name2ds.keys())
    
    
    datasets = st.multiselect('Datasets', names, [names[0], names[1]])
    
    
    if len(datasets) < 1:
        # run = st.button('Search')
        st.error('Select at least one datasets')
        
    else:
        common_celltypes = set.intersection(*[set(datacore.name2ds[ds].celltypes) for ds in datasets])    
        results = {} 
        
        my_bar = st.progress(0)
        
        
        
        with st.spinner(f'👀 Looking across {len(datasets)} datasets... (p-value < {thresh})'):
            for idx, ct in enumerate(common_celltypes):
                results[ct] = {}
                for ds in datasets:
                    results[ct][ds] = {}
                    if analysis_type == 'Transcription Factors':
                        candidates = list(datacore.name2ds[ds].transcription_factors(ct))
                    else:
                        candidates = list(datacore.name2ds[ds].surface_proteins(ct))
                    for tf in candidates:
                        # try:
                        if analysis_type == 'Transcription Factors':
                            pval, adj_pval = datacore.name2ds[ds].get_tf_pvals(ct, tf)
                        else:
                            pval, adj_pval = datacore.name2ds[ds].get_protein_pvals(ct, tf)
                        
                        if pval < thresh:
                            results[ct][ds][tf] = pval
                        
                my_bar.progress((idx)/(len(common_celltypes)-1))
                        
        for ct in common_celltypes:
            significant = set.intersection(*[set(results[ct][ds].keys()) for ds in datasets])
            if len(significant)>0:
                st.markdown('## '+ct)
                st.write(significant)
                
                cols = st.columns(len(datasets))
                
                for stf in significant:
                    for col, ds in zip(cols, datasets):
                        col.metric(label=f'{stf} | {ds}', value = round(results[ct][ds][stf], 4), help='p-value computed using wilcox test')
                        
                st.markdown('---')
                
        # st.write(results)
