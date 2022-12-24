import streamlit as st

from utils.api import get_pubmed_article
from utils.components import local_css, st_url, totags


def display_litcovid_data(litcovdf, option):
    
    show_abstracts = st.checkbox('Show abstracts')
    st.markdown(f'### 📄 Relevant publications for {option} based on LitCovid')
    
    
    if len(litcovdf) > 0: 
        for title in litcovdf.title:
            article = get_pubmed_article(title)
            if article:
                pmid = article.pubmed_id.split()[0]
                st_url(text=article.title, link=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}", tags='#####')
                if show_abstracts:
                    st.write(article.abstract)
                
                local_css()
                author = article.authors[0]
                author = f"{author['lastname']} et al. ({article.publication_date.year})"
                
                st.markdown(totags([author, article.journal, "PMID:"+pmid]), unsafe_allow_html=True)
                st.markdown('---')
    else:
        st.warning('No literature references found')
        
    