import streamlit as st
import requests
from bs4 import BeautifulSoup
import pandas as pd
import re
from pymed import PubMed

pubmed = PubMed(tool="MyTool", email="my@email.address")

def severus_pair(tf, protein):
    try:
        page = requests.get(f'http://severus.dbmi.pitt.edu/corona/index.php/search?q={tf}+and+{protein}')
        soup = BeautifulSoup(page.content, "html.parser")
        url = soup.findAll('a', 
                attrs={'href': 
                    re.compile("^http://severus.dbmi.pitt.edu/corona/index.php/search\?q=")})[0].get('href').replace(' ', '%20')
    
        return pd.read_html(url)[0]
    except:
        return pd.DataFrame([])


@st.cache
def load_tftargets(tf: str) -> pd.DataFrame:
    try:
        return pd.read_csv(f'./data/hTFtarget/{tf}.target.txt', sep='\t')
    except:
        return pd.DataFrame([])


@st.cache
def load_drug_info():
    drug_info = pd.read_parquet('./data/drugs/drugs_info.parquet')
    return drug_info

@st.experimental_singleton
def get_pubmed_article(title):
    articles = list(pubmed.query(title, max_results=1))
    if articles:
        article = articles[0]
        return article
    else:
        return None


@st.cache
def litcov_search(search):
    lines = []
    try:
        with open(f'./litcovid/tsv?text={search}', 'r') as f:
            lines = f.readlines()[33:]
        if lines:
            litcovdf = pd.DataFrame([x.split('\t') for x in lines[1:]], 
                    columns=['pmid', 'title', 'journal']).applymap(lambda x: x.replace('\n', ''))
            return litcovdf
        else:
            return pd.DataFrame([])
    
    except FileNotFoundError:
        return pd.DataFrame([])
