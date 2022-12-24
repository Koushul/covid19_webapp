import streamlit as st
from PIL import Image
import datetime

import streamlit_analytics


@st.cache(allow_output_mutation=True)
def load_home_image():
    return Image.open('./assets/images/home.png')

TITLE = 'COVID-19db linkage maps of cell surface proteins and transcription factors in immune cells'



def home_page():

    

    from texts.descriptions import Desc
    img = load_home_image()
    
    
    streamlit_analytics.start_tracking()    
    st.markdown('### ' + TITLE) 
    streamlit_analytics.stop_tracking()
    
    with open('./clock.time', 'r') as f:
        last_updated_on = f.readlines()[0]
    
    st.caption(last_updated_on)    
    st.image(img)      
    st.write(Desc.home)
    
        

    
    
    
if __name__ == '__main__':
    t = datetime.datetime.now()
    with open('./clock.time', 'w') as f:
        f.write(f'Last updated at {t.time():%H:%M} on {t.date():%Y-%m-%d}\n')
        print(f'Last updated at {t.time():%H:%M} on {t.date():%Y-%m-%d}\n')
    
    
    