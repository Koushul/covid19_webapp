import pandas as pd
import streamlit as st
from navigation.home import streamlit_analytics
import plotly.graph_objects as go


def stats_page():
    # cities = pd.read_table('https://gist.githubusercontent.com/weshi/03c92ae63ad5554baff487dfd866b0b5/raw/fd875dea90d3d3407d95e11bac19fc3b4d1f6bfc/worldcities.csv', sep=',')
    # cities = cities[['Latitude', 'Longitude']]
    # cities.columns = ['lat', 'lon']
    # st.map(cities, zoom=2)
    view_data = streamlit_analytics.main.counts['per_day']
    views = view_data['pageviews']
    start_time = streamlit_analytics.main.counts["start_time"]
        
        
    
    _, y, z, w = st.columns([2, 8, 1, 3])
    
    y.header(f'👤 COVID-19db daily users')
    w.metric('Pageviews', 
        views[-1],
        help = f"Number of times the app home page was visited since {start_time}")
        
           
    views_df = pd.DataFrame([view_data['days'], view_data['pageviews']]).T
    views_df.columns = ['date', 'pageviews']
            
    fig = go.Figure([go.Scatter(x=views_df['date'], y=views_df['pageviews'])])
    fig.update_xaxes(rangeslider_visible=True)
    
    st.plotly_chart(fig, use_container_width=True, theme='streamlit')

        
        
        
    
    
    
