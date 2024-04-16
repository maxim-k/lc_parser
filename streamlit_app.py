import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path
import os

from parser.chromatogram import Chromatogram

st.set_page_config(layout="wide")


def plot_processed_data(data, peaks):
    fig = px.line(data, x='Time (min)', y='Value (EU)', title='Raw Data')
    # Plotting peaks with the data from Peak instances
    for peak in peaks:
        peak_data = peak.data
        fig.add_scatter(x=peak_data['Time (min)'], y=peak_data['Value (EU)'], mode='markers', name=f'Peak at {peak.retention_time} min')
    return fig


# Streamlit app layout
st.title('Liquid Chromatography Data Viewer')
st.write('##### This app visualizes HPLC data and identified peaks.')

st.divider()

data_dir = 'data'
files = [f for f in os.listdir(data_dir) if f.endswith('.txt')]

settings, plot = st.columns([1, 2])

with settings:
    # File selector
    file_selector = st.selectbox('Select a data file:', files)
    # Sliders for peak detection parameters
    poly_degree = st.slider("Polynomial degree for the baseline correction", min_value=1, max_value=10, value=3, step=1,
                            help="Determines the degree of the polynomial used for the interpolation for the baseline correction.")
    min_height = st.slider("Minimum height of peaks", min_value=1.0, max_value=25.0, value=0.0, step=0.1,
                           help="Minimum height required to recognize a peak.")
    prominence = st.slider("Prominence of peaks", min_value=2.0, max_value=25.0, value=0.1, step=0.1,
                           help="Required prominence of peaks. This parameter helps to distinguish the main peak from the surrounding noise.")
    peak_window_length = st.slider("Window size for peaks", min_value=1, max_value=1000, value=100, step=1,
                              help="Window length in samples that limits the evaluated area for each peak.")
    sg_window_length = st.slider("Window size for Savitzky-Golay filter", min_value=5, max_value=101, value=25, step=2,
                              help="The length of the filter window. Must be a positive odd number.")

with plot:
    # File processing
    if file_selector:
        data_path = Path(data_dir) / file_selector
        chromatogram = Chromatogram(data_path)
        chromatogram.detect_peaks(poly_degree=poly_degree, min_height=min_height, prominence=prominence,
                                  peak_window_length=peak_window_length,
                                  sg_window_length=sg_window_length)

        st.write('### Peaks')
        processed_fig = plot_processed_data(chromatogram.raw_data, chromatogram.peaks)
        st.plotly_chart(processed_fig, use_container_width=True)
