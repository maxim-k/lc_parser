from parser.chromatogram import Chromatogram
from parser.peak import Peak
from pathlib import Path

import pandas as pd
import numpy as np


def test_chromatogram_init():
    # Assuming the file has been correctly placed and path is accurate for testing
    filepath = Path(__file__).parent.parent / "data" / "IgG Vtag 1_ACQUITY FLR ChA.txt"
    chrom = Chromatogram(filepath)
    assert isinstance(chrom, Chromatogram)
    assert bool(chrom.metadata['Injection Information'])
    assert not chrom.raw_data.empty


def test_parse_file():
    filepath = Path(__file__).parent.parent / "data" / "IgG Vtag 1_ACQUITY FLR ChA.txt"
    chrom = Chromatogram(filepath)

    # Metadata parsing
    assert "Injection Information" in chrom.metadata
    assert "Chromatogram Data Information" in chrom.metadata
    assert "Signal Parameter Information" in chrom.metadata

    assert chrom.metadata["Injection Information"]["Data Vault"] == "ChromeleonLocal"
    assert (
            chrom.metadata["Chromatogram Data Information"]["Time Min. (min)"] == "0.000000"
    )
    assert chrom.metadata["Signal Parameter Information"]["Signal Info"] == "techCHROM"

    # Raw data parsing
    assert not chrom.raw_data.empty
    assert isinstance(chrom.raw_data, pd.DataFrame)

    # Check first row
    assert chrom.raw_data.iloc[0]["Time (min)"] == 0.0
    assert np.isnan(chrom.raw_data.iloc[0]["Step (s)"])
    assert chrom.raw_data.iloc[0]["Value (EU)"] == 0.0


def test_detect_peaks():
    filepath = Path(__file__).parent.parent / "data" / "IgG Vtag 1_ACQUITY FLR ChA.txt"
    chrom = Chromatogram(filepath)
    chrom.detect_peaks()
    assert isinstance(chrom.peaks, list)
    assert all(isinstance(peak, Peak) for peak in chrom.peaks)
    assert len(chrom.peaks) > 0  # Expecting at least one peak detected


def test_get_peaks_df():
    filepath = Path(__file__).parent.parent / "data" / "IgG Vtag 1_ACQUITY FLR ChA.txt"
    chrom = Chromatogram(filepath)
    chrom.detect_peaks()
    assert isinstance(chrom.get_peaks_df(), pd.DataFrame)

def test_calculate_peak_area():
    # Create a simple peak with known area calculation
    data = pd.DataFrame({
        'Time (min)': [1, 2, 3],
        'Value (EU)': [0, 2, 0]
    })
    peak = Peak(left_base_idx=0, right_base_idx=2, height=2, retention_time=2, data=data)
    area = Chromatogram.calculate_peak_area(peak)
    assert round(area, 2) == 2.67  # Real area is 2, Simpson's area is 2.(6). It's a parabola - what can I do?


def test_calculate_elution_volume():
    data = pd.DataFrame({
        'Time (min)': [0, 1],
        'Value (EU)': [0, 1]
    })
    peak = Peak(left_base_idx=0, right_base_idx=1, height=1, retention_time=1, data=data)
    elution_volume = Chromatogram.calculate_elution_volume(peak, flow_rate=1.0)
    assert elution_volume == 1  # Elution volume = retention_time * flow_rate
