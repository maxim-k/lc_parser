from parser.chromatogram import Chromatogram
from pathlib import Path

import pandas as pd
import numpy as np


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
