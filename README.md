# Liquid Chromatography Data Viewer

## Overview
This application processes liquid chromatography data, specifically from the [ACQUITY UPLC FLR machine](https://www.waters.com/webassets/cms/support/docs/71500142403rb.pdf). It automates the parsing of metadata, detects peaks, interpolates the baseline, and calculates peak areas and elution volumes.

![Screenshot](screenshot.png)

## Features
### Baseline Interpolation and Detrending
Interpolates the baseline to detect peaks more accurately. The process uses a Savitzky-Golay filter to smooth the chromatogram data, helping with true baseline identification. This method reduces noise and improves the signal-to-noise ratio for better peak detection. Detrending involves subtracting the interpolated baseline from the original data to normalize the data around zero and make peaks more distinguishable.
#### Parameters
`poly_degree`: Polynomial degree for the interpolation of the baseline. Default value is `2`

### Peak Detection
Automatically detects peaks in chromatography data. Peaks are detected using the `find_peaks` function from `scipy.signal`, which identifies local maxima in the data.
#### Parameters:
`min_height`: Minimum height required for a peak to be considered. Default value is one standard deviation above the mean.
`prominence`: Vertical distance from the lowest contour line to the peak, indicating how much it stands out from the signal's baseline. Default value is `1`.
`peak_window_length`: Window length in samples that limits the evaluated area for each peak. Default value is `100`.
`sg_window_length`: Length of the window for the Savitzky-Golay filter. Default value is `25`.

### Peak Area Calculation
Calculates the area under the peak using Simpson's rule.

### Elution Volume Calculation
Computes the elution volume based on the flow rate and peak retention time.

## Usage 

### Prerequisites
- Place HPLC data files in `~/lc_parser/data/`.

### Run locally as a Streamlit app
- Use `streamlit run streamlit_app.py`
- To run in Docker use `docker build -t lc_parser:latest .` and `docker run -p 8501:8501 lc_parser:latest`

### Use the package in Python script

```bash
pip install -r requirements.txt
```

```python
from pathlib import Path
from parser.chromatogram import Chromatogram

# Load chromatogram data
file_path = Path("./data/sample_data.txt")
chrom = Chromatogram(file_path)

# Process data
chrom.detect_peaks()
for peak in chrom.peaks:
    area = chrom.calculate_peak_area(peak)
    volume = chrom.calculate_elution_volume(peak, flow_rate=1.0)
    print(f"Peak at {peak.retention_time} min: Area = {area}, Volume = {volume}")

```
Output:
```
Peak at 2.34 min: Area = 150.5, Volume = 2.34
Peak at 5.22 min: Area = 300.1, Volume = 5.22
```

### Use in a Jupyter Notebook
See [the demo notebook](https://github.com/maxim-k/lc_parser/blob/main/demo.ipynb) in this repo.


## Running tests
To run the tests, execute the following command:
```bash
pytest tests
```