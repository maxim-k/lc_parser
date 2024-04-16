from io import StringIO
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from parser.peak import Peak
from scipy.integrate import simpson
from scipy.signal import find_peaks, savgol_filter


class Chromatogram:
    """
    Handles the processing of chromatogram data from a file.

    :param filepath: The path to the chromatogram data file.
    :raises ValueError: If the filepath is not a pathlib.Path object.
    :raises FileNotFoundError: If no file is found at the provided filepath.
    """

    def __init__(self, filepath: Path):
        if not isinstance(filepath, Path):
            raise ValueError("filepath must be a pathlib.Path object")
        if not filepath.is_file():
            raise FileNotFoundError(f"No file found at {filepath}")

        self.metadata: Dict[str, Dict[str, str]] = {
            "Injection Information": {},
            "Chromatogram Data Information": {},
            "Signal Parameter Information": {},
            "Other": {},
        }
        self.raw_data: pd.DataFrame = pd.DataFrame()
        self.peaks: list[Peak] = []
        self._parse_file(filepath)

    def _parse_file(self, filepath: Path):
        """
        Parses the file at the given path to extract chromatogram data and metadata.

        This method reads a text file containing metadata and raw chromatogram data,
        separates the content into metadata and data sections, parses the metadata
        for different categories, and converts the raw data section into a pandas DataFrame.

        :param filepath: Path object pointing to the file to be parsed.
        :type filepath: Path

        :raises FileNotFoundError: If the file cannot be found or read.
        :raises ValueError: If the file content does not contain expected sections or format.
        :raises KeyError: If essential data columns are missing in the chromatogram data section.
        """
        try:
            with filepath.open("r") as file:
                content = file.read().strip()
        except IOError as e:
            raise FileNotFoundError(f"Error reading file: {e}")

        try:
            metadata, data = content.split("Chromatogram Data:\n")
        except ValueError:
            raise ValueError(
                "File format incorrect or 'Chromatogram Data:' section missing"
            )

        current_section = "Other"  # TODO What if it is in the middle

        for line in metadata.split("\n"):
            if line.endswith(":"):
                current_section = line[:-1]
            elif line:
                key, value = line.split("\t")
                self.metadata[current_section][key] = value

        self.raw_data = pd.read_csv(
            StringIO(data), sep="\t", na_values="n.a."
        )  # TODO check other NaN

        try:
            self.raw_data["Time (min)"].to_numpy()
            self.raw_data["Value (EU)"].to_numpy()
        except KeyError as e:
            raise ValueError(f"Expected column missing from the data: {e}")

    def detect_peaks(
        self,
        poly_degree: int = 3,
        min_height: float = None,
        prominence: float = 1,
        peak_window_length: int = 100,
        sg_window_length: int = 25,
    ):
        """
        Detects peaks in the chromatogram data.

        Applies polynomial interpolation for detrending, and a Savitzky-Golay filter to smooth the data. Uses the find_peaks method to identify peaks based on height and other characteristics.

        :param poly_degree: Polynomial degree for the interpolation of the baseline.
        :param min_height: Minimum height required for a peak to be considered.
        :param prominence: Vertical distance from the lowest contour line to the peak, indicating how much it stands out from the signal's baseline.
        :param peak_window_length: Window length in samples that limits the evaluated area for each peak.
        :param sg_window_length: Length of the window for the Savitzky-Golay filter.
        :raises ValueError: If the window size is too large for the dataset or the dataset is too short for effective smoothing.
        """
        values = self.raw_data["Value (EU)"].to_numpy()
        time = self.raw_data["Time (min)"].to_numpy()
        if len(values) < sg_window_length:
            raise ValueError(
                "Window size for Savitzky-Golay filter must be less than the length of the dataset"
            )
        if sg_window_length % 2 == 0:
            raise ValueError(
                "Window size for Savitzky-Golay filter must be must be an odd number"
            )

        # Polynomial fit for detrending
        coeffs = np.polyfit(time, values, deg=poly_degree)
        baseline_func = np.poly1d(coeffs)
        baseline = baseline_func(time)

        # Detrending the data by subtracting the fitted polynomial baseline
        detrended_data = values - baseline

        # Apply Savitzky-Golay filter to smooth the detrended data
        values = savgol_filter(detrended_data, sg_window_length, polyorder=3)

        # If not specified set min height to be more than one standard deviation above the mean
        if min_height is None:
            mean_emission = values.mean()
            std_emission = values.std()
            min_height = mean_emission + std_emission

        peaks, properties = find_peaks(
            values, height=min_height, prominence=prominence, wlen=peak_window_length
        )

        for i, peak in enumerate(peaks):
            peak_height = properties["peak_heights"][i]
            retention_time = time[peak]

            # Find the left and right base indices and slice the data for the peak
            left_base_idx = properties["left_bases"][i]
            right_base_idx = properties["right_bases"][i]
            peak_data = self.raw_data.iloc[left_base_idx : right_base_idx + 1]

            self.peaks.append(
                Peak(
                    left_base_idx,
                    right_base_idx,
                    peak_height,
                    retention_time,
                    peak_data,
                )
            )

    def calculate_peak_area(self, peak: Peak) -> float:
        """
        Calculates the area under the curve for a given peak using Simpson's rule.

        :param peak: The Peak object for which to calculate the area.
        :return: The calculated area under the curve.
        :raises ValueError: If the peak data is empty.
        """
        if peak.data.empty:
            raise ValueError("Peak data is empty, cannot calculate area.")
        area = simpson(
            peak.data["Value (EU)"], peak.data["Time (min)"]
        )
        return area

    def calculate_elution_volume(self, peak: Peak, flow_rate: float) -> float:
        """
        Calculates the elution volume for a given peak based on its retention time and the specified flow rate.

        :param peak: The Peak object for which to calculate the elution volume.
        :param flow_rate: The flow rate used in the experiment, must be a positive number.
        :return: The calculated elution volume.
        :raises ValueError: If the flow rate is not a positive number.
        """
        if flow_rate <= 0:
            raise ValueError("Flow rate must be a positive number.")
        elution_volume = peak.retention_time * flow_rate
        return elution_volume


if __name__ == "__main__":
    filepath = (
        Path(__file__).parent.parent / "data" / "IgG Vtag 1_ACQUITY FLR ChA.txt"
    )
    chrom = Chromatogram(filepath)
    chrom.detect_peaks()
    peaks = chrom.peaks
    print()
