from io import StringIO
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from peak import Peak
from scipy.integrate import simps
from scipy.signal import find_peaks


def calculate_elution_volume(peak: Peak, flow_rate: float):
    return peak.retention_time * flow_rate


class Chromatogram:
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
        It also initializes the baseline interpolation using the first and last data points.

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
            raise ValueError("File format incorrect or 'Chromatogram Data:' section missing")

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

    def detect_peaks(self, degree: int = 3, min_height: float = 0):
        values = self.raw_data["Value (EU)"].to_numpy()
        time = self.raw_data["Time (min)"].to_numpy()

        coeffs = np.polyfit(time, values, deg=degree)
        baseline_func = np.poly1d(coeffs)
        self.baseline = baseline_func(time)

        peaks, _ = find_peaks(values)
        for peak in peaks:
            left_thresh_idx = max(peak - 1, 0)
            right_thresh_idx = min(peak + 1, len(values) - 1)
            retention_time = time[peak]

            height_corrected = values[peak] - self.baseline[retention_time]
            peak_area = simps(
                y=values[left_thresh_idx:right_thresh_idx + 1] - self.baseline[retention_time],
                x=time[left_thresh_idx:right_thresh_idx + 1],
            )

            self.peaks.append(
                Peak(
                    time[left_thresh_idx],
                    time[right_thresh_idx],
                    height_corrected,
                    retention_time
                )
            )

    def calculate_peak_area(self, peak: Peak):
        values = self.raw_data["Value (EU)"].to_numpy()
        time = self.raw_data["Time (min)"].to_numpy()
        return simps(
                y=values[left_thresh_idx:right_thresh_idx + 1] - self.baseline[peak.retention_time],
                x=time[left_thresh_idx:right_thresh_idx + 1],
            )


if __name__ == "__main__":
    filepath = Path(__file__).parent.parent / "data" / "IgG Vtag 1_ACQUITY FLR ChA.txt"
    chrom = Chromatogram(filepath, 0.4, 1)
    chrom.detect_peaks("local")
    peaks1 = chrom.peaks
    chrom.detect_peaks("global")
    peaks2 = chrom.peaks
    print()
