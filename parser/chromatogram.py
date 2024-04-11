from io import StringIO
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from .peak import Peak
from scipy.interpolate import interp1d
from scipy.signal import find_peaks


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
        self.baseline: np.ndarray = np.array([])
        self._parse_file(filepath)

    def _parse_file(self, filepath: Path):
        with filepath.open("r") as file:
            content = file.read().strip()

        metadata, data = content.split("Chromatogram Data:\n")
        current_section = "Other"  # TODO What if it is in the middle

        for line in metadata.split("\n"):
            line = line
            if line in [
                "Injection Information:",
                "Chromatogram Data Information:",
                "Signal Parameter Information:",
            ]:
                current_section = line[:-1]
            elif line:
                key, value = line.split("\t")
                self.metadata[current_section][key] = value

        self.raw_data = pd.read_csv(
            StringIO(data), sep="\t", na_values="n.a."
        )  # TODO check other NaN
        time = self.raw_data["Time (min)"].to_numpy()
        values = self.raw_data["Value (EU)"].to_numpy()
        self.baseline = interp1d(
            [time[0], time[-1]], [values[0], values[-1]], fill_value="extrapolate"
        )

    def detect_peaks(self):
        """
        Detects peaks and calculates their corrected heights based on local baseline correction.

        The method detects peaks in the 'Value (EU)' column, correcting each peak's height by subtracting the
        baseline value at its position. Peaks are returned as a list of `Peak` objects, each containing detailed peak
        information.

        :return: A list of `Peak` objects, each representing a detected peak with properties for left threshold,
        right threshold, corrected height, and baseline value. :rtype: list[Peak]

        **Example**::

            chromatogram = Chromatogram(filepath)
            peaks = chromatogram.detect_peaks()
            for peak in peaks:
                print(peak)
        """
        values = self.raw_data["Value (EU)"].to_numpy()
        time = self.raw_data["Time (min)"].to_numpy()

        peaks, _ = find_peaks(values)

        for peak in peaks:
            left_thresh_idx = peak - 1 if peak > 0 else peak
            right_thresh_idx = peak + 1 if peak < len(values) - 1 else peak

            left_thresh = time[left_thresh_idx]
            right_thresh = time[right_thresh_idx]

            # peak_baseline = self.baseline[peak]  # for global baseline
            baseline_values = interp1d(
                [left_thresh, right_thresh],
                [values[left_thresh_idx], values[right_thresh_idx]],
                fill_value="extrapolate",
            )
            peak_baseline = baseline_values(time[peak])

            height_corrected = values[peak] - peak_baseline

            self.peaks.append(
                Peak(left_thresh, right_thresh, height_corrected, peak_baseline)
            )

    def calculate_elution_volume(self) -> None:
        pass


if __name__ == "__main__":
    filepath = Path(__file__).parent.parent / "data" / "IgG Vtag 1_ACQUITY FLR ChA.txt"
    chrom = Chromatogram(filepath)
    print()
