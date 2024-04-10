from typing import Dict
import pandas as pd
from pathlib import Path


class Chromatogram:
    def __init__(self, filepath: Path):
        if not isinstance(filepath, Path):
            raise ValueError("filepath must be a pathlib.Path object")
        if not filepath.is_file():
            raise FileNotFoundError(f"No file found at {filepath}")

        self.metadata: Dict[str, Dict[str, str]] = {
            'Injection Information': {},
            'Chromatogram Data Information': {},
            'Signal Parameter Information': {}
        }
        self.raw_data: pd.DataFrame = pd.DataFrame(columns=['Time (min)', 'Step (s)', 'Value (EU)'])
        self._parse_file(filepath)

    def _parse_file(self, filepath: Path) -> None:
        # Non None but whatever
        pass

    def detect_peaks(self) -> None:
        pass

    def calculate_elution_volume(self) -> None:
        pass

