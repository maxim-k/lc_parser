from io import StringIO
from pathlib import Path
from typing import Dict

import pandas as pd


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
        }
        self.raw_data: pd.DataFrame = pd.DataFrame()
        self._parse_file(filepath)

    def _parse_file(self, filepath: Path):
        with filepath.open("r") as file:
            content = file.read()

        metadata, data = content.split("Chromatogram Data:\n")
        current_section = "Other"

        for line in metadata.split("\n"):
            line = line.strip()
            if line in [
                "Injection Information:",
                "Chromatogram Data Information:",
                "Signal Parameter Information:",
            ]:
                current_section = line[:-1]
            elif line:
                key, value = line.split("\t")
                self.metadata[current_section][key] = value

        self.raw_data = pd.read_csv(StringIO(data), sep="\t")

    def detect_peaks(self) -> None:
        pass

    def calculate_elution_volume(self) -> None:
        pass
