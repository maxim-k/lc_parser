import pandas as pd


class Peak:
    """
    Represents a chromatographic peak.

    :param left_base_idx: The left index (start) of the detected peak.
    :param right_base_idx: The right index (end) of the detected peak.
    :param height: The height of the peak.
    :param retention_time: The retention time at which the peak occurs.
    :param data: A pandas DataFrame containing the data points of the peak.
    """

    def __init__(
        self,
        left_base_idx: float,
        right_base_idx: float,
        height: float,
        retention_time: float,
        data: pd.DataFrame,
    ):
        self.left_base_idx = left_base_idx
        self.right_base_idx = right_base_idx
        self.height = height
        self.retention_time = retention_time
        self.data = data

    def __repr__(self):
        return (
            f"Peak(left_thresh={self.left_thresh}, right_thresh={self.right_thresh}, "
            f"height={self.height}, retention_time={self.retention_time}, "
            f"data=[{self.data.to_string(index=False)}])"
        )
