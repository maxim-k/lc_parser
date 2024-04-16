import pandas as pd


class Peak:
    """
    Represents a chromatographic peak.

    :param left_thresh: The lower threshold for the peak detection.
    :param right_thresh: The upper threshold for the peak detection.
    :param height: The height of the peak.
    :param retention_time: The retention time at which the peak occurs.
    :param data: A pandas DataFrame containing the data points of the peak.
    """

    def __init__(
        self,
        left_thresh: float,
        right_thresh: float,
        height: float,
        retention_time: float,
        data: pd.DataFrame,
    ):
        self.left_thresh = left_thresh
        self.right_thresh = right_thresh
        self.height = height
        self.retention_time = retention_time
        self.data = data

    def __repr__(self):
        return (
            f"Peak(left_thresh={self.left_thresh}, right_thresh={self.right_thresh}, "
            f"height={self.height}, retention_time={self.retention_time}, "
            f"data=[{self.data.to_string(index=False)}])"
        )
