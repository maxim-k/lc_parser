class Peak:
    def __init__(
        self,
        left_thresh: float,
        right_thresh: float,
        height: float,
        retention_time: float,
    ):
        self.left_thresh = left_thresh
        self.right_thresh = right_thresh
        self.height = height
        self.retention_time = retention_time

    def __repr__(self):
        return (
            f"Peak(left_thresh={self.left_thresh}, right_thresh={self.right_thresh}, "
            f"height={self.height}, retention_time={self.retention_time}"
        )
