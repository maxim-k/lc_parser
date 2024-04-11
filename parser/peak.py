class Peak:
    def __init__(
        self,
        left_thresh: float,
        right_thresh: float,
        height: float,
        baseline: float,
        retention_time: float,
        elution_volume: float,
        area: float,
    ):
        self.left_thresh = left_thresh
        self.right_thresh = right_thresh
        self.height = height
        self.baseline = baseline
        self.retention_time = retention_time
        self.elution_volume = elution_volume
        self.area = area

    def __repr__(self):
        return (
            f"Peak(left_thresh={self.left_thresh}, right_thresh={self.right_thresh}, "
            f"height={self.height}, baseline={self.baseline}, retention_time={self.retention_time}, "
            f"elution_volume={self.elution_volume}, area={self.area})"
        )
