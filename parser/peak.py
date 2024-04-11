class Peak:
    def __init__(
        self, left_thresh: float, right_thresh: float, height: float, baseline: float
    ):
        self.left_thresh = left_thresh
        self.right_thresh = right_thresh
        self.height = height
        self.baseline = baseline

    def __repr__(self):
        return (
            f"Peak(left_thresh={self.left_thresh}, right_thresh={self.right_thresh}, "
            f"height={self.height}, baseline={self.baseline})"
        )
