from parser.peak import Peak

import pandas as pd

sample_data = pd.DataFrame({
    'Time (min)': [0.0, 0.5, 1.0, 1.5, 2.0],
    'Value (EU)': [0, 10, 20, 10, 0]
})


def test_peak_init():
    peak = Peak(left_thresh=0.0, right_base_idx=2.0, height=20, retention_time=1.0, data=sample_data)
    assert peak.left_thresh == 0.0
    assert peak.right_thresh == 2.0
    assert peak.height == 20
    assert peak.retention_time == 1.0
    assert peak.data.equals(sample_data)
