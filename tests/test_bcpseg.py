from linear_segment import __version__, segment
import numpy as np


def test_version():
    assert __version__ == "1.0.0"


def test_segment():
    test_array = np.random.random(1000)
    test_array[100:200] += 1
    labels = np.repeat("a", len(test_array)).astype(str)

    segs = segment(test_array, labels)
    assert np.sum(segs.starts == np.array([0,100,200])) == 3
    assert np.sum(segs.ends == np.array([100,200,1000])) == 3