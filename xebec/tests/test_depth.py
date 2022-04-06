import numpy as np

from xebec.src._depth import map_quartiles


def test_map_quartiles():
    depths = np.linspace(0, 100, 100)
    quartiles = map_quartiles(depths)
    for q in ["q1", "q2", "q3", "q4"]:
        assert quartiles.count(q) == 25
