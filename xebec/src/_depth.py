import biom
import numpy as np
import pandas as pd


def map_quartiles(depths):
    q1, q2, q3 = np.quantile(depths, [0.25, 0.50, 0.75])
    quartiles = []
    for value in depths:
        if value < q1:
            q = "q1"
        elif q1 <= value < q2:
            q = "q2"
        elif q2 <= value < q3:
            q = "q3"
        else:
            q = "q4"
        quartiles.append(q)
    return quartiles
