from ..gas import Gas
import numpy as np


def source(source_strength: float):
    P = np.zeros((10, 5))
    P[-1, :] = source_strength
    return P


def sink(sink_strength: float):
    S = np.zeros((10, 5))
    S[0, :] = sink_strength
    return S
