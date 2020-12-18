from ..gas import Gas
import numpy as np


def source(gas: Gas):
    P = np.zeros((gas.na, gas.nr))
    P[-1, :] = 1.5
    return P


def sink(gas: Gas):
    S = np.zeros((gas.na, gas.nr))
    S[0, :] = 2
    return S
