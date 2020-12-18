from ..tree import Tree
from ..gas import Gas
import numpy as np


def source(gas: Gas, tree: Tree):
    P = np.zeros((gas.na, gas.nr))
    P[-1, :] = 1.5
    return P
