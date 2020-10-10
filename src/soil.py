import numpy as np


class Soil:
    """"
        Model of a soil. Each array should have the same length except for the depth
        array whose length should be one higher than all the other arrays
    """

    def __init__(self, depth, pressure, hydraulic_conductivity):
        self.depth: np.ndarray = np.asarray(depth).reshape(len(depth), 1)
        self.num_elements: float = len(self.depth)-1
        self.pressure: np.ndarray = np.asarray(pressure).reshape(self.num_elements, 1)
        self.hydraulic_conductivity: np.ndarray = np.asarray(hydraulic_conductivity).reshape(self.num_elements, 1)

    def layer_thickness(self) -> np.ndarray:
        return np.diff(self.depth, axis=0).reshape(self.num_elements, 1)
