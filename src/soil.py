import numpy as np


class Soil:
    """"Model of a soil. Each array should have the same length except for the depth
        array whose length should be one higher than all the other arrays
    """

    def __init__(self, depth, water_potential, hydraulic_conductivity):
        self.depth = depth
        self.water_potential = water_potential
        self.hydraulic_conductivity = hydraulic_conductivity
        self.num_elements = len(self.thickness)
        self.pressure = np.asarray([0 for i in range(self.num_elements)]).reshape(self.num_elements, 1)

    def layer_thickness(self):
        return np.diff(self.depth)
