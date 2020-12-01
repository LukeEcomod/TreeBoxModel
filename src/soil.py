import numpy as np


class Soil:
    """
        Model of a soil. Each array should have the same length except for the depth
        array whose length should be one higher than all the other arrays

        Args:
            layer_thickness (List[float] or numpy.ndarray): Thickness of each layer in units :math:`m`.
            pressure (List[float] or numpy.ndarray): Pressure (soil water potential) in every layer in units :math:`Pa`
            hydraulic_conductivity (List[float] or numpy.ndarray): Hydraulic conductivity in each soil layer in units
                :math:`\\frac{m}{s}`

        Attributes:
            layer_thickness (numpy.ndarray(dtype=float, ndim=2) [self.num_elements, 1]): Thickness of each layer in
                units :math:`m`.
            pressure (numpy.ndarray(dtype=float, ndim=2) [self.num_elements, 1]): Pressure (soil water potential) in
                every layer in units :math:`Pa`
            hydraulic_conductivity (numpy.ndarray(dtype=float, ndim=2) [self.num_elements, 1]): Hydraulic conductivity
                in each soil layer in units :math:`\\frac{m}{s}`
            num_elements (float): Number of soil elements. Calculated from the length of layer_thickness argument.
    """

    def __init__(self, layer_thickness, pressure, hydraulic_conductivity):
        self.layer_thickness: np.ndarray = np.asarray(layer_thickness).reshape(len(layer_thickness), 1)
        self.num_elements: float = len(self.layer_thickness)
        self.pressure: np.ndarray = np.asarray(pressure).reshape(self.num_elements, 1)
        self.hydraulic_conductivity: np.ndarray = np.asarray(hydraulic_conductivity).reshape(self.num_elements, 1)

    def depth(self) -> np.ndarray:
        """ Returns the midpoint of every soil element """
        dz = np.concatenate(([0], self.layer_thickness.reshape(len(self.layer_thickness),)))
        cumulative_sum = np.cumsum(dz).reshape(len(dz), 1)
        result = self.layer_thickness/2 + cumulative_sum[:-1]
        return result.reshape(self.num_elements, 1)
