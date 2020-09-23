import numpy as np
from .soil import Soil


class Roots:
    """ Model of the root system of the tree. The number of elements (layers) is calculated from the length
        of the area density variable. All the arrays should have the same length.

    Args:
        rooting_depth (float): depth of the root system in the soil (:math:`m`)
        area_density (List[float] or numpy.ndarray): surface area of roots in given layer per layer
            volume (:math:`\\frac{m^2}{m^3}`)
        effective radius (List[float] or numpy.ndarray): effective horizontal root radius in a layer
        conductance_to_soil (List[float] or numpy.ndarray): typical soil conductance in a layer CHECK UNIT!
        soil (Soil): instance of the soil class where the roots lie

    Attributes:
        rooting_depth (float): depth of the root system in the soil (:math:`m`)
        area_density (numpy.ndarray(dtype=float, ndim=2) [self.num_elements, 1]):
            surface area of roots in given layer per layer volume (:math:`\\frac{m^2}{m^3}`)
        effective radius (numpy.ndarray(dtype=float, ndim=2) [self.num_elements, 1]):
            effective horizontal root radius in a layer
        conductance_to_soil (numpy.ndarray(dtype=float, ndim=2) [self.num_elements, 1]):
            typical soil conductance in a layer CHECK UNIT!
        soil (Soil): instance of the soil class where the roots are

    """

    def __init__(self, rooting_depth, area_density, effective_radius, conductance_to_soil, soil):
        self.rooting_depth: float = rooting_depth
        self.num_elements = len(area_density)
        self.area_density: np.ndarray = area_density
        self.effective_radius: np.ndarray = effective_radius
        self.conductance_to_soil: np.ndarray = conductance_to_soil
        self.soil: Soil = soil
        self.pressure = np.asarray([0 for i in range(self.num_elements)]).reshape(self.num_elements, 1)

    def root_area_index(self) -> np.ndarray:
        """ Calculates the root area index i.e.

            :math:`RAI = \\sum_{i=1}^{N} B_i \\Delta z_i`

            where :math:`B_i` is the root area density and :math:`\\Delta z_i` is the thickness of layer i

        Returns:
            numpy.ndarray(dtype=float, ndim=2) [self.soil.num_elements, 1]: The root area index
            (unitless)
        """
        return np.sum(self.area_density*self.soil.layer_thickness())
