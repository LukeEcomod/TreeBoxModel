from typing import List
import numpy as np
from .soil import Soil
from .solute import Solute
from .constants import M_SUCROSE, RHO_SUCROSE


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
            typical soil conductance in a layer (:math:`s`)

    """

    def __init__(self, rooting_depth: float,
                 area_density: np.ndarray,
                 effective_radius: np.ndarray,
                 conductance_to_soil: float):

        self.num_elements = len(area_density)
        self.area_density: np.ndarray = area_density
        self.effective_radius: np.ndarray = effective_radius
        self.conductance_to_soil: np.ndarray = conductance_to_soil
        self.rooting_depth = rooting_depth

    def root_area_index(self, soil: Soil) -> float:
        """ Calculates the root area index i.e.

            :math:`RAI = \\sum_{i=1}^{N} B_i \\Delta z_i`

            where :math:`B_i` is the root area density and :math:`\\Delta z_i` is the thickness of layer i
        Args:
            soil (Soil): Instance of the soil unit class. The soil layer thicknesses is used in calculation.

        Returns:
            float: The root area index (:math:`\\frac{m^2}{m^2}`)
        """
        return np.sum(self.area_density*soil.layer_thickness())

    def root_conductance(self, soil: Soil, ind: List[int] = None) -> np.ndarray:
        """ calculates to conductance from soil-root interface up to the xylem of the tree which is
            :math:`k_{r,i}` in Volpe et al., (2013).

            :math:`k_{r,i} = \\frac{B_i \\Delta z_i}{\\beta}`

            where :math:`B_i` is the root area density, :math:`\\Delta z_i` is the thickness of layer i
            and :math:`\\beta` is "typical root to soil conductance" in Volpe et. al., (2013)

        Args:
            soil (Soil): Instance of the soil unit class. The soil layer thicknesses is used in calculation.
            ind (List[int] or numpy.ndarray(dtype=int, ndim=1), optional): the indices of the elements
                for which conductance is calculated. If no ind is given, the
                cross-sectional area is calculated for every element.

        Returns:
            numpy.ndarray(dtype=float, ndim=2) [len(ind) or self.num_elements, 1]: conductance from
            root-soil interface to root xylem (:math:`s^{-1}`)

        References:
            Volpe, V. et. al., "Root controls on water redistribution and carbon uptake in the soil–plant
            system under current and future climate", Advances in Water Resources, 60, 110-120, 2013.
        """
        dz = soil.layer_thickness()

        if ind is None:
            ind = []

        if(len(ind) > 0):
            return self.area_density[ind] * dz[ind]/self.conductance_to_soil\
                .reshape(len(ind), 1)
        else:
            return self.area_density * dz / self.conductance_to_soil\
                .reshape(self.num_elements, 1)
