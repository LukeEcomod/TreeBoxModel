from typing import List
from .roots import Roots
import numpy as np


class Soil:
    """"
        Model of a soil. Each array should have the same length except for the depth
        array whose length should be one higher than all the other arrays
    """

    def __init__(self, depth, water_potential, hydraulic_conductivity):
        self.depth = depth
        self.water_potential = water_potential
        self.hydraulic_conductivity = hydraulic_conductivity
        self.num_elements = len(self.layer_thickness())

    def layer_thickness(self):
        return np.diff(self.depth)

    def soil_conductance(self, roots: Roots, ind: List[int] = None) -> np.ndarray:
        """ Calcualtes the conductance in layer i which is :math:`k_{s,i}` in Volpe et al., (2013)
            which is the horizontal conductance in soil to the soil-root interface.

            :math:`k_{s,i} = \\alpha K_i B_i`

            where :math:`K_i` is the water hydraulic conductivity in layer i,
            :math:`B_i` is the root area density in layer i and

            :math:`\\alpha = \\left(\\frac{L}{2RAI \\cdot r_i}\\right)^{1/2}`

            where :math:`L` is the rooting depth, RAI is the root area index and :math:`r_i`
            is the effective root radius in layer i.
            TODO: check that effective radius is the average root radius

        Args:
            roots (Roots): instance of the root class. Root area index, rooting depth, root effective radius and
                root area density are used in the calculations.

        Returns:
            numpy.ndarray(dtype=float, ndim=2) [len(ind) or self.num_elements, 1]: conductance from
            root-soil interface to root xylem (:math:`s^{-1}`)

        References:
            Volpe, V. et. al., "Root controls on water redistribution and carbon uptake in the soil–plant
            system under current and future climate", Advances in Water Resources, 60, 110-120, 2013.
        """

        if ind is None:
            ind = []

        if len(ind) > 0:
            return ((roots.rooting_depth/(2*roots.root_area_index()*roots.effective_radius[ind]))**(1/2)
                    * self.hydraulic_conductivity[ind]*roots.area_density[ind]).reshape(len(ind), 1)
        else:
            return ((roots.rooting_depth/(2*roots.root_area_index()*roots.effective_radius))**(1/2)
                    * self.hydraulic_conductivity*roots.area_density).reshape(self.num_elements, 1)
