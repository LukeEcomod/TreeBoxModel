from typing import List
import numpy as np
from .soil import Soil


class Roots:
    """ Model of the root system of the tree. The number of elements (layers) is calculated from the length
        of the area density variable. All the arrays should have the same length.

    Args:
        rooting_depth (float): depth of the root system in the soil (:math:`m`)
        area_density (List[float] or numpy.ndarray): surface area of roots in given layer per tree
        and per soil layer thickness (:math:`\\frac{m^2}{m}`)
        effective radius (List[float] or numpy.ndarray): effective horizontal root radius in a layer
        soil_conductance_scale (List[float] or numpy.ndarray): typical soil conductance in a layer.
        total_root_area (float): Lateral surface area of roots. Used in calculating root area density.
        num_elements (int): number of elements in the root zone.
        RAI (float):  Root area index of the stand (:math:`m^2 m^{-2}`)

    Attributes:
        rooting_depth (float): depth of the root system in the soil (:math:`m`)
        area_density (numpy.ndarray(dtype=float, ndim=2) [self.num_elements, 1]):
            surface area of roots in given layer per tree and per soil layer thickness (:math:`\\frac{m^2}{m}`)
        effective radius (numpy.ndarray(dtype=float, ndim=2) [self.num_elements, 1]):
            effective horizontal root radius in a layer
        soil_conductance_scale (numpy.ndarray(dtype=float, ndim=2) [self.num_elements, 1]):
            typical soil to root xylem conductance in a layer (:math:`s`)
        RAI (float): Root area index of the stand (:math:`m^2 m^{-2}`)
    """

    def __init__(self, rooting_depth: float,
                 area_density: np.ndarray,
                 effective_radius: np.ndarray,
                 soil_conductance_scale: float,
                 total_root_area: float,
                 num_elements: int,
                 RAI: float ):

        self.num_elements = num_elements
        self.area_density: np.ndarray = area_density
        self.effective_radius: np.ndarray = effective_radius
        self.soil_conductance_scale: float = soil_conductance_scale
        self.total_root_area: float = total_root_area
        self.rooting_depth: float = rooting_depth
        self.RAI = RAI

    def root_area_index(self, soil: Soil) -> float:
        """ Calculates the root area index i.e.

            .. math::
                RAI = \\frac{\\sum_{i=1}^{N} B_i \\Delta z_i}}

            where :math:`B_i` is the root area density, :math:`\\Delta z_i` is the thickness of layer i.

            NB! the soil depth must match the depth for which Roots.area_density is given.
            NB2! This function is no longer used
        Args:
            soil (Soil): Instance of the soil unit class. The soil layer thicknesses is used in calculation.

        Returns:
            float: The root area index (:math:`\\frac{m^2}{m^2}`)
        """
        ind, _ = np.where(soil.depth() <= self.rooting_depth)
        dz = soil.layer_thickness[ind].reshape(len(ind), 1)
        return np.sum(self.area_density*dz)

    def root_conductance(self, soil: Soil, ind: List[int] = None) -> np.ndarray:
        """ calculates to conductance from soil-root interface up to the xylem of the tree which is
            :math:`k_{r,i}` in Volpe et al., (2013).

            .. math::
                k_{r,i} = \\frac{B_i \\Delta z_i}{\\beta}

            where :math:`B_i` is the root area density, :math:`\\Delta z_i` is the thickness of layer i,
            :math:`\\beta` is typical root to soil conductance as presented in Volpe et. al., (2013).

        Args:
            soil (Soil): Instance of the soil unit class. The soil layer thicknesses is used in calculation.
            ind (List[int] or numpy.ndarray(dtype=int, ndim=1), optional): the indices of the elements
                for which conductance is calculated. If no ind is given, the
                cross-sectional area is calculated for every element.

        Returns:
            numpy.ndarray(dtype=float, ndim=2) [len(ind) or self.num_elements, 1]: conductance from
            root-soil interface to root xylem (:math:`m^2 s^{-1}`)

        References:
            Volpe, V. et. al., "Root controls on water redistribution and carbon uptake in the soil-plant
            system under current and future climate", Advances in Water Resources, 60, 110-120, 2013.
        """

        if ind is None:
            ind = []

        if len(ind) > 0:
            conductance = (self.area_density[ind] \
                * soil.layer_thickness[ind]/self.soil_conductance_scale)\
                .reshape(len(ind), 1)
        else:
            root_ind, _ = np.where(soil.depth() < self.rooting_depth)
            conductance = (self.area_density * soil.layer_thickness[root_ind].reshape(len(root_ind), 1)
                    / self.soil_conductance_scale).reshape(self.num_elements, 1)
        return conductance
    def soil_root_conductance(self, soil: Soil, ind: List[int] = None) -> np.ndarray:
        """ Calculates the conductance in layer i which is :math:`k_{s,i}` in Volpe et al., (2013)
            which is the horizontal conductance in soil to the soil-root interface.

            .. math::
                k_{s,i} = \\alpha K_i B_i

            where :math:`K_i` is the water hydraulic conductivity in layer i,
            :math:`B_i` is the root area density in layer i and

            :math:`\\alpha = \\left(\\frac{L}{2RAI \\cdot r_i}\\right)^{1/2}`

            where :math:`L` is the rooting depth, RAI is the root area index, :math:`r_i`
            is the effective root radius in layer i,

        Args:
            soil (Soil): Instance of the soil class. The soil layer thicknesses is used in calculation.
            ind (List[int] or numpy.ndarray(dtype=int, ndim=1), optional): the indices of the elements
                for which conductance is calculated. If no ind is given, the
                cross-sectional area is calculated for every element.

        Returns:
            numpy.ndarray(dtype=float, ndim=2) [len(ind) or self.num_elements, 1]: conductance from
            root-soil interface to root xylem (:math:`m^2 s^{-1}`)

        References:
            Volpe, V. et. al., "Root controls on water redistribution and carbon uptake in the soil-plant
            system under current and future climate", Advances in Water Resources, 60, 110-120, 2013.
        """

        if ind is None:
            ind = []

        if len(ind) > 0:
            result = ((self.rooting_depth /
                       (2*self.RAI*self.effective_radius[ind]))**(1/2)
                      * soil.hydraulic_conductivity[ind]*self.area_density[ind]).reshape(len(ind), 1)
        else:
            root_ind, _ = np.where(soil.depth() < self.rooting_depth)

            result = ((self.rooting_depth /
                       (2*self.RAI*self.effective_radius))**(1/2)
                      * soil.hydraulic_conductivity[root_ind].reshape(len(root_ind), 1)
                      * self.area_density).reshape(self.num_elements,1)

        result[np.isnan(result)] = 0.0
        return result

    def conductivity(self, soil: Soil) -> np.ndarray:
        """ Calculates the total conductivity :math:`g_i` in each root/soil layer.

        The total conductivity is calculated according to Volpe et al., (2013)

        .. math::
            g_i = \\frac{k_{s,i} k_{r,i}}{k_{s,i}+k_{r,i}}

        where
        * :math:`k_{s,i}`: Conductance from soil to near the root system in root layer i. See
            [soil_root_conductance](index.html#src.roots.Roots.soil_root_conductance) for details
        * :math:`k_{r,i}`: Conductance from near the root system to the root xylem in root layer i. See
            [root_conductance](index.html#src.roots.Roots.root_conductance) for details.

        Args:
            soil (Soil): Instance of the soil class.

        Returns:
            numpy.ndarray (dtype=float, ndim=2)[self.tree.num_elements, 1]: Total conductance from soil to
            the root xylem in all the root layers in units :math:`\\frac{m^2}{s}`.

        References:
            Volpe, V. et. al., "Root controls on water redistribution and carbon uptake in the soil-plant
            system under current and future climate", Advances in Water Resources, 60, 110-120, 2013.

        """

        ks = self.soil_root_conductance(soil)
        kr = self.root_conductance(soil)
        result = np.zeros((self.num_elements, 1))
        divisor = kr*ks
        divider = kr+ks
        ind = np.where(divider != 0)
        result[ind] = divisor[ind]/divider[ind]
        return result

    def layer_thickness(self, soil: Soil) -> np.ndarray:
        """ Returns the layer thickness from soil surface until the self.rooting_depth.

        Args:
            soil (Soil): Instance of the soil class.

        Returns:
            numpy.ndarray(dtype=float, ndim=2)[n,1]: Layer thickness in units :math:`m`. n
            is the amount of soil layers whose depth is smaller than the rooting depth of the root system.
        """

        root_ind, _ = np.where(soil.depth() < self.rooting_depth)
        return soil.layer_thickness[root_ind].reshape(len(root_ind), 1)

    def layer_depth(self, soil: Soil) -> np.ndarray:
        """ Returns the midpoint depth of every layer that has roots.

        Args:
            soil (Soil): Instance of the soil class.

        Returns:
            numpy.ndarray(dtype=float, ndim=2)[self.num_elements+1,1]: Layer depth relative to surface
            in units :math:`m`.
        """
        dz = self.layer_thickness(soil)
        length = np.concatenate(([0], dz.reshape(self.num_elements, )))
        cumulative_sum = np.cumsum(length).reshape(len(length), 1)
        return dz/2 + cumulative_sum[:-1]

    def soil_elements(self, soil: Soil) -> np.ndarray:
        """ Returns array of elements in the soil object whose depth is above rooting depth

        Args:
            soil (Soil): Instance of the soil class.

        Returns:
            numpy.ndarray(dtype=int, ndim=1)[n,]: array of elements in the soil object for which the layer depth
            is smaller than the rooting depth.
        """
        root_ind, _ = np.where(soil.depth() < self.rooting_depth)
        return root_ind
