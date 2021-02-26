from .roots import Roots
import numpy as np
import math
from typing import List
from .solute import Solute
from .constants import MOLAR_GAS_CONSTANT, M_SUCROSE, RHO_SUCROSE,\
    MAX_ELEMENT_COLUMNS, TEMPERATURE, VISCOSITY_WATER


class Tree:
    """ Model of a tree.

    Provides properties and functionality for saving and editing the modelled tree. Arguments whose
    type is List[float] or List[List[float]] are converted to numpy.ndarray with numpy.asarray method.
    Thus, also numpy.ndarray is a valid type for these arguments.

    For arguemnts whose type is List[float] (except for initial_radius) the length of the arguments must
    be equal to num_elements. The order of the list should be from the top of the tree (the first item) to the
    bottom of the tree (the last item)

    For arguments whose type is List[List[float]] the length of the arguemnts must be equal to num_elements
    and each sub list must contain two elements, one for the xylem and one for the phloem in this order. The
    order of the sub lists should be from the top of the tree (the first sub list) to the bottom of the tree
    (the last sub list).

    Args:
        height (float): total tree height (:math:`m`)
        initial_radius (List[float] or numpy.ndarray): the radius of the heartwood, xylem and phloem (:math:`m`)
            in this order. See from the [modelled system](modelled_system.html), how the radii should be given.
            Only three values can be given and the radius of each element is set to be the same in the tree
            initialization.
        num_elements (int): number of vertical elemenets in the tree.
            The height of an element is determined by
            :math:`\\text{element height} = \\frac{\\text{tree height}}{\\text{number of elements}}`
        transpiration_profile (List[float] or numpy.ndarray): The rate of transpiration (:math:`\\frac{kg}{s}`) in the
            xylem. The length of the list must be equal to num_elements and the order is from the top of the tree
            (first value) in the list to the bottom of the tree (last value in the list).
        photosynthesis_profile (List[float]): The rate of photosynthesis (:math:`\\frac{mol}{s}`). Currently this
            variable is not used and the rate of photosynthesis should be equal to the sugar_loading_profile.
        sugar_profile (List[float]] or numpy.ndarray): The initial sugar (sucrose) concentration in the phloem
            (:math:`\\frac{mol}{m^3}`)
        sugar_loading_profile (List[List[float]] or numpy.ndarray): the rate at which sugar concentration increases
            in each phloem element (:math:`\\frac{mol}{s}`)
        sugar_unloading_profile (List[float] or numpy.ndarray): The initial sugar unloading rate (the rate at which the
            sugar concentration decreases in a given phloem element) (:math:`\\frac{mol}{s}`). The unloading rate is
            updated in [src.odefun.odefun](index.html#src.odefun.odefun).
        sugar_target_concentration (float): the target concentration after which the sugar unloading
            starts (:math:`\\frac{mol}{m^3}`)
        sugar_unloading_slope (float): the slope parameter for unloading (see
            [Nikinmaa et. al., (2014)](https://academic.oup.com/aob/article/114/4/653/2769025)).
        axial_permeability_profile (List[List[float]] or numpy.ndarray): axial permeabilities of both xylem and phloem
            (:math:`m^2`)
        radial_hydraulic_conductivity_profile (List[float]] or numpy.ndarray): radial hydraulic conductivity between the
            xylem and the phloem (:math:`\\frac{m}{Pa \\: s}`)
        elastic_modulus_profile (List[List[float]] or numpy.ndarray): Elastic modulus of every element (:math:`Pa`).
        ground_water_potential (float): The water potential in the soil. This is used to calculate the sap flux between
            soil and the bottom xylem element.

    Attributes:
        height (float): total tree height (:math:`m`)
        num_elements (float): number of vertical elemenets in the tree.
        transpiration_rate (numpy.ndarray(dtype=float, ndim=2) [tree.num_elements, 1]): The rate of transpiration
            (:math:`\\frac{kg}{s}`) in the xylem.
        photosynthesis_rate (numpy.ndarray(dtype=float, ndim=2) [tree.num_elements, 1]): The rate of photosynthesis
            (:math:`\\frac{mol}{s}`). Currently this variable is not used.
        sugar_loading_rate (numpy.ndarray(dtype=float, ndim=2) [tree.num_elements, 1]): The rate at which sugar
            concentration increases in each phloem element (:math:`\\frac{mol}{s}`).
        sugar_unloading_rate (numpy.ndarray(dtype=float, ndim=2) [tree.num_elements, 1]): The rate at which the
            sugar concentration decreases in a given phloem element (:math:`\\frac{mol}{s}`).
        sugar_target_concentration (float): The target concentration after which the sugar unloading
            starts (:math:`\\frac{mol}{m^3}`).
        sugar_unloading_slope (float): The slope parameter for unloading (see
            [Nikinmaa et. al., (2014)](https://academic.oup.com/aob/article/114/4/653/2769025)).
        solutes (numpy.ndarray(dtype=src.solute.Solute, ndim=2) [tree.num_elements, 2]): Array of
            src.solute.Solute which contain the solutes in the sap of xylem and phloem.
        axial_permeability (numpy.ndarray(dtype=float, ndim=2) [tree.num_elements, 2]): Axial permeabilities of both
            xylem and phloem (:math:`m^2`).
        radial_hydraulic_conductivity (numpy.ndarray(dtype=float, ndim=2) [tree.num_elements, 1]): Radial hydraulic
            conductivity between the xylem and the phloem (:math:`\\frac{m}{Pa \\: s}`).
        elastic_modulus (numpy.ndarray(dtype=float, ndim=2) [tree.num_elements, 2]): Elastic modulus of every element
            (:math:`Pa`).
        ground_water_potential (float): The water potential in the soil.
        pressure (numpy.ndarray(dtype=float, ndim=2) [tree.num_elements, 2]): Pressure of each element (:math:`Pa`)
        element_radius (numpy.ndarray(dtype=float, ndim=2) [tree.num_elements, 3]): Radius of each element (:math:`m`)
        element_height (numpy.ndarray(dtype=float, ndim=2) [tree.num_elements, 2]): Height of each element (:math:`m`)
        viscosity (numpy.ndarray(dtype=float, ndim=2) [tree.num_elements, 2]): The dynamic viscosity of each element
            (:math:`Pa \\: s`)

    """

    def __init__(self, height: float, element_height: List[float], initial_radius: List[float],
                 num_elements: int, transpiration_profile: List[float],
                 photosynthesis_profile: List[float], sugar_profile: List[float],
                 sugar_loading_profile: List[float], sugar_unloading_profile: List[float],
                 sugar_target_concentration: float, sugar_unloading_slope: float,
                 axial_permeability_profile: List[List[float]],
                 radial_hydraulic_conductivity_profile: List[float],
                 elastic_modulus_profile: List[List[float]],
                 roots: Roots):

        # for all arrays/lists the first column is for xylem and the second for phloem
        # all the arrays/lists have num_elements rows

        self.height: float = height  # unit: m

        self.num_elements: int = num_elements  # number of total elements

        self.transpiration_rate: np.ndarray = np.asarray(transpiration_profile).reshape(self.num_elements, 1)

        self.photosynthesis_rate: np.ndarray = np.asarray(photosynthesis_profile).reshape(self.num_elements, 1)

        self.solutes: np.ndarray = np.asarray([[Solute('', 0, 0, 0),
                                                Solute('Sucrose', M_SUCROSE, RHO_SUCROSE, s_conc)]
                                               for s_conc in sugar_profile])

        self.sugar_loading_rate: np.ndarray = np.asarray(sugar_loading_profile).reshape(self.num_elements, 1)

        self.sugar_unloading_rate: np.ndarray = np.asarray(sugar_unloading_profile).reshape(self.num_elements, 1)

        self.sugar_target_concentration: float = sugar_target_concentration

        self.sugar_unloading_slope: float = sugar_unloading_slope

        self.axial_permeability: np.ndarray = np.asarray(axial_permeability_profile).reshape(self.num_elements, 2)

        self.radial_hydraulic_conductivity: np.ndarray = np.asarray(radial_hydraulic_conductivity_profile)\
            .reshape(self.num_elements, 1)

        self.elastic_modulus: np.ndarray = np.asarray(elastic_modulus_profile).reshape(self.num_elements, 2)

        self.roots = roots

        # set root elements as the roots.num_elements lowest elements

        self.root_elements = np.arange(self.num_elements - self.roots.num_elements, self.num_elements)
        self.tree_elements = np.arange(0, self.num_elements-self.roots.num_elements)
        # initialize pressure to be 0
        self.pressure = np.asarray([0 for i in range(self.num_elements)]).reshape(self.num_elements, 1)\
            .repeat(2, axis=1)

        # self.pressure[:, 1] = (self.sugar_concentration_as_numpy_array()*MOLAR_GAS_CONSTANT*TEMPERATURE
        #                        ).reshape(self.num_elements,)

        # calculate radius and height for every element in the tree
        self.element_radius: np.ndarray = np.asarray([initial_radius]*self.num_elements)

        self.element_height: np.ndarray = element_height.reshape(self.num_elements, 1)

        # initialize viscosity to be viscosity of water
        self.viscosity: np.ndarray = np.asarray([[VISCOSITY_WATER, VISCOSITY_WATER]]*self.num_elements)

        # update phloem viscosity due to sucrose
        self.update_sap_viscosity()

    def sugar_concentration_as_numpy_array(self) -> np.ndarray:
        """ Transforms the phloem sugar concentration in [self.solutes](index.html#src.tree.Tree.solutes)
        into numpy.ndarray.

        Returns:
            numpy.ndarray(dtype=float, ndim=2) [self.num_elements, 1]: The sugar concentration in the phloem.
            (:math:`\\frac{mol}{m^3}`)
        """
        get_concentration = np.vectorize(lambda s: s.concentration)
        return get_concentration(self.solutes[:, 1]).reshape(self.num_elements, 1)

    def update_sugar_concentration(self, new_concentration: np.ndarray) -> None:
        """ Sets the sugar concentration in [self.solutes](index.html#src.tree.Tree.solutes) to new_concentration.

        Args:
            new_concentration (numpy.ndarray(dtype=float, ndim=2)[self.num_elements,1]): new concentration values.
                the order is from top of the tree (first element, new_concentration[0]) to bottom of the tree
                (last element, new_concentration[self.num_elements-1]) (:math:`\\frac{mol}{m^3}`)
        """

        def update_concentration(ele, new_value):
            ele.concentration = new_value
        update = np.vectorize(update_concentration)
        update(self.solutes[:, 1], new_concentration.reshape(self.num_elements,))

    def element_area(self, ind: List[int] = None, column: int = 0) -> np.ndarray:
        """ Calculates the base area of the xylem or the phloem.

        Args:
            ind (List[int] or numpy.ndarray(dtype=int, ndim=1), optional): the indices of the elements
                for which the base area is calculated. If no ind is given, the base area is calculated
                for every element.
            column (int, optional): The column in the tree grid for which the base area is calculated.
                use column=0 for the xylem and column=1 for the phloem. If not column is given returns
                the base area for the xylem.

        Returns:
            numpy.ndarray(dtype=float, ndim=2) [len(ind) or self.num_elements, 1]: Base area of either
            the xylem or the phloem (:math:`m^2`)
        """
        num_elements = self.num_elements
        if ind is None:
            ind = []
        radii = self.element_radius
        if len(ind) > 0:
            radii = radii[ind, :]
            num_elements = len(ind)
        if column < MAX_ELEMENT_COLUMNS:
            radii = radii[:, 0:column+2]

        # calculate areas for every element column until desired column
        total_area: np.ndarray = (np.sum(radii, axis=1))**2
        total_area = total_area.reshape(num_elements, 1)

        inner_radius: np.ndarray = np.sum(radii[:, 0:column+1], axis=1)
        inner_radius = inner_radius.reshape(num_elements, 1)
        return math.pi*(total_area - inner_radius**2)

    def element_volume(self, ind: List[int] = None, column: int = 0) -> np.ndarray:
        """ Calculates the volume of the xylem or the phloem.

        Args:
            ind (List[int] or numpy.ndarray(dtype=int, ndim=1), optional): the indices of the elements
                for which the volume is calculated. If no ind is given, the volume is calculated
                for every element.
            column (int, optional): The column in the tree grid for which the volume is calculated.
                use column=0 for the xylem and column=1 for the phloem. If not column is given returns
                the volume for the xylem.

        Returns:
            numpy.ndarray(dtype=float, ndim=2) [len(ind) or self.num_elements, 1]: Volume of either
            the xylem or the phloem (:math:`m^3`)
        """
        # TODO: refactor the self.element_radius finding to own function (used multiple times)
        if ind is None:
            ind = []
        heights = self.element_height.reshape(self.num_elements, 1)

        if len(ind) > 0:
            heights = heights[ind, :]

        return self.element_area(ind, column) * heights

    def cross_sectional_area(self, ind: List[int] = None) -> np.ndarray:
        """ Calculates the cross-sectional area between the xylemn and the phloem.

        The cross sectional area is equal to lateral surface area of the xylem.

        Args:
            ind (List[int] or numpy.ndarray(dtype=int, ndim=1), optional): the indices of the elements
                for which the cross-sectinoal area is calculated. If no ind is given, the
                cross-sectional area is calculated for every element.

        Returns:
            numpy.ndarray(dtype=float, ndim=2) [len(ind) or self.num_elements, 1]: Cross-sectional area
            between the xylem and phloem elements (:math:`m^2`)
        """
        if ind is None:
            ind = []

        element_heights = self.element_height
        element_radii = np.sum(self.element_radius[:, 0:2], axis=1).reshape(self.num_elements, 1)
        if len(ind) > 0:
            element_heights = element_heights[ind]
            element_radii = element_radii[ind]

        return element_heights*2.0*math.pi*element_radii

    def update_sap_viscosity(self) -> None:
        """ Calculates and sets the viscosity in the phloem according to the sugar concenration.

            The sap viscosity is calculated according to Morrison (2002)

            .. math::
                \\eta = \\eta_w \\exp{\\frac{4.68 \\cdot 0.956 \\Phi_s}{1-0.956 \\Phi_s}}

            where
            * :math:`\\eta_w`: Dynamic viscosity of water (:math:`\\eta_w \\approx 0.001`)
            * :math:`\\Phi_s`: Volume fraction of sugar (sucrose) in the phloem sap.

            References:
                Morison, Ken R. "Viscosity equations for sucrose solutions: old and new 2002."
                Proceedings of the 9th APCChE Congress and CHEMECA. 2002.
        """
        # TODO: refactor into solution class

        # calculate sugar volume in sap
        sugar_volume_fraction: np.ndarray = self.sugar_concentration_as_numpy_array() * M_SUCROSE / RHO_SUCROSE
        sugar_volume_fraction = sugar_volume_fraction.reshape(self.num_elements, 1)
        sugar_volume_fraction = np.minimum(sugar_volume_fraction, 0.7)
        viscosity: np.ndarray = VISCOSITY_WATER * np.exp(4.68 * 0.956 * sugar_volume_fraction /
                                                         (1 - 0.956 * sugar_volume_fraction))

        self.viscosity[:, 1] = viscosity.reshape(self.num_elements,)
