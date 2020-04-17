import numpy as np
import math
from typing import List
from .solute import Solute
from .constants import M_SUCROSE, PHLOEM_RADIUS, RHO_SUCROSE, RHO_WATER,\
    GRAVITATIONAL_ACCELERATION, HEARTWOOD_RADIUS, MAX_ELEMENT_COLUMNS, VISCOSITY_WATER,\
    XYLEM_PHLOEM_CONTACT_ANGLE, XYLEM_RADIUS, TEMPERATURE, MOLAR_GAS_CONSTANT


class Tree:
    def __init__(self, height, initial_radius, num_elements, transpiration_profile,
                 photosynthesis_profile, sugar_profile,
                 sugar_loading_profile, sugar_unloading_profile, sugar_target_concentration,
                 axial_permeability_profile,
                 radial_hydraulic_conductivity_profile,
                 elastic_modulus_profile,
                 ground_water_potential):

        # for all arrays/lists the first column is for xylem and the second for phloem
        # all the arrays/lists have num_elements rows

        self.height: float = height  # unit: m

        self.num_elements: int = num_elements  # number of tree elements

        self.initial_radius: np.ndarray = np.asarray(initial_radius)  # initial radius of the xylem and phloem

        # (num_elements,1) array of transpiration rates in the xylem
        # unit: kg/s
        self.transpiration_rate: np.ndarray = np.asarray(transpiration_profile).reshape(40, 1)

        # (num_elements,1) array of photosynth. rate in the phloem
        # unit: mol/s
        self.photosynthesis_rate: np.ndarray = np.asarray(photosynthesis_profile).reshape(40, 1)

        # (num_elements,1) array of sugar conc. at t=0s in phloem
        # unit: mol/m3
        self.solutes: np.ndarray = np.asarray([[Solute('', 0, 0, 0),
                                                Solute('Sucrose', M_SUCROSE, RHO_SUCROSE, s_conc)]
                                               for s_conc in sugar_profile])

        # (num_elements,1) array of sugar loading rates in phloem
        # unit: mol/s
        self.sugar_loading_rate: np.ndarray = np.asarray(sugar_loading_profile).reshape(40, 1)

        # (num_elements,1) array of sugar unloading rates in phloem
        # unit: mol/s
        self.sugar_unloading_rate: np.ndarray = np.asarray(sugar_unloading_profile).reshape(40, 1)

        self.sugar_target_concentration: float = sugar_target_concentration

        self.axial_permeability: np.ndarray = np.asarray(axial_permeability_profile)

        self.radial_hydraulic_conductivity: np.ndarray = np.asarray(radial_hydraulic_conductivity_profile)

        self.elastic_modulus: np.ndarray = np.asarray(elastic_modulus_profile)

        self.ground_water_potential: float = ground_water_potential

        # calculate what the pressure in the xylem will be
        # assuming that at the base of the tree water potential
        # equals ground water potential
        # order 1 = top of the tree, N = base of the tree

        self.pressure = np.asarray([self.ground_water_potential
                                    - i*RHO_WATER*GRAVITATIONAL_ACCELERATION*self.height/self.num_elements
                                    for i in range(self.num_elements)]).reshape(self.num_elements, 1)
        # reverse pressure so the order is correct (elemenent N has pressure equal to ground water potential)
        self.pressure = np.concatenate((np.flip(self.pressure),
                                        np.flip(self.pressure)),
                                       axis=1)

        # calculate radius and height for every element in the tree
        self.element_radius: np.ndarray = np.asarray([initial_radius]*self.num_elements)

        self.element_height: np.ndarray = np.asarray(
            [self.height/self.num_elements]*self.num_elements).reshape(self.num_elements, 1)

        # initialize viscosity to be viscosity of water
        self.viscosity: np.ndarray = np.asarray([[VISCOSITY_WATER, VISCOSITY_WATER]]*self.num_elements)
        # update phloem viscosity due to sucrose
        self.update_sap_viscosity()

    def sugar_concentration_as_numpy_array(self) -> np.ndarray:
        get_concentration = np.vectorize(lambda s: s.concentration)
        return get_concentration(self.solutes[:, 1])

    def update_sugar_concentration(self, new_concentration: np.ndarray) -> None:

        def update_concentration(ele, new_value):
            ele.concentration = new_value
        update = np.vectorize(update_concentration)
        update(self.solutes[:, 1], new_concentration.reshape(self.num_elements,))

    def element_area(self, ind: List[int] = None, column: int = 0) -> np.ndarray:
        """ returns element areas specified in parameter ind and column of self.elements.

            If no ind is given returns the areas for every element.
            If no column is given returns the areas in column 0 of self.elements (the xylem)
        """
        if ind is None:
            ind = []
        radii = self.element_radius
        if len(ind) > 0:
            radii = radii[ind, :]
        if column < MAX_ELEMENT_COLUMNS:
            radii = radii[:, 0:column+1]

        # calculate areas for every element column until desired column
        total_area: np.ndarray = (np.sum(radii[:, 0:column+1], axis=1)+HEARTWOOD_RADIUS)**2
        total_area = total_area.reshape(self.num_elements, 1)
        if column == 0:
            inner_radius: np.ndarray = np.repeat(HEARTWOOD_RADIUS, repeats=self.num_elements)\
                .reshape(self.num_elements, 1)
        else:
            inner_radius: np.ndarray = radii[:, 0:column] + HEARTWOOD_RADIUS
        return math.pi*(total_area - inner_radius**2)

    def element_volume(self, ind: List[int] = None, column: int = 0) -> np.ndarray:
        """ returns element volumes specified in parameter ind and column of self.elements.

            If no ind is given returns the volume of every element.
            If no column is given returns the volumes in column 0 of self.elements (the xylem)
        """
        # TODO: refactor the self.element_radius finding to own function (used multiple times)
        if ind is None:
            ind = []
        heights = self.element_height.reshape(self.num_elements, 1)

        if len(ind) > 0:
            heights = heights[ind, :]

        return self.element_area(ind, column) * heights

    def cross_sectional_area(self, ind: List[int] = None) -> np.ndarray:
        """ calculates the cross sectional area between xylem and phloem

            if no ind is given returns the cross sectional area for every axial element
        """
        if ind is None:
            ind = []
        phloem_element_height = self.element_height
        if len(ind) > 0:
            phloem_element_height = phloem_element_height[ind]

        return XYLEM_PHLOEM_CONTACT_ANGLE/(2.0*math.pi)*phloem_element_height

    def update_sap_viscosity(self) -> None:
        """ Updates the viscosity in column 1 of self.elements (the phloem)
        """
        # TODO: refactor into solution class

        # calculate sugar volume in sap
        sugar_volume_fraction: np.ndarray = np.sum(self.sugar_concentration_as_numpy_array()) * M_SUCROSE / RHO_SUCROSE

        viscosity: np.ndarray = VISCOSITY_WATER * np.exp(4.68 * 0.956 * sugar_volume_fraction /
                                                         (1 - 0.956 * sugar_volume_fraction))
        self.viscosity[:, 1] = viscosity
