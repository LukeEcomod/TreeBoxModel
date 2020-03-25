import numpy as np
import math
from typing import List
from .treeelement import TreeElement
from .solute import Solute
from .constants import M_SUCROSE, RHO_SUCROSE, RHO_WATER,\
    GRAVITATIONAL_ACCELERATION, HEARTWOOD_RADIUS, MAX_ELEMENT_COLUMNS


class Tree:
    def __init__(self, height, initial_radius, num_elements, transpiration_profile,
                 photosynthesis_profile, sugar_profile,
                 sugar_loading_profile, sugar_unloading_profile,
                 axial_permeability_profile,
                 radial_hydraulic_conductivity_profile,
                 elastic_modulus_profile,
                 ground_water_potential):

        # TODO: remove the variables from class definition that are not needed

        self.height: float = height  # unit: m

        self.num_elements: int = num_elements  # number of tree elements

        self.initial_radius = initial_radius  # initial radius of the xylem and phloem

        # (num_elements,1) array of transpiration rates in the xylem
        # unit: kg/s
        self.transpiration_profile: List[float] = transpiration_profile

        # (num_elements,1) array of photosynth. rate in the phloem
        # unit: mol/s
        self.photosynthesis_profile: List[float] = photosynthesis_profile

        # (num_elements,1) array of sugar conc. at t=0s in phloem
        # unit: mol/m3
        self.sugar_profile: List[float] = sugar_profile

        # (num_elements,1) array of sugar loading rates in phloem
        # unit: mol/s
        self.sugar_loading_profile: List[float] = sugar_loading_profile

        # (num_elements,1) array of sugar unloading rates in phloem
        # unit: mol/s
        self.sugar_unloading_profile: List[float] = sugar_unloading_profile

        self.axial_permeability_profile: List[float] = axial_permeability_profile

        self.radial_hydraulic_conductivity_profile: List[float] = radial_hydraulic_conductivity_profile

        self.elastic_modulus_profile: List[float] = elastic_modulus_profile

        self.ground_water_potential: float = ground_water_potential

        # calculate what the pressure in the xylem will be
        # assuming that at the base of the tree water potential
        # equals ground water potential
        # order 1 = top of the tree, N = base of the tree

        pressures = [self.ground_water_potential - i*RHO_WATER*GRAVITATIONAL_ACCELERATION*self.height/self.num_elements
                     for i in range(num_elements)]
        pressures.reverse()  # reverse so that the order is correct

        # set up dim = (num_elements, 2) list of TreeElements
        # columm 1 = xylem
        # column 2 = phloem
        self.elements: List[List[TreeElement]] = []
        for pressure,\
            transpiration_rate,\
            photosynthesis_rate,\
            sugar_concentration,\
            sugar_loading_rate,\
            sugar_unloading_rate,\
            axial_permeability,\
            radial_hydraulic_conductivity,\
            elastic_modulus in zip(
                pressures,
                self.transpiration_profile,
                self.photosynthesis_profile,
                self.sugar_profile,
                self.sugar_loading_profile,
                self.sugar_unloading_profile,
                self.axial_permeability_profile,
                self.radial_hydraulic_conductivity_profile,
                self.elastic_modulus_profile):

            sugar: Solute = Solute(M_SUCROSE, RHO_SUCROSE, sugar_concentration)
            xylem_temp: TreeElement = TreeElement(pressure=pressure,
                                                  solutes=[],
                                                  viscosity=1e-3,
                                                  permeability=axial_permeability[0],
                                                  elastic_modulus=elastic_modulus[0],
                                                  hydraulic_conductivity=radial_hydraulic_conductivity,
                                                  height=self.height/self.num_elements,
                                                  radius=self.initial_radius[0],
                                                  transpiration_rate=transpiration_rate,
                                                  photosynthesis_rate=0,
                                                  sugar_loading_rate=0,
                                                  sugar_unloading_rate=0)

            phloem_temp: TreeElement = TreeElement(pressure=pressure,
                                                   solutes=[sugar],
                                                   viscosity=0,  # The viscosity is set below
                                                   permeability=axial_permeability[1],
                                                   elastic_modulus=elastic_modulus[1],
                                                   hydraulic_conductivity=radial_hydraulic_conductivity,
                                                   height=self.height/self.num_elements,
                                                   radius=self.initial_radius[1],
                                                   transpiration_rate=0,
                                                   photosynthesis_rate=photosynthesis_rate,
                                                   sugar_loading_rate=sugar_loading_rate,
                                                   sugar_unloading_rate=sugar_unloading_rate)
            self.elements.append([xylem_temp, phloem_temp])
        self.update_sap_viscosity()  # update viscosity

    def element_property_as_numpy_array(self, property: str) -> np.ndarray:
        xylem_values: np.ndarray = np.asarray([self.elements[i][0].__getattribute__(property)
                                               for i in range(self.num_elements)])
        phloem_values: np.npdarray = np.asarray([self.elements[i][1].__getattribute__(property)
                                                 for i in range(self.num_elements)])
        return np.transpose(np.concatenate([[xylem_values], [phloem_values]]))

    def sugar_concentration_as_numpy_array(self) -> np.ndarray:
        # TODO: make this search dynamic not "solutes[0]"
        return np.asarray([self.elements[i][1].solutes[0].concentration
                           for i in range(self.num_elements)])

    def element_area(self, ind: List[int] = [], column: int = 0) -> np.ndarray:
        """ returns element areas specified in parameter ind and column of self.elements.

            If no ind is given returns the areas for every element.
            If no column is given returns the areas in column 0 of self.elements (the xylem)
        """
        radii: np.ndarray = self.element_property_as_numpy_array('radius')

        if len(ind) > 0:
            radii = radii[ind, :]
        if column < MAX_ELEMENT_COLUMNS:
            radii = radii[:, 0:column+1]

        areas: np.ndarray = np.sum(radii[:, 0:column+1]**2, axis=1)  # every element column until desired column
        areas = areas - HEARTWOOD_RADIUS**2  # subtract the area of heartwood
        if column > 0:  # subtract the area of every element column before desired column
            for i in range(column):
                areas = areas - radii[:, i]**2
        return areas*math.pi

    def element_volume(self, ind: List[int] = [], column: int = 0) -> float:
        """ returns element volumes specified in parameter ind and column of self.elements.

            If no ind is given returns the volume of every element.
            If no column is given returns the volumes in column 0 of self.elements (the xylem)
        """
        # TODO: refactor the radii finding to own function (used multiple times)
        heights: np.ndarray = self.element_property_as_numpy_array('height')
        heights = heights[:, column]

        if len(ind) > 0:
            heights = heights[ind]

        return self.element_area(ind, column) * heights

    def update_sap_viscosity(self) -> None:
        """ Updates the viscosity in column 1 of self.elements (the phloem)
        """
        # TODO: refactor into solution class

        # calculate sugar volume in sap
        sugar_volume_fraction: np.ndarray = np.sum(self.sugar_concentration_as_numpy_array()) * M_SUCROSE / RHO_SUCROSE

        viscosity: np.ndarray = 1e-3 * np.exp(4.68 * 0.956 * sugar_volume_fraction /
                                              (1 - 0.956 * sugar_volume_fraction))
        # TODO: make a function to tree class that allows setting all the viscosity values at once
        for i in range(self.num_elements):
            self.elements[i][1].viscosity = viscosity
