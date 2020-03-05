from typing import List
from .treeelement import TreeElement
from .solute import Solute
from .constants import RHO_WATER, GRAVITATIONAL_ACCELERATION


class Tree:
    def __init__(self, height, num_elements, transpiration_profile,
                 photosynthesis_profile, sugar_profile,
                 sugar_loading_profile, sugar_unloading_profile,
                 axial_permeability_profile,
                 radial_hydraulic_conductivity_profile,
                 elastic_modulus_profile,
                 ground_water_potential):

        self.height: float = height  # unit: m

        self.num_elements: int = num_elements  # number of tree elements

        # (num_elements,1) array of transpiration rates in the xylem
        # unit: m3/s
        self.transpiration_profile: List[float] = transpiration_profile

        # (num_elements,1) array of photosynth. rate in the phloem
        # unit: m3/s
        self.photosynthesis_profile: List[float] = photosynthesis_profile

        # (num_elements,1) array of sugar conc. at t=0s in phloem
        # unit: mol/m3
        self.sugar_profile: List[float] = sugar_profile

        # (num_elements,1) array of sugar loading rates in phloem
        # unit: m3/s
        self.sugar_loading_profile: List[float] = sugar_loading_profile

        # (num_elements,1) array of sugar unloading rates in phloem
        # unit: m3/s
        self.sugar_unloading_profile: List[float] = sugar_unloading_profile

        self.axial_permeability_profile: List[float] = axial_permeability_profile

        self.radial_hydraulic_conductivity_profile: List[float] = radial_hydraulic_conductivity_profile

        self.elastic_modulus_profile: List[float] = elastic_modulus_profile

        self.ground_water_potential: float = ground_water_potential

        # calculate what the pressure in the xylem will be
        # assuming that at the base of the tree water potential
        # equals ground water potential
        # order 1 = top of the tree, N = base of the tree

        pressures = [ground_water_potential - i*RHO_WATER*GRAVITATIONAL_ACCELERATION*self.height/self.num_elements
                     for i in range(num_elements)]
        pressures.reverse()  # reverse so that the ordering is correct
        # set up dim = (num_elements, 2) list of TreeElements
        # columm 1 = xylem
        # column 2 = phloem
        self.elements: List[List[TreeElement]] = []
        i = 1
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

            sugar: Solute = Solute(0.3423, 1590, sugar_concentration)
            xylem_temp: TreeElement = TreeElement(pressure=pressure,
                                                  solutes=[],
                                                  viscosity=1e-3,
                                                  permeability=axial_permeability[0],
                                                  elastic_modulus=elastic_modulus[0],
                                                  hydraulic_conductivity=radial_hydraulic_conductivity,
                                                  height=self.height/self.num_elements,
                                                  transpiration_rate=transpiration_rate,
                                                  photosynthesis_rate=photosynthesis_rate,
                                                  sugar_loading_rate=0,
                                                  sugar_unloading_rate=0)

            phloem_temp: TreeElement = TreeElement(pressure=pressure,
                                                   solutes=[sugar],
                                                   viscosity=1e-3,  # TODO: calculate viscosity dynamically
                                                   permeability=axial_permeability[1],
                                                   elastic_modulus=elastic_modulus[1],
                                                   hydraulic_conductivity=radial_hydraulic_conductivity,
                                                   height=self.height/self.num_elements,
                                                   transpiration_rate=transpiration_rate,
                                                   photosynthesis_rate=photosynthesis_rate,
                                                   sugar_loading_rate=0,
                                                   sugar_unloading_rate=0)
            self.elements.append([xylem_temp, phloem_temp])
