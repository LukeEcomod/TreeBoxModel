from src.treeelement import TreeElement
from src.solute import Solute
from src.constants import RHO_WATER, GRAVITATIONAL_ACCELEREATION


class Tree:
    def __init__(self, height, num_elements, transpiration_profile,
                 photosynthesis_profile, sugar_profile,
                 sugar_loading_profile, sugar_unloading_profile,
                 axial_permeability_profile,
                 radial_hydraulic_conductivity_profile,
                 elastic_modulus_profile,
                 ground_water_potential):

        self.height = height  # unit: m

        self.num_elements = num_elements  # number of tree elements

        # (num_elements,1) array of transpiration rates in the xylem
        # unit: m3/s
        self.transpiration_profile = transpiration_profile

        # (num_elements,1) array of photosynth. rate in the phloem
        # unit: m3/s
        self.photosynthesis_profile = photosynthesis_profile

        # (num_elements,1) array of sugar conc. at t=0s in phloem
        # unit: mol/m3
        self.sugar_profile = sugar_profile

        # (num_elements,1) array of sugar loading rates in phloem
        # unit: m3/s
        self.sugar_loading_profile = sugar_loading_profile

        # (num_elements,1) array of sugar unloading rates in phloem
        # unit: m3/s
        self.sugar_unloading_profile = sugar_unloading_profile

        self.axial_permeability_profile = axial_permeability_profile

        self.radial_hydraulic_conductivity_profile = radial_hydraulic_conductivity_profile

        self.elastic_modulus_profile = elastic_modulus_profile

        self.ground_water_potential = ground_water_potential

        # calculate what the pressure in the xylem will be
        # assuming that at the base of the tree water potential
        # equals ground water potential
        # order 1 = top of the tree, N = base of the tree

        pressures = [ground_water_potential - i*RHO_WATER*GRAVITATIONAL_ACCELEREATION*self.height/num_elements for i in range(num_elements)]
        pressures = pressures.reverse()  # reverse so that the ordering is correct

        # set up dim = (num_elements, 2) list of TreeElements
        # columm 1 = xylem
        # column 2 = phloem

        self.elements = []
        for pressure, transpiration_rate, photosynthesis_rate, sugar_concentration,\
            sugar_loading_rate, sugar_unloading_rate in zip(
                pressures,
                transpiration_profile, photosynthesis_profile, sugar_profile,
                sugar_loading_profile, sugar_unloading_profile):

            sugar = Solute(0.3423, 1590, sugar_concentration)
            self.elements.append(sugar)  # NB! this is just for test
