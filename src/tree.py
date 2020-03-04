from src.treeelement import TreeElement
from src.solute import Solute


class Tree:
    def __init__(self, height, num_elements, transpiration_profile,
                 photosynthesis_profile, sugar_profile,
                 sugar_loading_profile, sugar_unloading_profile,
                 axial_permeability_profile,
                 radial_hydraulic_conductivity_profile,
                 elastic_modulus_profile):

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
        # set up dim = (num_elements, 2) list of TreeElements
        # columm 1 = xylem
        # column 2 = phloem
        # order 1 = top of the tree, N = base of the tree
        self.elements = []
        for transpiration_rate, photosynthesis_rate, sugar_concentration,\
            sugar_loading_rate, sugar_unloading_rate in zip(
                 transpiration_profile, photosynthesis_profile, sugar_profile,
                 sugar_loading_profile, sugar_unloading_profile):

            sugar = Solute(0.3423, 1590, sugar_concentration)
            self.elements.append(sugar)
