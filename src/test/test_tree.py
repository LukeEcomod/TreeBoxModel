import pytest

from src.tree import Tree

# values from Hölttä et al. (2006)
# except for transp/photosynth/loading/unloading profile
height = 12

num_elements = 40

transpiration_profile = [0 for i in range(num_elements)]
transpiration_profile[0] = 1

photosynth_profile = [0 for i in range(num_elements)]
photosynth_profile[0:3] = [1, 1, 1]

sugar_profile = [0 for i in range(0, 10)]
sugar_loading_profile = photosynth_profile

sugar_unloading_profile = [0, 0, 0, 0, 0, 0, 0, 0, 0, 2]

axial_permeability_profile = [[1.5e-12, 6.0e-12]] * num_elements

radial_hydr_conductivity = [1e-13] * num_elements

elastic_modulus_profile = [[1000e6, 30e6]] * num_elements

test_tree = Tree(height=height,
                 num_elements=num_elements,
                 transpiration_profile=transpiration_profile,
                 photosynthesis_profile=photosynth_profile,
                 sugar_profile=sugar_profile,
                 sugar_loading_profile=sugar_loading_profile,
                 sugar_unloading_profile=sugar_unloading_profile,
                 axial_permeability_profile=axial_permeability_profile,
                 radial_hydraulic_conductivity_profile=radial_hydr_conductivity,
                 elastic_modulus_profile=elastic_modulus_profile)


def test_tree_init():

    assert test_tree.height == 12
    assert test_tree.elements[1].concentration == 0
