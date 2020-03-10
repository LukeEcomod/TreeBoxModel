from ..constants import HEARTWOOD_RADIUS, PHLOEM_RADIUS, XYLEM_RADIUS
import pytest
import math
from ..tree import Tree
from typing import List


@pytest.fixture(scope="function")
def test_tree():

    # values from Hölttä et al. (2006)
    # except for transp/photosynth/loading/unloading profile
    height = 12.0

    num_elements = 40

    radii = [1, 2]
    transpiration_profile: List[float] = [0 for i in range(num_elements)]
    transpiration_profile[0] = 0.9*1e-6  # m3/s

    photosynth_profile: List[float] = [0 for i in range(num_elements)]
    photosynth_profile[0:3] = [1e-6, 1e-6, 1e-6]

    sugar_profile = [1e1 for i in range(num_elements)]
    sugar_loading_profile = photosynth_profile

    sugar_unloading_profile = [0 for i in range(num_elements)]
    sugar_unloading_profile[-1] = 2

    axial_permeability_profile = [[1.5e-12, 6.0e-12]]*num_elements

    radial_hydr_conductivity = [1e-13] * num_elements

    elastic_modulus_profile = [[1000e6, 30e6]] * num_elements

    ground_water_potential = 0.0

    return Tree(height=height,
                num_elements=num_elements,
                initial_radius=radii,
                transpiration_profile=transpiration_profile,
                photosynthesis_profile=photosynth_profile,
                sugar_profile=sugar_profile,
                sugar_loading_profile=sugar_loading_profile,
                sugar_unloading_profile=sugar_unloading_profile,
                axial_permeability_profile=axial_permeability_profile,
                radial_hydraulic_conductivity_profile=radial_hydr_conductivity,
                elastic_modulus_profile=elastic_modulus_profile,
                ground_water_potential=ground_water_potential)


def test_tree_init(test_tree):

    assert test_tree.height == 12

    # Test that pressure at ground level (Nth element)
    # is equal to ground water potential
    assert test_tree.elements[-1][0].pressure == 0.0
    for i in range(test_tree.num_elements):
        assert test_tree.elements[i][1].solutes[0].concentration == test_tree.sugar_profile[i]
    # TODO: write more tests in initialization


def test_element_area(test_tree):

    assert test_tree.element_area([1], 0)[0] == math.pi*(1 - HEARTWOOD_RADIUS**2)
    assert test_tree.element_area([1], 1)[0] == pytest.approx(math.pi*(2**2 - HEARTWOOD_RADIUS**2), rel=1e-6)


def test_element_volume(test_tree):

    assert test_tree.element_volume([1], 0)[0] == pytest.approx(math.pi*(1 - HEARTWOOD_RADIUS**2)*12/40, rel=1e-6)
    assert test_tree.element_volume([1], 1)[0] == pytest.approx(math.pi*(2**2 - HEARTWOOD_RADIUS**2)*12/40, rel=1e-6)
