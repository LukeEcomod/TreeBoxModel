from ..constants import M_SUCROSE, RHO_SUCROSE
from ..constants import VISCOSITY_WATER
import pytest
import math
import numpy as np
from ..tree import Tree
from ..roots import Roots
from ..soil import Soil
from .test_soil import test_soil
from .test_roots import test_roots
from typing import List


@pytest.fixture(scope="function")
def test_tree(test_roots, test_soil):

    # values from Hölttä et al. (2006)
    # except for transp/photosynth/loading/unloading/sugar profile
    height = 12.0
    element_height = np.concatenate((np.repeat(0.3, repeats=40),
                                     test_roots.layer_thickness(test_soil).reshape(5,)))
    num_elements = 45

    radii = [0.5, 1, 2]
    transpiration_profile: List[float] = [0 for i in range(num_elements)]
    transpiration_profile[0] = 0.9*1e-6  # m3/s

    photosynth_profile: List[float] = [0 for i in range(num_elements)]
    photosynth_profile[0:3] = [1e-6, 1e-6, 1e-6]

    sugar_profile = [10.0 for i in range(num_elements)]
    sugar_loading_profile = photosynth_profile

    sugar_unloading_profile = [0.0 for i in range(num_elements)]
    sugar_unloading_profile[-1] = 2

    axial_permeability_profile = [[1.5e-12, 6.0e-12]]*num_elements

    radial_hydr_conductivity = [1e-13] * num_elements

    elastic_modulus_profile = [[1000e6, 30e6]] * num_elements
    return Tree(height=height,
                num_elements=num_elements,
                element_height=element_height,
                initial_radius=radii,
                transpiration_profile=transpiration_profile,
                photosynthesis_profile=photosynth_profile,
                sugar_profile=sugar_profile,
                sugar_loading_profile=sugar_loading_profile,
                sugar_unloading_profile=sugar_unloading_profile,
                axial_permeability_profile=axial_permeability_profile,
                radial_hydraulic_conductivity_profile=radial_hydr_conductivity,
                elastic_modulus_profile=elastic_modulus_profile,
                sugar_target_concentration=0.9,
                sugar_unloading_slope=1,
                roots=test_roots)


def test_tree_init(test_tree):

    assert test_tree.height == 12

    # Test that pressure in xylem at ground level (Nth element)
    # is equal to ground water potential
    assert test_tree.pressure[-1, 0] == 0.0
    for i in range(test_tree.num_elements):
        assert test_tree.solutes[i, 1].concentration == 10

    assert np.array_equal(test_tree.root_elements, np.arange(40, 45))
    assert np.array_equal(test_tree.tree_elements, np.arange(0, 40))


def test_element_area(test_tree):

    assert test_tree.element_area([1], 0)[0] == math.pi*(1.5**2 - 0.5**2)
    assert test_tree.element_area([1], 1)[0] == pytest.approx(math.pi*(3.5**2 - 1.5**2), rel=1e-6)


def test_element_volume(test_tree):

    assert test_tree.element_volume([1], 0)[0] == pytest.approx(math.pi*(1.5**2 - 0.5**2)*12/40, rel=1e-6)
    assert test_tree.element_volume([1], 1)[0] == pytest.approx(math.pi*(3.5**2 - 1.5**2)*12/40, rel=1e-6)


def test_viscosity_update(test_tree):
    #  check that the viscosity in xylem is 1e-3 (water viscosity) in the start
    assert test_tree.viscosity[0, 0] == 1e-3
    # check that the viscosity in phloem is calculated correctly after initialization
    sugar_volume_fraction = test_tree.sugar_concentration_as_numpy_array() * M_SUCROSE / RHO_SUCROSE
    assert test_tree.viscosity[30, 1] == pytest.approx(VISCOSITY_WATER *
                                                       np.exp(4.68*0.956*sugar_volume_fraction[30] /
                                                              (1-0.956*sugar_volume_fraction[30])), rel=1e-6)
    # update sugar concentartion of all elements
    for i in range(test_tree.num_elements):
        test_tree.solutes[i, 1].concentration = 100
    test_tree.update_sap_viscosity()
    sugar_volume_fraction = test_tree.sugar_concentration_as_numpy_array() * M_SUCROSE / RHO_SUCROSE
    # check that the viscosity is calculated correctly after the concentration changes
    assert test_tree.viscosity[30, 1] == pytest.approx(
        VISCOSITY_WATER * np.exp(4.68*0.956*sugar_volume_fraction[30] /
                                 (1-0.956*sugar_volume_fraction[30])), rel=1e-6)


def test_cross_sectional_area(test_tree):
    cross_sectional_area = test_tree.cross_sectional_area()

    assert cross_sectional_area[0, 0] == pytest.approx(math.pi*2.0*(0.5+1)*0.3)
