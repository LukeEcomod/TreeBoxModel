import numpy as np
import pytest
from ..soil import Soil
from ..roots import Roots
from .test_soil import test_soil


@pytest.fixture(scope="function")
def test_roots():
    num_elements = 5
    rooting_depth = 5

    area_density = np.zeros((5, 1))
    area_density[:2] = 10

    effective_radius = np.zeros((5, 1))
    effective_radius[:2] = 0.1

    soil_conductance_scale = 0.01
    area_per_tree = 1
    return Roots(rooting_depth=rooting_depth,
                 area_density=area_density,
                 effective_radius=effective_radius,
                 soil_conductance_scale=soil_conductance_scale,
                 area_per_tree=area_per_tree,
                 num_elements=num_elements)


def test_init(test_roots):
    assert test_roots.soil_conductance_scale == 0.01
    assert any(test_roots.area_density == 10)


def test_root_conductance(test_roots, test_soil):
    result = test_roots.root_conductance(test_soil)
    assert all(result[:2] == 1000)

    result_short = test_roots.root_conductance(test_soil, [1, 2])

    assert np.array_equal(result_short, np.asarray([1000.0, 0.0]).reshape(2, 1))


def test_soil_root_conductance(test_roots, test_soil):
    result = test_roots.soil_root_conductance(test_soil)
    alpha = (test_roots.rooting_depth/test_roots.root_area_index(test_soil))**(1/2)\
        * (2*test_roots.effective_radius)**(-1/2)
    assert result[0] == alpha[0]*test_soil.hydraulic_conductivity[0]*test_roots.area_density[0]
    assert all([not ele for ele in np.isnan(result)])


def test_conductivity(test_roots, test_soil):
    result = test_roots.conductivity(test_soil)
    ks = test_roots.soil_root_conductance(test_soil)
    kr = test_roots.root_conductance(test_soil)
    print(result)
    correct_result = (ks*kr)/(ks+kr).reshape(len(ks), 1)
    correct_result[np.isnan(correct_result)] = 0.0
    assert np.array_equal(result, correct_result)
