import numpy as np
import pytest
from src.roots import Roots
from src.test.test_soil import test_soil


@pytest.fixture(scope="function")
def test_roots():
    num_elements = 5
    rooting_depth = 5

    area_density = np.zeros((5, 1)) + 1

    effective_radius = np.zeros((5, 1))+0.1

    soil_conductance_scale = 0.01

    RAI = 3

    total_root_area = np.sum(area_density * rooting_depth/num_elements)
    return Roots(rooting_depth=rooting_depth,
                 area_density=area_density,
                 effective_radius=effective_radius,
                 soil_conductance_scale=soil_conductance_scale,
                 total_root_area=total_root_area,
                 num_elements=num_elements,
                 RAI=RAI)


def test_init(test_roots):
    assert test_roots.soil_conductance_scale == 0.01
    assert all(test_roots.area_density == 1)


def test_root_conductance(test_roots, test_soil):
    result = test_roots.root_conductance(test_soil)
    assert all(result[:2] == 100)

    result_short = test_roots.root_conductance(test_soil, [1, 2])

    print(result_short)

    assert np.array_equal(result_short, np.asarray([100.0, 100.0]).reshape(2, 1))


def test_soil_root_conductance(test_roots, test_soil):
    result = test_roots.soil_root_conductance(test_soil)
    alpha = (test_roots.rooting_depth/test_roots.RAI)**(1/2)\
        * (2*test_roots.effective_radius)**(-1/2)
    assert result[0] == pytest.approx(alpha[0]*test_soil.hydraulic_conductivity[0]*test_roots.area_density[0], 1e-6)
    assert all(not ele for ele in np.isnan(result))


def test_conductivity(test_roots, test_soil):
    result = test_roots.conductivity(test_soil)
    ks = test_roots.soil_root_conductance(test_soil)
    kr = test_roots.root_conductance(test_soil)
    correct_result = (ks*kr)/(ks+kr).reshape(len(ks), 1)
    correct_result[np.isnan(correct_result)] = 0.0
    assert np.array_equal(result, correct_result)
