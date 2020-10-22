import numpy as np
from ..soil import Soil
import pytest


@pytest.fixture(scope="function")
def test_soil():
    layer_thickness = np.diff(np.linspace(0, 100, 101))
    hydraulic_conductivity = np.repeat([2], repeats=100)
    pressure = np.linspace(-2, -0.02, 100)

    return Soil(layer_thickness=layer_thickness,
                hydraulic_conductivity=hydraulic_conductivity,
                pressure=pressure)


def test_init(test_soil):
    assert test_soil.num_elements == 100


def test_depth(test_soil):
    result = test_soil.depth()
    assert np.array_equal(result[0:4], np.array([0.5, 1.5, 2.5, 3.5]).reshape(4, 1))
