import numpy as np
from ..soil import Soil
import pytest


@pytest.fixture(scope="function")
def test_soil():
    depth = np.linspace(0, 100, 101)
    hydraulic_conductivity = np.repeat([2], repeats=100)
    pressure = np.linspace(-2, -0.02, 100)

    return Soil(depth=depth, hydraulic_conductivity=hydraulic_conductivity, pressure=pressure)


def test_init(test_soil):
    assert test_soil.num_elements == 100


def test_layer_thickness(test_soil):
    assert all(test_soil.layer_thickness() == 1)
