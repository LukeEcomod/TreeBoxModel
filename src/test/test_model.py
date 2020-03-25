from ..constants import GRAVITATIONAL_ACCELERATION, HEARTWOOD_RADIUS, RHO_WATER
import pytest
import math
import numpy as np
from ..tree import Tree
from ..model import Model
from .test_tree import test_tree
from typing import List


@pytest.fixture(scope="function")
def test_model(test_tree):
    # change all viscosities to be 1
    axial_permeability = [10, 20]
    transpiration = [1, 0]
    for i in range(test_tree.num_elements):
        for j in range(2):
            test_tree.elements[i][j].viscosity = 1
            test_tree.elements[i][j].permeability = axial_permeability[j]
            test_tree.elements[i][j].height = 1
            test_tree.elements[i][j].pressure = i
            test_tree.elements[i][j].transpiration_rate = transpiration[j]
    return Model(test_tree)


def test_model_init(test_model):
    assert test_model.tree.element_property_as_numpy_array('pressure')[-1][0] == 39
    assert test_model.tree.element_property_as_numpy_array('viscosity')[5][1] == 1
    assert test_model.tree.element_property_as_numpy_array('permeability')[28][1] == 20


def test_axial_fluxes(test_model):
    for i in range(test_model.tree.num_elements):
        assert test_model.axial_fluxes()[i][0] == pytest.approx(10*math.pi*(1-HEARTWOOD_RADIUS**2)-1, rel=1e6)

    assert test_model.axial_fluxes()[-1][1] == pytest.approx(10*math.pi*(1-HEARTWOOD_RADIUS**2)*(-39)-1, rel=1e6)

    # The phloem axial fluxes should sum to zero
    sumflux: np.ndarray = np.sum(test_model.axial_fluxes(), axis=0)
    assert sumflux[1] == 0

    # The xylem axial fluxes should be equal to the out/ingoing flux between the last element and soil
    # and the sum of transpiration
    lastpressure: float = test_model.tree.element_property_as_numpy_array('pressure')[-1][0]
    transport_coeff: float = 10*math.pi*(1-HEARTWOOD_RADIUS**2)
    assert sumflux[0] == pytest.approx(
        (test_model.tree.ground_water_potential - lastpressure)*transport_coeff - 40, rel=1e-6)