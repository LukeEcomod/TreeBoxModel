import numpy as np
from ..gas import Gas
import pytest


@pytest.fixture(scope="function")
def test_gas():
    nr = 5
    na = 10
    element_radius = np.repeat(1, 50).reshape(na, nr)
    element_height = np.repeat(0.1, 50).reshape(na, nr)
    D = np.repeat(1e-11, 50).reshape(na, nr)
    eq_rate = 0.005
    velocity = np.repeat(1e-4, 50).reshape(na, nr)
    air_fraction = np.repeat(0.35, 50).reshape(na, nr, 1)
    water_fraction = np.repeat(0.5, 50).reshape(na, nr, 1)
    cell_fraction = np.repeat(0.15, 50).reshape(na, nr, 1)
    space_division = np.concatenate((air_fraction, water_fraction, cell_fraction), axis=2)
    kh = 0.342
    temperature = 298.15
    ambient_concentration = 1
    air_concentration = np.repeat(ambient_concentration, 50).reshape(na, nr, 1)
    water_concentration = np.repeat(kh*ambient_concentration, 50).reshape(na, nr, 1)
    concentration = np.concatenate((air_concentration, water_concentration), axis=2)

    return Gas(num_radial_elements=nr,
               num_axial_elements=na,
               element_radius=element_radius,
               element_height=element_height,
               diffusion_coefficients=D,
               equilibration_rate=eq_rate,
               velocity=velocity,
               space_division=space_division,
               concentration=concentration,
               henrys_law_coefficient=kh,
               temperature=temperature,
               ambient_concentration=ambient_concentration)


def test_init(test_gas):
    assert test_gas.temperature == 298.15


def test_radius_from_pith(test_gas):
    assert 100 == 100


def test_head_area(test_gas):
    assert 100 == 100
