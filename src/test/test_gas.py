import numpy as np
import pytest
from src.gas import Gas
from src.test.source_sink_test import source


@pytest.fixture(scope="function")
def test_gas():
    nr = 5
    na = 10
    element_radius = np.repeat(1, 50).reshape(na, nr)
    element_height = np.repeat(0.1, 50).reshape(na, nr)
    D = np.repeat(1e-11, 50).reshape(na, nr)
    velocity = np.repeat(1e-4, 50).reshape(na, nr)
    air_fraction = np.repeat(0.35, 50).reshape(na, nr, 1)
    water_fraction = np.repeat(0.5, 50).reshape(na, nr, 1)
    cell_fraction = np.repeat(0.15, 50).reshape(na, nr, 1)
    space_division = np.concatenate((air_fraction, water_fraction, cell_fraction), axis=2)
    ambient_concentration = np.repeat(1, 10).reshape(na, 1)
    air_concentration = np.arange(0, 50).reshape(na, nr, 1)
    water_concentration = (0.342*air_concentration).reshape(na, nr, 1)
    concentration = np.concatenate((air_concentration, water_concentration), axis=2)
    axial_water_uptake = np.zeros((na, 1))
    axial_water_uptake[-1] = 1
    soil_gas_aq_concentration = np.zeros((na, 1))
    soil_gas_aq_concentration[-1] = 0.5
    sources_and_sinks = np.zeros((na, nr, 2))
    return Gas(num_radial_elements=nr,
               num_axial_elements=na,
               element_radius=element_radius,
               element_height=element_height,
               diffusion_coefficients=D,
               equilibration_rate=0.005,
               velocity=velocity,
               space_division=space_division,
               concentration=concentration,
               henrys_law_coefficient=0.342,
               temperature=298.15,
               ambient_concentration=ambient_concentration,
               axial_water_uptake=axial_water_uptake,
               soil_gas_aq_concentration=soil_gas_aq_concentration,
               sources_and_sinks=sources_and_sinks)


def test_init(test_gas):
    assert test_gas.temperature == 298.15
    assert test_gas.concentration.shape == (10, 5, 2)


def test_radius_from_pith(test_gas):
    r = test_gas.radius_from_pith
    assert all(a == b for a, b in zip(r[0, :], np.array([1, 2, 3, 4, 5])))


def test_head_area(test_gas):
    A = test_gas.head_area
    assert A[0, 3] == np.pi*(4**2-3**2)


def test_axial_fluxes(test_gas):
    Q = test_gas.axial_fluxes()

    # Test that the fluxes cancel out
    assert np.sum(np.sum(Q)) == pytest.approx(0.0, rel=1e-10)
    # since the concentration is highest at the bottom and lowest at the top
    # test that the first row (top of the tree) gets positive amount of concentration
    # and that the last row (bottom of the tree) loses concentration
    assert all(a > 0 for a in Q[0, :])
    assert all(a < 0 for a in Q[-1, :])


def test_radial_fluxes(test_gas):
    Q, _, _ = test_gas.radial_fluxes()
    cbark = np.arange(4, 50, 5)
    outflux = -2.0*np.pi*5*0.1*1e-11*(cbark-1)/0.5  # outflux = A_(stem_surface) * D *Delta(conc)/dr_atm (5-4.5)

    # Test that the sum flux equals flux out of the tree
    assert all(a == pytest.approx(b, rel=1e-15) for a, b in zip(np.sum(Q, axis=1), outflux))


def test_air_water_fluxes(test_gas):
    Q = test_gas.air_water_fluxes()

    assert np.sum(Q) == pytest.approx(0, rel=1e-15)

