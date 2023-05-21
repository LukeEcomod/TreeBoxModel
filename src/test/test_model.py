from ..constants import GRAVITATIONAL_ACCELERATION, MOLAR_GAS_CONSTANT, RHO_WATER
import pytest
import numpy as np
from src.model import Model
from .test_tree import test_tree
from .test_soil import test_soil
from. test_roots import test_roots


@pytest.fixture(scope="function")
def test_model(test_tree, test_soil):
    # change all viscosities to be 1
    axial_permeability = [10, 20]
    for i in range(test_tree.num_elements):
        for j in range(2):
            test_tree.viscosity[i, j] = 1
            test_tree.axial_permeability[i, j] = axial_permeability[j]
            test_tree.element_height[i] = 1
            test_tree.pressure[i, j] = i
            test_tree.transpiration_rate[i] = 1
    return Model(test_tree, test_soil)


def test_model_init(test_model):
    assert test_model.tree.pressure[test_model.tree.tree_elements[-1], 0] == 39
    assert test_model.tree.viscosity[5, 1] == 1
    assert test_model.tree.axial_permeability[28, 1] == 20


def test_axial_fluxes(test_model):
    flux, flux_up, flux_down = test_model.axial_fluxes()
    # length = test_model.tree.element_height
    # l= np.concatenate(([0], length.reshape(len(length),)))
    # cumulative_sum = np.cumsum(l).reshape(len(l), 1)
    # l_midpoints = length/2 + cumulative_sum[:-1]
    # dl_midpoints = np.diff(l_midpoints, axis=0)
    # print(dl_midpoints)
    #print(test_model.tree.element_height)
    transport_ax = np.asarray([10, 20]*test_model.tree.num_elements).reshape(
        test_model.tree.num_elements, 2)/1/test_model.tree.element_height*RHO_WATER*np.concatenate(
            [test_model.tree.element_area([], 0),
             test_model.tree.element_area([], 1)], axis=1)

    # test flux up
    for ind, f in enumerate(flux_up):
        if ind == 0:
            assert f[0] == 0
        else:

            assert f[0] == transport_ax[ind, 0]*(test_model.tree.pressure[ind-1, 0]-test_model.tree.pressure[ind, 0]
                                                 + RHO_WATER*GRAVITATIONAL_ACCELERATION
                                                 * test_model.tree.element_height[ind, 0])
    # test flux down

    for ind, f in enumerate(flux_down):
        if ind == flux_down.shape[0]-1:
            assert f[1] == 0
            assert f[0] == 0

        else:
            assert f[0] == transport_ax[ind, 0]*(test_model.tree.pressure[ind+1, 0]-test_model.tree.pressure[ind, 0]
                                                 - RHO_WATER*GRAVITATIONAL_ACCELERATION
                                                 * test_model.tree.element_height[ind, 0])

    # test full flux now that up and down are correct
    Qroot = test_model.root_fluxes()
    for ind, f in enumerate(flux):
        assert f[0] == flux_down[ind, 0] + flux_up[ind, 0] - test_model.tree.transpiration_rate[ind, 0] + Qroot[ind, 0]


def test_radial_fluxes(test_model):

    flux = test_model.radial_fluxes()
    N = test_model.tree.num_elements
    C = test_model.tree.sugar_concentration_as_numpy_array()
    transport_rad = test_model.tree.radial_hydraulic_conductivity.reshape(N, 1)\
        * test_model.tree.element_height.reshape(N, 1)\
        * test_model.tree.cross_sectional_area()*RHO_WATER

    for ind, f in enumerate(flux):

        assert f[0] == pytest.approx(transport_rad[ind, 0]
                                     * ((test_model.tree.pressure[ind, 1] - test_model.tree.pressure[ind, 0])
                                        - C[ind, 0]*MOLAR_GAS_CONSTANT*298.15), rel=1e-6)

        assert f[1] == pytest.approx(-f[0], rel=1e-6)

def test_mass_balance(test_model):
    ''' With these values there should be no storage change in the tree i.e., roots should take everything that is lost by transpiration '''
    # advance model 1 hour to let it equilibrate
    test_model.run_scipy(time_start=1e-10, time_end=3600)

    _, Q_ax_down, Q_ax_up = test_model.axial_fluxes()

    Q_ax = Q_ax_up + Q_ax_down

    assert np.sum(Q_ax) == pytest.approx(0, 1e-5)
