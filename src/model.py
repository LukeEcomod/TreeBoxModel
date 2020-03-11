import numpy as np
import math
from .tree import Tree
from .constants import HEARTWOOD_RADIUS, XYLEM_RADIUS, PHLOEM_RADIUS


class Model:
    def __init__(self, tree: Tree):
        self.tree = tree

    def axial_fluxes(self) -> np.ndarray:
        """ calculates change in sap mass of every element

        variable names in this function's scope are the same as in Hölttä et al. (2006)
        shorter variable names are used to make errror prone calculations more readable
        """
        # TODO:figure out why ndarray type hint does not raise error in pyright
        # FIXME: refactor into axial and radial fluxes
        # TODO: think is this the smartest way to do the calculations
        pressures: np.ndarray = self.tree.element_property_as_numpy_array('pressure')
        L: np.ndarray = self.tree.element_property_as_numpy_array('hydraulic_conductivity')
        k: np.ndarray = self.tree.element_property_as_numpy_array('permeability')
        eta: np.ndarray = self.tree.element_property_as_numpy_array('viscosity')
        length: np.ndarray = self.tree.element_property_as_numpy_array('height')
        E: np.ndarray = self.tree.element_property_as_numpy_array('transpiration_rate')
        P: np.ndarray = self.tree.element_property_as_numpy_array('photosynthesis_rate')
        L: np.ndarray = self.tree.element_property_as_numpy_array('sugar_loading_rate')
        U: np.ndarray = self.tree.element_property_as_numpy_array('sugar_unloading_rate')
        C: np.ndarray = self.tree.sugar_concentration_as_numpy_array()

        # calculate transport coefficients
        transport_ax: np.ndarray = k/eta/length*np.transpose(
            np.array([self.tree.element_area([], 0),
                      self.tree.element_area([], 1)]))
        # calculate upward fluxes
        Q_ax_down: np.ndarray = np.zeros((self.tree.num_elements, pressures.shape[1]))
        Q_ax_down[0:-1, :] = np.diff(pressures, axis=0)*transport_ax[0:-1, :]
        Q_ax_down[-1, 0] = (self.tree.ground_water_potential-pressures[-1, 0])*transport_ax[-1, 0]
        Q_ax_down[-1, 1] = 0  # no flux from phloem to soil

        Q_ax_up: np.ndarray = np.zeros((self.tree.num_elements, pressures.shape[1]))
        # TODO: Think if there is a better way to achieve this without flipping twice
        Q_ax_up[1:len(Q_ax_up), :] = np.flip(
            np.diff(
                np.flip(
                    pressures, axis=0), axis=0), axis=0) * transport_ax[1:len(Q_ax_up), :]

        Q_ax_up[0, :] = 0  # the upward flux is handled in transpiration rate for the highest element

        Q_ax: np.ndarray = Q_ax_up + Q_ax_down - E
        return Q_ax
