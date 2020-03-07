import numpy as np
import math
from .tree import Tree
from .constants import HEARTWOOD_RADIUS, XYLEM_RADIUS, PHLOEM_RADIUS


def calculate_axial_fluxes(tree: Tree) -> np.ndarray:
    """ calculates change in sap mass of every element

    variable names in this function's scope are the same as in Hölttä et al. (2006)
    shorter variable names are used to make errror prone calculations more readable
    """
    # TODO:figure out why ndarray type hint does not raise error in pyright
    # FIXME: refactor into axial and radial fluxes
    # TODO: think is this the smartest way to do the calculations
    pressures = tree.element_property_as_numpy_array('pressure')
    L = tree.element_property_as_numpy_array('hydraulic_conductivity')
    k = tree.element_property_as_numpy_array('permeability')
    eta = tree.element_property_as_numpy_array('viscosity')
    length = tree.element_property_as_numpy_array('height')
    E = tree.element_property_as_numpy_array('transpiration_rate')
    P = tree.element_property_as_numpy_array('photosynthesis_rate')
    L = tree.element_property_as_numpy_array('sugar_loading_rate')
    U = tree.element_property_as_numpy_array('sugar_unloading_rate')
    C = tree.sugar_concentration_as_numpy_array()

    xylem_area = math.pi*((HEARTWOOD_RADIUS + XYLEM_RADIUS)**2 - HEARTWOOD_RADIUS**2)
    phloem_area = math.pi*((HEARTWOOD_RADIUS + XYLEM_RADIUS + PHLOEM_RADIUS)**2
                           - (HEARTWOOD_RADIUS+XYLEM_RADIUS)**2)
    # calculate transport coefficients
    transport_ax = k/eta/length*xylem_area
    # calculate fluxes
    Q_ax = np.diff(pressures, axis=0)*transport_ax[0:-1, :]
    Q_rad = 0

    # update phloem sap viscosity
    tree.update_sap_viscosity()

    return 4
