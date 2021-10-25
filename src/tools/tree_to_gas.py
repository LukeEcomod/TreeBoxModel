from typing import Dict, Tuple
import numpy as np


def convert_tree_to_gas_properties(model, gas_dims: Tuple, c_gas_soil=0.0) -> Dict:
    '''Convert tree properties in model.py to gas properties to be used in gas.py'''

    r, r_mask = convert_tree_radii_to_gas(model, gas_dims)
    h, h_mask = convert_tree_height_to_gas(model, gas_dims)
    v = convert_tree_flux_to_velocity(model)
    n = gas_dims[0]*gas_dims[1]

    params = {}
    params['radius'] = np.repeat(r, repeats=n).reshape(gas_dims)
    params['height'] = np.repeat(h, repeats=n).reshape(gas_dims)
    params['velocity'] = np.zeros(gas_dims, dtype=np.float64)
    params['root_uptake'] = convert_root_fluxes_to_source_term(model, c_gas_soil)

    # set the velocity to correct values
    for row, _ in enumerate(params['velocity']):
        params['velocity'][row, :] = np.array([0.0 if i == 0
                                               else v[0][h_mask[row], 0] if i == 1
                                               else v[1][h_mask[row], 0] for i in r_mask])

    return params


def convert_tree_flux_to_velocity(model):
    '''Convert axial upward flux to sap flow velocity to be used in gas.py'''
    RHO_WATER = 1000
    _, flux, _ = model.axial_fluxes()
    flux = -1.0*flux/RHO_WATER  # flux in m3/s
    velocity_xylem = flux[:, 0].reshape(model.tree.num_elements, 1) / model.tree.element_area([], 0)
    velocity_phloem = flux[:, 1].reshape(model.tree.num_elements, 1) / model.tree.element_area([], 1)

    return (velocity_xylem, velocity_phloem)


def convert_root_fluxes_to_source_term(model, c_gas_soil: float):
    ''' Convert root fluxes to be used as a source term in sources_and_sinks_func in gas.py'''
    RHO_WATER = 1000
    return model.root_fluxes()/RHO_WATER*c_gas_soil


def convert_tree_radii_to_gas(model, gas_dims: Tuple) -> Tuple:
    """Converts the tree heartwood, xylem and phloem radii to equally spaced gas element radii

    Args:
        tree (Tree): Instance of the tree class
        gas_dims (Tuple): Number of axial and radial elements in the gas in this order.

    Returns:
        Tuple where the first element is the radii for every element and the second list
        contains mask whether the element is heartwood (0), sapwood (1) or phloem(2)
    """

    r = (model.tree.element_radius[0, :])
    percentage_r = np.cumsum(r/np.sum(r))
    r = np.sum(r)/gas_dims[1]
    mask = [0 if i/gas_dims[1] < percentage_r[0]
            else 1 if (i/gas_dims[1] < percentage_r[1] and i/gas_dims[1] >= percentage_r[0])
            else 2 for i in range(gas_dims[1])]
    return (r, mask)


def convert_tree_height_to_gas(model, gas_dims: Tuple) -> Tuple:
    """ Converts tree element heights to equally spaced heights such that
        .. math::
            n*h_{new} = \\sum_{i=1}^m h_i

        where
        * :math:`n`: gas_dims[0] (number of axial elements in gas).
        * :math:`m`: number of axial elements in the tree.
        * :math:`h_{new}`: new element height.
        * :math:`h_i`: element heights in the tree.

    Args:
        tree (Tree): Instance of the tree class
        gas_dims (Tuple): Number of axial and radial elements in the gas in this order.
    """
    h = np.sum(model.tree.element_height)
    h_percentage = np.cumsum(model.tree.element_height/h)
    mask = np.zeros((gas_dims[0], 1), dtype='int')
    gas_h_percentage = np.array([i/gas_dims[0] for i in range(gas_dims[0])])
    gas_h_percentage[-1] = 1.0
    for (ind, _) in enumerate(h_percentage):
        if ind == 0:
            mask[np.where(gas_h_percentage <= h_percentage[0])] = ind
        else:
            mask[np.where((gas_h_percentage <= h_percentage[ind]) & (gas_h_percentage > h_percentage[ind-1]))] = ind
    h = h/gas_dims[0]
    return (h, mask)
