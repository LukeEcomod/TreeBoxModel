import numpy as np
import math
from src.constants import HEARTWOOD_RADIUS, RHO_WATER, MAX_ELEMENT_COLUMNS


def odefun(t, y, model):
    """ Calculates the right hand side of the model ODEs."""

    # Update model tree parameters
    pressures = y[0:model.tree.num_elements*2].reshape(model.tree.num_elements, MAX_ELEMENT_COLUMNS, order='F')
    sugar_concentration = y[model.tree.num_elements*2:model.tree.num_elements*3].reshape(model.tree.num_elements, 1,
                                                                                         order='F')
    element_radius = y[model.tree.num_elements*3:].reshape(model.tree.num_elements, MAX_ELEMENT_COLUMNS, order='F')

    model.tree.pressure = pressures
    model.tree.update_sugar_concentration(sugar_concentration)
    model.tree.element_radius = element_radius
    model.tree.update_sap_viscosity()

    model.tree.sugar_unloading_rate = model.tree.cross_sectional_area()\
        * np.max(np.concatenate([np.zeros((model.tree.num_elements, 1)), model.tree.sugar_unloading_slope
                                 * (model.tree.sugar_concentration_as_numpy_array()
                                    - model.tree.sugar_target_concentration)], axis=1), axis=1).reshape(
        model.tree.num_elements, 1)
    model.tree.sugar_unloading_rate[-1] *= 10
    model.tree.sugar_unloading_rate[0:20] = np.zeros((20, 1))
    dmdt_ax: np.ndarray = model.axial_fluxes()
    dmdt_rad: np.ndarray = model.radial_fluxes()

    dydt = [np.zeros((model.tree.num_elements, MAX_ELEMENT_COLUMNS)),
            np.zeros((model.tree.num_elements, 1)),
            np.zeros((model.tree.num_elements, MAX_ELEMENT_COLUMNS))]

    dydt[0] = model.tree.elastic_modulus/(np.transpose(
        np.array([model.tree.element_volume([], 0),
                  model.tree.element_volume([], 1)])) * RHO_WATER) * (dmdt_ax + dmdt_rad)

    dydt[1] = 1.0 / model.tree.element_volume([], 1).reshape(model.tree.num_elements, 1)\
        * (dmdt_ax[:, 1].reshape(model.tree.num_elements, 1)
           * model.tree.sugar_concentration_as_numpy_array()
           .reshape(model.tree.num_elements, 1)
           / RHO_WATER
           + model.tree.sugar_loading_rate.reshape(model.tree.num_elements, 1)
           - model.tree.sugar_unloading_rate.reshape(model.tree.num_elements, 1))

    dmdt = (dmdt_ax + dmdt_rad)
    dydt[2][:, 0] = ((dmdt[:, 0].reshape(model.tree.num_elements, 1))
                     / (2.0*math.pi*RHO_WATER * model.tree.element_height
                        * (model.tree.element_radius[:, 0].reshape(model.tree.num_elements, 1)
                           + HEARTWOOD_RADIUS))).reshape(model.tree.num_elements,)

    dydt[2][:, 1] = ((dmdt[:, 1].reshape(model.tree.num_elements, 1))
                     / (2.0*math.pi*RHO_WATER*model.tree.element_height) *
                     ((HEARTWOOD_RADIUS*model.tree.element_radius[:, 0].reshape(model.tree.num_elements, 1)
                       - model.tree.element_radius[:, 0].reshape(model.tree.num_elements, 1))
                      / ((HEARTWOOD_RADIUS+model.tree.element_radius[:, 0].reshape(model.tree.num_elements, 1))
                         * (HEARTWOOD_RADIUS+model.tree.element_radius[:, 0].reshape(model.tree.num_elements, 1)
                            + model.tree.element_radius[:, 1].reshape(model.tree.num_elements, 1))))
                     ).reshape(model.tree.num_elements,)

    dydt = np.concatenate([dydt[0].reshape(model.tree.num_elements*2, order='F'),
                           dydt[1].reshape(model.tree.num_elements, order='F'),
                           dydt[2].reshape(model.tree.num_elements*2, order='F')])
    return dydt
