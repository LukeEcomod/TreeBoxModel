import numpy as np
import math
from src.constants import RHO_WATER, MAX_ELEMENT_COLUMNS


def odefun(t, y, model):
    # Update model tree parameters
    pressures = y[0:model.tree.num_elements*2].reshape(model.tree.num_elements, MAX_ELEMENT_COLUMNS, order='F')
    sugar_concentration = y[model.tree.num_elements*2:model.tree.num_elements*3].reshape(model.tree.num_elements, 1,
                                                                                         order='F')
    element_radius = y[model.tree.num_elements*3:].reshape(model.tree.num_elements, MAX_ELEMENT_COLUMNS, order='F')

    model.tree.pressure = pressures
    model.tree.update_sugar_concentration(sugar_concentration)
    model.tree.element_radius = element_radius
    model.tree.update_sap_viscosity()
    # model.tree.sugar_unloading_rate[-5:] = (sugar_concentration[-5:] - model.tree.sugar_target_concentration)*1e-9

    dmdt_ax: np.ndarray = model.axial_fluxes()
    dmdt_rad: np.ndarray = model.radial_fluxes()

    dydt = [np.zeros((model.tree.num_elements, MAX_ELEMENT_COLUMNS)),
            np.zeros((model.tree.num_elements, 1)),
            np.zeros((model.tree.num_elements, MAX_ELEMENT_COLUMNS))]

    dydt[0] = model.tree.elastic_modulus/(np.transpose(
        np.array([model.tree.element_volume([], 0),
                  model.tree.element_volume([], 1)])) * RHO_WATER) * (dmdt_ax + dmdt_rad)
    dydt[1] = 1.0 / model.tree.element_volume([], 1).reshape(40, 1) * (dmdt_ax[:, 1].reshape(40, 1)
                                                                       * model.tree.sugar_concentration_as_numpy_array()
                                                                       .reshape(40, 1)
                                                                       / RHO_WATER
                                                                       + model.tree.sugar_loading_rate.reshape(40, 1)
                                                                       - model.tree.sugar_unloading_rate.reshape(40, 1))

    dydt[2] = (dmdt_ax+dmdt_rad) / (math.pi*RHO_WATER * np.repeat(model.tree.element_height, 2, axis=1)
                                    * model.tree.element_radius)
    dydt = np.concatenate([dydt[0].reshape(model.tree.num_elements*2, order='F'),
                           dydt[1].reshape(model.tree.num_elements, order='F'),
                           dydt[2].reshape(model.tree.num_elements*2, order='F')])
    return dydt
