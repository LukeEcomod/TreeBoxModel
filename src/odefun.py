import numpy as np
import math
from .constants import RHO_WATER, MAX_ELEMENT_COLUMNS


def odefun(t: float, y: np.ndarray, model) -> np.ndarray:
    """ Calculates the right hand side of the model ODEs.

    The modelled systen and the ODEs are described in the [modelled system](modelled_system.html).
    The scipy.solve_ivp() function in src.model.Model.run_scipy() method calls this function during
    the simulation.

    Args:
      t (float): time in the model simulation
      y (numpy.ndarray(dtype=float, ndims=1)[5 :math:`\\cdot` model.tree.num_elements,]): 1D array where elements
         0:2 :math:`\\cdot` model.tree.num_elements are the pressures in the xylem and phloem of the tree, elements
         2 :math:`\\cdot` model.tree.num_elements:3 :math:`\\cdot` model.tree.num_elements are for the sucrose
         concentration in the phloem and elements 3 :math:`\\cdot` model.tree.num_elements:5 :math:`\\cdot`
         model.tree.num_elements are for the element radii both in the xylem and the phloem.
      model (src.model.Model): Instace of the model class

    Returns:
      (numpy.ndarray(dtype=float,ndims=1)[5 :math:`\\cdot` model.tree.num_elements,]): 1D array of the right hand side
      values of the model ODEs where the elements 0:2 :math:`\\cdot` model.tree.num_elements are
      :math:`\\frac{\\text{d(pressure)}}{\\text{dt}}`, elements 2 :math:`\\cdot` model.tree.num_elements:
      3 :math:`\\cdot` model.tree.num_elements are :math:`\\frac{\\text{d[C(sucrose)]}}{\\text{dt}}` and elements
      3 :math:`\\cdot` model.tree.num_elements:5 :math:`\\cdot` model.tree.num_elements are
      :math:`\\frac{\\text{d(radius)}}{\\text{dt}}`

    """

    # Update model tree parameters
    pressures = y[0:model.tree.num_elements*2].reshape(model.tree.num_elements, MAX_ELEMENT_COLUMNS, order='F')
    sugar_concentration = y[model.tree.num_elements*2:model.tree.num_elements*3].reshape(model.tree.num_elements, 1,
                                                                                         order='F')
    element_radius = y[model.tree.num_elements*3:].reshape(model.tree.num_elements, MAX_ELEMENT_COLUMNS+1, order='F')

    model.tree.pressure = pressures
    model.tree.update_sugar_concentration(sugar_concentration)
    model.tree.element_radius = element_radius
    model.tree.update_sap_viscosity()

    model.tree.sugar_unloading_rate = model.tree.cross_sectional_area()\
        * np.max(np.concatenate([np.zeros((model.tree.num_elements, 1)), model.tree.sugar_unloading_slope
                                 * (model.tree.sugar_concentration_as_numpy_array()
                                    - model.tree.sugar_target_concentration)], axis=1), axis=1).reshape(
        model.tree.num_elements, 1)
    model.tree.sugar_unloading_rate[model.tree.root_elements] *= 10
    dmdt_ax: np.ndarray = model.axial_fluxes()[0]
    dmdt_rad: np.ndarray = model.radial_fluxes()

    dydt = [np.zeros((model.tree.num_elements, MAX_ELEMENT_COLUMNS)),
            np.zeros((model.tree.num_elements, 1)),
            np.zeros((model.tree.num_elements, MAX_ELEMENT_COLUMNS+1))]

    dydt[0] = model.tree.elastic_modulus/(np.transpose(
        np.array([model.tree.element_volume([], 0),
                  model.tree.element_volume([], 1)])) * RHO_WATER) * (dmdt_ax + dmdt_rad)

    # add flux from soil to root xylem (Pa/s)
    dydt[0][model.tree.root_elements] += model.tree.roots.conductivity(model.soil) *\
        (model.soil.pressures - model.tree.pressures[model.root_elements])

    dydt[1] = 1.0 / model.tree.element_volume([], 1).reshape(model.tree.num_elements, 1)\
        * (dmdt_ax[:, 1].reshape(model.tree.num_elements, 1)
           * model.tree.sugar_concentration_as_numpy_array()
           .reshape(model.tree.num_elements, 1)
           / RHO_WATER
           + model.tree.sugar_loading_rate.reshape(model.tree.num_elements, 1)
           - model.tree.sugar_unloading_rate.reshape(model.tree.num_elements, 1))

    dmdt = (dmdt_ax + dmdt_rad)
    dydt[2][:, 0] = 0.0
    dydt[2][:, 1] = ((dmdt[:, 0].reshape(model.tree.num_elements, 1))
                     / (2.0*math.pi*RHO_WATER * model.tree.element_height
                        * (model.tree.element_radius[:, 1].reshape(model.tree.num_elements, 1)
                           + model.tree.element_radius[:, 0].reshape(model.tree.num_elements, 1))
                        )).reshape(model.tree.num_elements,)

    dydt[2][:, 2] = ((dmdt[:, 1].reshape(model.tree.num_elements, 1))
                     / (2.0*math.pi*RHO_WATER*model.tree.element_height) *
                     ((model.tree.element_radius[:, 0].reshape(model.tree.num_elements, 1)
                       * model.tree.element_radius[:, 1].reshape(model.tree.num_elements, 1)
                       - model.tree.element_radius[:, 1].reshape(model.tree.num_elements, 1))
                      / ((model.tree.element_radius[:, 0].reshape(model.tree.num_elements, 1)
                          + model.tree.element_radius[:, 1].reshape(model.tree.num_elements, 1))
                         * (model.tree.element_radius[:, 0].reshape(model.tree.num_elements, 1)
                            + model.tree.element_radius[:, 1].reshape(model.tree.num_elements, 1)
                            + model.tree.element_radius[:, 2].reshape(model.tree.num_elements, 1))))
                     ).reshape(model.tree.num_elements,)
    # no radius change in the roots
    dydt[2][model.tree.root_elements, :] = 0.0
    dydt = np.concatenate([dydt[0].reshape(model.tree.num_elements*2, order='F'),
                           dydt[1].reshape(model.tree.num_elements, order='F'),
                           dydt[2].reshape(model.tree.num_elements*3, order='F')])
    return dydt
