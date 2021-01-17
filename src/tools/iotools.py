import numpy as np
import os.path
from ..tree import Tree
from ..gas import Gas
from typing import Dict
from netCDF4 import Dataset, Variable
# from ..model import Model # TODO: check why this import fails
from ..model_variables import (index_dim_vars, soil_dim_vars, root_dim_vars, axial_layer_dim_vars,
                               gas_three_dim_vars, gas_index_dim_vars, gas_axial_dim_vars)


def initialize_netcdf(filename: str,
                      dimensions: Dict,
                      variables: Dict) -> Dataset:
    """ Initializes a netcdf file to be ready for saving simulation results.

        Args:
            filename (str): name of the NETCDF4 file that is created
            dimensions (Dict): dimensions of the ncf file to be created.
                Key = name of the dimension, value=dimension length.
            variables (Dict): The names, descriptions, units, dimensions and precision of each variable that is saved
                to the netcdf file. The key of each dictionary element is used to label the variables.
                The value of each dictionary element needs to be a list where <br />
                list[0] (str) Description of the variable <br />
                list[1] (str): unit of the variable  <br />
                list[2] (Tuple): NETCDF dimensions of the variable  <br />
                list[3] (str): Datatype (precision) of the variable  <br />
                The possible NETCDF dimensions are "index", "radial_layers" and "axial_layers".

        Returns:
            (NETCDF4.Dataset): A NETCDF4 file where the simulation results can be saved.

        """
    # create the file if it does not exist
    if os.path.isfile(filename):
        ncf = Dataset(filename, "a")
    else:
        ncf = Dataset(filename, "w", format="NETCDF4")

    # create dimensions
    ncf.createDimension("index", 0)
    for (key, value) in dimensions.items():
        ncf.createDimension(key, value)

    # create the index variable
    ind: Variable = ncf.createVariable('index', 'i4', ('index'))
    ind.description = 'NETCDF index variable'
    # create the other variables
    for var_name in variables.keys():
        description: str = variables[var_name][0]
        units: str = variables[var_name][1]
        dimensions: str = variables[var_name][2]
        datatype: str = variables[var_name][3]
        new_var: Variable = ncf.createVariable(var_name, datatype, dimensions)
        new_var.units = units
        new_var.description = description

    return ncf


def write_netcdf(ncf: Dataset, results: Dict) -> None:
    """ Write a simulation result dictionary into a netcdf file.

    The variables that can be written are defined in the src.model_variables file

    Args:
        ncf (netCDF4.Dataset): the netcdf file where the results are written
        results (Dict): the results dictionary. Use the
            [tree_properties_to_dict](index.html#src.tools.iotools.tree_properties_to_dict) function
            to create the dictionary.

    """
    ind: int = ncf['index'].shape[0]
    ncf['index'][ind] = ind
    for key in results.keys():
        if (key in index_dim_vars
                or key in gas_index_dim_vars):
            # variable has dimension ("index")
            ncf[key][ind] = results[key]
        elif (key in root_dim_vars
              or key in soil_dim_vars
              or key in axial_layer_dim_vars
              or key in gas_axial_dim_vars):
            ncf[key][ind, :] = results[key].reshape(len(results[key]),)
        elif key in gas_three_dim_vars:
            ncf[key][ind, :, :, :] = results[key]
        else:
            # the variable has dimension
            # ("index", "radial_layers", "axial_layers")
            ncf[key][ind, :, :] = results[key]


def tree_properties_to_dict(tree: Tree) -> Dict:
    """ Transfers tree properties into a dictionary.

    Args:
        tree (Tree): Instance of the tree class.

    Returns:
        (Dict): Dictionary of the tree properties.

    """
    properties = {}
    properties['height'] = tree.height
    properties['num_elements'] = tree.num_elements
    properties['transpiration_rate'] = tree.transpiration_rate
    properties['photosynthesis_rate'] = tree.photosynthesis_rate
    sugar_conc = tree.sugar_concentration_as_numpy_array().reshape(tree.num_elements, 1)
    zeros = np.zeros((tree.num_elements, 1))
    properties['sugar_concentration'] = np.concatenate((zeros, sugar_conc), axis=1)
    properties['sugar_loading_rate'] = tree.sugar_loading_rate
    properties['sugar_unloading_rate'] = tree.sugar_unloading_rate
    properties['axial_permeability'] = tree.axial_permeability
    properties['radial_hydraulic_conductivity'] = np.repeat(
        tree.radial_hydraulic_conductivity.reshape(tree.num_elements, 1), 2, axis=1)
    properties['viscosity'] = tree.viscosity
    properties['elastic_modulus'] = tree.elastic_modulus
    properties['pressure'] = tree.pressure
    properties['radius'] = tree.element_radius[:, 1:]
    properties['area'] = np.concatenate([tree.element_area([], 0), tree.element_area([], 1)], axis=1)
    properties['volume'] = np.concatenate([tree.element_volume([], 0), tree.element_volume([], 1)], axis=1)

    return properties


def gas_properties_to_dict(gas: Gas) -> Dict:
    """ Transfers gas properties into a dictionary.

    Args:
        gas (Gas): Instance of the Gas class.

    Returns:
        (Dict): Dictionary of the gas properties.

    """
    properties = {}
    properties['gas_concentration'] = gas.concentration
    properties['gas_space_division'] = gas.space_division
    properties['gas_element_radii'] = gas.element_radius
    properties['gas_element_height'] = gas.element_height
    properties['gas_diffusion_coef'] = gas.diffusion_coefficients
    properties['gas_eq_rate'] = gas.equilibration_rate
    properties['gas_velocity'] = gas.velocity
    properties['gas_henry_coef'] = gas.kh
    properties['gas_temperature'] = gas.temperature
    properties['gas_concentration_ambient'] = gas.ambient_concentration
    properties['gas_moles_out'] = gas.n_out

    return properties
