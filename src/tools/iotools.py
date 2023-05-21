import os.path
from typing import Dict
import numpy as np
from netCDF4 import Dataset, Variable
from ..model_variables import (index_dim_vars, soil_dim_vars, root_dim_vars, axial_layer_dim_vars,
                               gas_three_dim_vars, gas_index_dim_vars, gas_axial_dim_vars)


def initialize_netcdf(filename: str,
                      dimensions: Dict,
                      variables: Dict) -> Dataset:
    """ Initializes a netcdf file to be ready for saving simulation properties.

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
            (NETCDF4.Dataset): A NETCDF4 file where the simulation properties can be saved.

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


def write_netcdf(ncf: Dataset, properties: Dict) -> None:
    """ Write a simulation result dictionary into a netcdf file.

    The variables that can be written are defined in the src.model_variables file

    Args:
        ncf (netCDF4.Dataset): the netcdf file where the properties are written
        properties (Dict): the properties dictionary. Use the
            [tree_properties_to_dict](index.html#src.tools.iotools.tree_properties_to_dict) function
            to create the dictionary.

    """
    ind: int = ncf['index'].shape[0]
    ncf['index'][ind] = ind
    for key in properties.keys():
        if (key in index_dim_vars
                or key in gas_index_dim_vars):
            # variable has dimension ("index")
            ncf[key][ind] = properties[key]
        elif (key in root_dim_vars
              or key in soil_dim_vars
              or key in axial_layer_dim_vars
              or key in gas_axial_dim_vars):
            ncf[key][ind, :] = properties[key].reshape(len(properties[key]),)
        elif key in gas_three_dim_vars:
            ncf[key][ind, :, :, :] = properties[key]
        else:
            # the variable has dimension
            # ("index", "radial_layers", "axial_layers")
            ncf[key][ind, :, :] = properties[key]


def compound_properties_to_dict(gas) -> Dict:
    """ Transfers compound properties into a dictionary.

    Args:
        gas (Gas): Instance of the Gas class.

    Returns:
        (Dict): Dictionary of the gas properties.

    """
    properties = {}
    properties['cmp_concentration'] = gas.concentration
    properties['cmp_space_division'] = np.stack((gas.space_division[:, :, 0],
                                                 np.sum(gas.space_division[:, :, 1:], axis=2)), axis=2)
    properties['cmp_element_radii'] = gas.element_length
    properties['cmp_element_height'] = gas.element_height
    properties['cmp_head_area'] = gas.head_area
    properties['cmp_element_volume'] = gas.element_volume
    properties['cmp_element_volume_air'] = gas.element_volume_air
    properties['cmp_element_volume_water'] = gas.element_volume_water
    properties['cmp_element_volume_cell'] = gas.element_volume_cell
    properties['cmp_gas_radial_diffusion_coefficient'] = gas.radial_diffusion_coefficient
    properties['cmp_gas_axial_diffusion_coefficient'] = gas.axial_diffusion_coefficient
    properties['cmp_eq_rate'] = gas.equilibration_rate
    properties['cmp_velocity'] = gas.velocity
    properties['cmp_henry_coef'] = gas.kh
    properties['cmp_temperature'] = gas.temperature
    properties['cmp_ambient_concentration'] = gas.ambient_concentration
    properties['cmp_moles_out'] = gas.n_out
    return properties
