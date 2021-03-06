from typing import Dict
from netCDF4 import Dataset, Variable
# from ..model import Model # TODO: check why this import fails
from ..constants import MAX_ELEMENT_COLUMNS
from ..model_variables import index_dim_vars, soil_dim_vars, root_dim_vars, axial_layer_dim_vars
from ..tree import Tree
import os.path
import numpy as np


def initialize_netcdf(filename: str,
                      axial_elements: int,
                      soil_elements: int,
                      root_elements: int,
                      variables: Dict) -> Dataset:
    """ Initializes a netcdf file to be ready for saving simulation results.

        Args:
            filename (str): name of the NETCDF4 file that is created
            axial_elements (int): Number of axial elements in the Tree
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
    ncf.createDimension("radial_layers", MAX_ELEMENT_COLUMNS)
    ncf.createDimension("axial_layers", axial_elements)
    ncf.createDimension("soil_elements", soil_elements)
    ncf.createDimension("root_elements", root_elements)

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
        if key in index_dim_vars:
            # variable has dimension ("index")
            ncf[key][ind] = results[key]
        elif key in root_dim_vars or key in soil_dim_vars or key in axial_layer_dim_vars:
            ncf[key][ind, :] = results[key].reshape(len(results[key]),)
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
