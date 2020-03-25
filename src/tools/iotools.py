from typing import Dict
from netCDF4 import Dataset, Variable
# from ..model import Model # TODO: check why this import fails
from ..constants import MAX_ELEMENT_COLUMNS
from ..model_variables import index_dim_vars
from ..tree import Tree
import os.path
import numpy as np


def initialize_netcdf(model, variables: Dict) -> Dataset:

    # create the file if it does not exist
    if os.path.isfile(model.outputfile):
        ncf = Dataset(model.outputfile, "a")
    else:
        ncf = Dataset(model.outputfile, "w", format="NETCDF4")

    # create dimensions
    ncf.createDimension("index", 0)
    ncf.createDimension("radial_layers", MAX_ELEMENT_COLUMNS)
    ncf.createDimension("axial_layers", model.tree.num_elements)

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
    ind: int = ncf['index'].shape[0]
    ncf['index'][ind] = ind
    for key in results.keys():
        if key in index_dim_vars:
            # variable has dimension ("index")
            ncf[key][ind] = results[key]
        else:
            # the variable has dimension
            # ("index", "radial_layers", "axial_layers")
            ncf[key][ind, :, :] = results[key]


def tree_properties_to_dict(tree: Tree) -> Dict:
    """ cast tree properties to dictionary for saving"""
    properties = {}
    properties['height'] = tree.height
    properties['num_elements'] = tree.num_elements
    properties['ground_water_potential'] = tree.ground_water_potential
    properties['transpiration_rate'] = tree.element_property_as_numpy_array('transpiration_rate')
    properties['photosynthesis_rate'] = tree.element_property_as_numpy_array('photosynthesis_rate')
    sugar_conc = tree.sugar_concentration_as_numpy_array().reshape(40, 1)
    zeros = np.zeros((40, 1))
    properties['sugar_concentration'] = np.concatenate((zeros, sugar_conc), axis=1)
    properties['sugar_loading_rate'] = tree.element_property_as_numpy_array('sugar_loading_rate')
    properties['sugar_unloading_rate'] = tree.element_property_as_numpy_array('sugar_unloading_rate')
    properties['axial_permeability'] = tree.element_property_as_numpy_array('permeability')
    properties['radial_hydraulic_conductivity'] = tree.element_property_as_numpy_array('hydraulic_conductivity')
    properties['viscosity'] = tree.element_property_as_numpy_array('viscosity')
    properties['elastic_modulus'] = tree.element_property_as_numpy_array('elastic_modulus')
    properties['pressure'] = tree.element_property_as_numpy_array('pressure')
    properties['radius'] = tree.element_property_as_numpy_array('radius')

    return properties
