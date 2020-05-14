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
    properties['transpiration_rate'] = tree.transpiration_rate
    properties['photosynthesis_rate'] = tree.photosynthesis_rate
    sugar_conc = tree.sugar_concentration_as_numpy_array().reshape(40, 1)
    zeros = np.zeros((40, 1))
    properties['sugar_concentration'] = np.concatenate((zeros, sugar_conc), axis=1)
    properties['sugar_loading_rate'] = tree.sugar_loading_rate
    properties['sugar_unloading_rate'] = tree.sugar_unloading_rate
    properties['axial_permeability'] = tree.axial_permeability
    properties['radial_hydraulic_conductivity'] = np.repeat(tree.radial_hydraulic_conductivity.reshape(40, 1),
                                                            2, axis=1)
    properties['viscosity'] = tree.viscosity
    properties['elastic_modulus'] = tree.elastic_modulus
    properties['pressure'] = tree.pressure
    properties['radius'] = tree.element_radius
    properties['area'] = np.concatenate([tree.element_area([], 0), tree.element_area([], 1)], axis=1)
    properties['volume'] = np.concatenate([tree.element_volume([], 0), tree.element_volume([], 1)], axis=1)

    return properties
