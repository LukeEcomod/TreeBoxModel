from typing import Dict
from netCDF4 import Dataset
from ..model import Model
from ..constants import MAX_ELEMENT_COLUMNS
import os.path


def initialize_netcdf(model: Model, variables: Dict) -> Dataset:

    # create the file if it does not exist
    if os.path.isfile(model.outputfile):
        ncf = Dataset(model.outputfile, "a")
    else:
        ncf = Dataset(model.outputfile, "w", format="NETCDF4")

    # create dimensions
    ncf.createDimension("index", 0)
    ncf.createDimension("radial_layers", MAX_ELEMENT_COLUMNS)
    ncf.createDimension("axial_layers", model.tree.num_elements)

    # create variables
    for var_name in variables.keys():
        description: str = variables[var_name][0]
        units: str = variables[var_name][1]
        dimensions: str = variables[var_name][2]
        datatype: str = variables[var_name][3]
        new_var = ncf.createVariable(var_name, datatype, dimensions)
        new_var.units = units
        new_var.description = description

    return ncf


def write_netcdf(ncf: Dataset, results: Dict):
    ind: int = ncf['index'].shape[0] - 1
    index_dim_vars: List = ['num_elements', 'simTime', 'index']
    for key in results.keys():
        if key in index_dim_vars:
            # variable has dimension ("index")
            ncf[key][ind] = results[key]
        else:
            # the variable has dimension
            # ("index", "radial_layers", "axial_layers")
            ncf[key][ind, :, :] = results[key]
