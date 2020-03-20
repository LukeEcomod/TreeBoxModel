from .test_model import test_model
from .test_tree import test_tree
from typing import List
from ..tools.iotools import initialize_netcdf, write_netcdf
import os
import pytest


@pytest.fixture(scope="function")
def iomodel(test_model):

    test_model.outputfile = "iotools_test.nc"
    return test_model


def test_ncf_inif(iomodel):
    variables = {'radius': ['radius of an element', 'm', ("index", "radial_layers", "axial_layers"), 'f4'],
                 'pressures': ['pressure of an element', 'Pa', ("index", "radial_layers", "axial_layers"), 'f4'],
                 'elastic modulus':
                 ['elastic modulus of an element', 'Pa', ("index", "radial_layers", "axial_layers"), 'f4'],
                 'num_elements':
                 ['number of tree elements', 'pcs.', ("index", "radial_layers", "axial_layers"), 'i4'],
                 'simTime': ['Time in simulation', 's', ("index", "radial_layers", "axial_layers"), 'f4']}
    ncf = initialize_netcdf(iomodel, variables)
    assert ncf.variables['simTime'].units == 's'
    ncf.close()
    # os.remove(iomodel.outputfile)
