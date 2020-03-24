from .test_model import test_model
from .test_tree import test_tree
from typing import List
from ..tools.iotools import initialize_netcdf, write_netcdf
import os
import pytest
import xarray as xr
import numpy as np


@pytest.fixture(scope="function")
def ncf(test_model):

    test_model.outputfile = "iotools_test.nc"
    variables = {'radius': ['radius of an element', 'm', ("index", "radial_layers", "axial_layers"), 'f4'],
                 'pressure': ['pressure of an element', 'Pa', ("index", "radial_layers", "axial_layers"), 'f4'],
                 'elastic modulus':
                 ['elastic modulus of an element', 'Pa', ("index", "radial_layers", "axial_layers"), 'f4'],
                 'num_elements':
                 ['number of tree elements', 'pcs.', ("index"), 'i4'],
                 'simTime': ['Time in simulation', 's', ("index"), 'f4']}

    ncf = initialize_netcdf(test_model, variables)
    yield ncf
    print('\n Cleaning up NetCDF4 file')
    ncf.close()
    os.remove(test_model.outputfile)
    print('\n Cleanup finished')


def test_ncf_inif(ncf):
    assert ncf.variables['simTime'].units == 's'


def test_ncf_writer(ncf):
    # add random data to pressures and simTime
    random_pressure_data = np.random.uniform(size=(2, 40, 2))
    random_simTime_data = np.random.uniform(size=(2, 1))
    results = {'simTime': random_simTime_data[0], 'pressure': random_pressure_data[0, :, :].reshape(40, 2)}
    write_netcdf(ncf, results)
    # write again
    results = {'simTime': random_simTime_data[1], 'pressure': random_pressure_data[1, :, :].reshape(40, 2)}
    write_netcdf(ncf, results)
    # read the netcdf file with xarray
    data = xr.open_dataset('iotools_test.nc')
    # check that the sums match
    assert float(np.sum(data['pressure'])) == pytest.approx(np.sum(random_pressure_data), rel=1e-6)
    assert float(np.sum(data['simTime'])) == pytest.approx(np.sum(random_simTime_data), rel=1e-6)
