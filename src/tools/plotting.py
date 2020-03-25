from typing import List
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


def single_property_plot(filename: str, property: str, ind: List) -> None:
    if(not(isinstance(ind[1], int) and isinstance(ind[2], int))):
        raise ValueError('ind[1] and ind[2] cannot be longer than 1')
    # load the property
    data = xr.open_dataset(filename)
    propertydata = data.__getattr__(property).isel(index=ind[0], axial_layers=ind[1], radial_layers=ind[2])
    indexdata = np.asarray(data.simulation_time.isel(index=ind[0]))

    labels = ['{:03.1f}'.format(index) for index in indexdata]

    propertydata.plot.line()
    plt.xticks(ticks=data.index.isel(index=ind[0]), labels=labels)
    plt.xlabel('Time [s]')
