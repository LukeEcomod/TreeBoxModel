import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


data = xr.open_dataset('testnr10_noflux.nc')
na = len(data.axial_layers)
nr = len(data.radial_layers)
fig = plt.figure()
ax = plt.axes()
start_data = np.zeros((na, nr))
im = ax.imshow(start_data, vmin=-6, vmax=2, cmap='jet')


def init():
    start_data = np.zeros((na, nr))
    im.set_data(start_data)
    return [im]


def animate(i):
    slide_data = np.log10(np.array(data.gas_concentration[i, :, :, 1]))
    im.set_array(slide_data)
    # print(slide_data[-1, :])
    return [im]


anim = FuncAnimation(fig, animate, init_func=init,
                     frames=200, interval=20, blit=True)


anim.save('conc_test.gif', writer='imagemagick')
