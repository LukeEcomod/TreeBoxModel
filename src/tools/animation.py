import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


data = xr.open_dataset('conceptual_gas_data_root_stem.nc')
na = len(data.axial_layers)
nr = len(data.radial_layers)
plt.xkcd()
fig = plt.figure()
ax = plt.axes()
start_data = np.zeros((na, nr))
im = ax.imshow(start_data, vmin=-6, vmax=-2, cmap='jet')
cbar = fig.colorbar(im, ticks=[-6, -2])
cbar.ax.set_yticklabels(['Ambient', 'Higher'])
cbar.ax.set_ylabel('CH$_{4}$ Concentration')
plt.xlabel('Stem radius (cm)')
plt.ylabel('Tree height (m)')
plt.xticks([0, 4, 9, 14, 19, 24], labels=[0, 2, 4, 6, 8, 10])
plt.yticks([0, 19, 39, 59, 79, 99], labels=[10, 8, 6, 4, 2, 0])
len_data = len(data.index)
frames = 200
increment = int(len_data/frames)
indeces_plot = np.arange(0, len_data, increment)


def init():
    '''Initialize image data.'''
    im.set_data(start_data)
    return [im]


def animate(i):
    '''Function for updating the image plot'''
    slide_data = np.log10(np.array(data.gas_concentration[indeces_plot[i], :, :, 1]))
    im.set_array(slide_data)
    plt.title(f'Time= {np.round(indeces_plot[i]*30*60/86400, 1)} days')
    return [im]


anim = FuncAnimation(fig, animate, init_func=init,
                     frames=frames, interval=10, blit=True)


anim.save('conceptual_test_root_stem_xkcd.gif', writer='imagemagick')
