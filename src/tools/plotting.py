from typing import List
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

font = {'size': 12}
plt.rc('font', **font)


def plot_xylem_pressure_top_bottom(filename: str):
    data = xr.open_dataset(filename)
    p_top = np.asarray(data.pressure[:, 0, 0])*1e-6
    p_bottom = np.asarray(data.pressure[:, 39, 0])*1e-6

    sim_time = np.asarray(data.simulation_time)

    plt.plot(sim_time/3600, p_top, 'k-', label='Top')
    plt.plot(sim_time/3600, p_bottom, 'r-', label='Bottom')
    plt.legend()
    plt.xlabel('Time (h)')
    plt.ylabel('Pressure (MPa)')
    plt.title('Xylem pressure')
    plt.xticks(np.linspace(0, 24, 7))
    filename = filename.split('.')[0]+'_xylem_pressure.png'
    plt.savefig(fname=filename, format='png')
    plt.close()


def plot_phloem_pressure_top_bottom(filename: str):
    data = xr.open_dataset(filename)
    p_top = np.asarray(data.pressure[:, 0, 1])*1e-6
    p_bottom = np.asarray(data.pressure[:, 39, 1])*1e-6

    sim_time = np.asarray(data.simulation_time)

    plt.plot(sim_time, p_top, 'k-', label='Top')
    plt.plot(sim_time, p_bottom, 'r-', label='Bottom')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (MPa)')
    plt.title('Phloem pressure')
    filename = filename.split('.')[0]+'_phloem_pressure.png'
    plt.savefig(fname=filename, format='png')
    plt.close()
