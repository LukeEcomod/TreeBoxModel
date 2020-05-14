from typing import Dict, List
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

    plt.plot(sim_time/3600, p_top, 'k-', label='Top')
    plt.plot(sim_time/3600, p_bottom, 'r-', label='Bottom')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (MPa)')
    plt.title('Phloem pressure')
    plt.xticks(np.linspace(0, 24, 7))
    filename = filename.split('.')[0]+'_phloem_pressure.png'
    plt.savefig(fname=filename, format='png')
    plt.close()


def plot_variable_vs_time(filename: str, params: Dict = None):
    if(params is None):
        params = {}
    data: xr.Dataset = xr.open_dataset(filename)

    variable: xr.DataArray = data[params['variable_name']]
    time: xr.DataArray = data['simulation_time']
    if(any(s == 'cut' for s in params.keys())):
        # cut needs to be a dictionary with at least one key from
        # keys: index, axial_layers and radial_layers
        variable = variable[params['cut']]
        time = time[{key: params['cut'][key] for key in ['index']}]
        time = time - data['simulation_time'][params['cut']['index'][0]]
    if(not any(s == 'time_divide' for s in params.keys())):
        params['time_divide'] = 1
    if(not any(s == 'variable_divide' for s in params.keys())):
        params['variable_divide'] = 1
    fig, ax = plt.subplots(figsize=(15, 8))
    # assume that the first axis in 'variable' is time
    # and there will be alltogether variable.shape[1]*variable.shape[2]
    # lines
    variable_shape = variable.shape
    for row_number in np.arange(variable_shape[1]):
        for col_number in np.arange(variable_shape[2]):
            ax.plot(time/params['time_divide'], variable[:, row_number, col_number]/params['variable_divide'],
                    label=params['labels'][row_number][col_number], color=params['line_colors'][row_number][col_number],
                    linewidth=params['line_widths'][row_number][col_number])

    ax.legend()
    plt.title(params['title'])
    plt.xlabel(params['xlabel'])
    plt.ylabel(params['ylabel'])
    if(any(s == 'xticks' for s in params.keys())):
        plt.xticks(params['xticks'])

    if(not any(s == 'folder' for s in params.keys())):
        filename = filename.split('.')[0]+'_' + params['filename_ending'] + '.png'
    else:
        filename = params['folder'] + params['filename_ending'] + '.png'

    plt.savefig(fname=filename, format='png')
    plt.close()


def plot_simulation_results(filename: str, foldername: str):
    params = {'variable_name': '',
              'time_divide': 86400,
              'labels': [['top'], ['bottom'], ['middle']],
              'line_colors': [['k'], ['r'], ['b']],
              'line_widths': [[3], [3], [3]],
              'cut': {'index': range(1500),
                      'axial_layers': [0, 39],
                      'radial_layers': [1]},
              'folder': foldername,
              'xticks': np.linspace(0, 3, 7),
              'xlabel': 'Time (d)'}

    plt.rcParams.update({'font.size': 22})

    variable_names = ['pressure', 'pressure', 'sugar_concentration', 'dqax', 'dqax', 'dqrad']
    filename_endings = ['xylem_pressure', 'phloem_pressure', 'sugar_concentration', 'xylem_dqax', 'phloem_dqax',
                        'phloem_dqrad']
    titles = ['Xylem pressure', 'Phloem pressure', 'Sugar concentration', 'Axial sap flux in the xylem',
              'Axial sap flux in the phloem', 'Radial sap flux in the phloem']
    ylabels = ['Pressure (Mpa)', 'Pressure (Mpa)', 'Sugar concentration (mol/L)', 'Flux (g/s)', 'Flux (g/s)',
               'Flux (g/s)']
    variable_divides = [1e6, 1e6, 1e3, 1e-3, 1e-3, 1e-3]
    end_index = 1153
    start_index = 288
    cuts = [{'index': range(start_index, end_index, 1),
             'axial_layers': [0, 39, 19],
             'radial_layers': [0]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [0, 39, 19],
             'radial_layers': [1]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [0, 39, 19],
             'radial_layers': [1]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [0, 39, 19],
             'radial_layers': [0]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [0, 39, 19],
             'radial_layers': [1]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [0, 39, 19],
             'radial_layers': [1]}]

    for (ind, var) in enumerate(variable_names):
        params['variable_name'] = var
        params['title'] = titles[ind]
        params['ylabel'] = ylabels[ind]
        params['variable_divide'] = variable_divides[ind]
        params['cut'] = cuts[ind]
        params['filename_ending'] = filename_endings[ind]
        plot_variable_vs_time(filename, params)
