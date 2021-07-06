from typing import Dict, List
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

font = {'size': 12}
plt.rc('font', **font)


def plot_xylem_pressure_top_bottom(filename: str) -> None:
    """ Plot the xylem pressure at the top and bottom of the tree.

        The figure is saved in the current working directory with the same
        name as the filename with "_xylem_pressure.png" appended.

        Args:
            filename (str): name of the NetCDF file in the current directory that
                includes the simulation results.
    """
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


def plot_phloem_pressure_top_bottom(filename: str) -> None:
    """ Plot the phloem pressure at the top and bottom of the tree.

    The figure is saved in the current working directory with the same
    name as the filename with "_phloem_pressure.png" appended.

    Args:
        filename (str): name of the NetCDF file in the current directory that
            includes the simulation results.
    """
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


def plot_variable_vs_time(filename: str, params: Dict = None) -> None:
    """Plot any variable as a function of time from the NetCDF file that contains the simulation results.

    Args:
        filename (str): name of the NetCDF file in the current directory that
            includes the simulation results.
        params (Dict, optional): Parameters and instructions for making the plot (see Examples)

    Examples:
        An example of the params dictionary is
        ```python
        params = {'variable_name': 'pressure',
        'time_divide': 86400,
        'variable_divide': 1e6
        'labels': [['top'], ['bottom'], ['middle']],
        'line_colors': [['k'], ['r'], ['b']],
        'line_widths': [[3], [3], [3]],
        'cut': {'index': range(1500),
                'axial_layers': [0, 39],
                'radial_layers': [1]},
        'xticks': np.linspace(0, 10, 11),
        'xlabel': 'Time (d)'
        'ylabel': 'Pressure (MPa)',
        'title': 'Pressure in the xylem'
        'folder', 'figure/',
        'filename_ending': 'xylem_pressure.png'}
        ```
        * Variable name is the name of the variable that is to be plotted. The name must equal to a variable
        in the NetCDF file
        * time_divide is a float which is used to divide the time vector (x-axis) in the plot. e.g., 86400 converts
        seconds to days
        * variable_divide is a float which is used to divide the variable vector (y-axis) in the plot
        * labels are the labels of each line that is drawn
        * line_colors are the colors of each line that is drawn
        * line_widths are the line widths of each line that is drawn
        * cut contains the data range in the NetCDF file that are plotted. The function can handle only
        3-dimensional data. If the variable that is plotted has shape (1500,2,1) like in this example,
        the function expects that there will be 2\\*1=3 lines in the plot. Consequently, the labels, line_colors,
        and line_widhths values in the params dictionary need to have 3 elements like in this example.
        * xticks are the matplotlib.pyplot.xticks function argument
        * xlabel is the matplotlib.pyplot.xlabel function argument
        * ylabel is the matplotlib.pyplot.ylabel function argument
        * title is the matplotlib.pyplot.title function argument
        * folder refers to the folder starting from the current working directory where the figure is saved
        * filename_ending is a string that is appended to the filename when the figure is saved
        * If no folder is specified the argument filename and filename_ending is used to create the name of the figure
        that is saved.
    """
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


def plot_simulation_results(filename: str, foldername: str) -> None:
    """Plot sugar concentration, fluxes and pressures from one simulation.

    Currently you need edit tha start and end indeces in the code to capture
    the correct time window in the simulation results. The indices refer to
    the index dimension in the NetCDF files.

    Args:
        filename (str):  name of the NetCDF file in the current directory that
            includes the simulation results.
        foldername (str): name of the folder where the figures are saved
    """
    data = xr.open_dataset(filename)
    params = {'variable_name': '',
              'time_divide': 86400,
              'labels': [['top'], ['bottom'], ['middle']],
              'line_colors': [['k'], ['r'], ['b']],
              'line_widths': [[3], [3], [3]],
              'cut': {'index': range(1500),
                      'axial_layers': [0, 39],
                      'radial_layers': [1]},
              'folder': foldername,
              'xticks': np.linspace(0, 10, 11),
              'xlabel': 'Time (d)'}

    plt.rcParams.update({'font.size': 22})

    variable_names = ['pressure', 'pressure', 'sugar_concentration', 'dqax', 'dqax', 'dqax_up',
                      'dqax_up', 'dqax_down', 'dqax_down', 'dqrad', 'transpiration_rate',
                      'photosynthesis_rate']
    filename_endings = ['xylem_pressure', 'phloem_pressure', 'sugar_concentration', 'xylem_dqax', 'phloem_dqax',
                        'xylem_dqax_up', 'phloem_dqax_up', 'xylem_dqax_down', 'phloem_dqax_down',
                        'phloem_dqrad', 'transpiration_rate', 'photosynthesis_rate']
    titles = ['Xylem pressure', 'Phloem pressure', 'Sugar concentration', 'Axial sap flux in the xylem',
              'Axial sap flux in the phloem', 'Axial upward sap flux in the xylem',
              'Axial upward sap flux in the phloem',
              'Axial downward sap flux in the xylem', 'Axial downward sap flux in the phloem',
              'Radial sap flux in the phloem', 'Transpiration rate',
              'Photosynthesis rate']
    ylabels = ['Pressure (Mpa)', 'Pressure (Mpa)', 'Sugar concentration (mol/L)', 'Flux (g/s)', 'Flux (g/s)',
               'Flux (g/s)', 'Flux (g/s)', 'Flux (g/s)', 'Flux (g/s)', 'Flux (g/s)',
               'Transpiration rate (g/s)', r'Photosynthesis rate ($\mu$mol/s)']
    variable_divides = [1e6, 1e6, 1e3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-6]
    end_index = 433
    start_index = 10
    bottom = int(data.num_elements[0]-1)
    middle = int((data.num_elements[0]-1)/2)
    top = 1
    cuts = [{'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [0]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [1]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [1]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [0]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [1]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [0]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [1]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [0]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [1]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [1]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [0]},
            {'index': range(start_index, end_index, 1),
             'axial_layers': [top, bottom, middle],
             'radial_layers': [1]}]

    for (ind, var) in enumerate(variable_names):
        params['variable_name'] = var
        params['title'] = titles[ind]
        params['ylabel'] = ylabels[ind]
        params['variable_divide'] = variable_divides[ind]
        params['cut'] = cuts[ind]
        params['filename_ending'] = filename_endings[ind]
        plot_variable_vs_time(filename, params)


def plot_ax_up_change(filenames):
    colors = ['k', 'r', 'b']
    lines = ['-', '--', '-.']
    legends = ['base case $k_{base} = 1.5 \\cdot 10^{-12}$ at bottom',
               filenames[1][15:-3].replace('_', ' '),
               filenames[2][15:-3].replace('_', ' ')]
    for (ind, f) in enumerate(filenames):
        data = xr.open_dataset(f)
        time = np.asarray(data.simulation_time[10:]/86400)
        middle_ind = int(data.num_elements[0]/2)
        dqax_up = np.abs(np.asarray(data.dqax_up[10:, middle_ind, 0])*1e3)
        plt.plot(time, dqax_up, color=colors[ind], linestyle=lines[ind], label=legends[ind])

    data = xr.open_dataset(filenames[0])
    time = np.asarray(data.simulation_time[10:]/86400)
    e = np.sum(np.asarray(data.transpiration_rate[10:, :, 0]), axis=1)*1e3
    plt.plot(time, e, marker='o', mec='k', mfc='k', markersize=5, linestyle='', label='Total transpiration')
    plt.title('$|Q_{ax,up}|$ in the xylem')
    plt.ylabel('Flux (g/s)')
    plt.xlabel('Time (d)')
    plt.legend()

# def plot_two_simulations_vs_time()

# def plot_two_simulation_comparison()
