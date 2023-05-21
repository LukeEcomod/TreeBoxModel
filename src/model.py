import datetime
from typing import Dict
import numpy as np
from netCDF4 import Dataset
from scipy.integrate import solve_ivp
from .odefun_tree import odefun_tree
from .model_variables import all_variables
from .tree import Tree
from .soil import Soil
from .constants import (GRAVITATIONAL_ACCELERATION,
                        RHO_WATER,
                        MOLAR_GAS_CONSTANT,
                        MAX_ELEMENT_COLUMNS)
from .tools.iotools import (initialize_netcdf,
                            write_netcdf)
from .tools.tree_to_gas import convert_tree_flux_to_velocity


class Model:
    """ Calculates the next time step for given tree and saves the tree stage.

    Provides functionality for solving the ordinary differential equations (ODE) describing
    the behaviour of the [modelled system](modelled_system.html).

    Args:
        tree (Tree): instance of the tree class for which the ODEs are solved
        outputfile (str): name of the file where the NETCDF4 output is written

    Attributes:
        tree (Tree): instance of the tree class for which the ODEs are solved
        ncf (netCDF4.Dataset): the output file
    """

    def __init__(self, tree: Tree, soil: Soil, outputfile: str = ''):
        self.tree: Tree = tree
        self.soil = soil
        self.outputfile = outputfile
        if len(outputfile) != 0:
            dims = {'axial_layers': tree.num_elements, 'radial_layers': 2,
                    'root_elements': tree.roots.num_elements, 'soil_elements': soil.num_elements}
            self.ncf: Dataset = initialize_netcdf(outputfile, dims, all_variables)

    def axial_fluxes(self) -> np.ndarray:
        """Calculates axial sap mass flux for every element.

        The axial flux in the xylem and phloem are calculated independently
        from the sum of bottom and top fluxes.

        .. math::
            Q_{ax,i} = Q_{ax,bottom,i} + Q_{ax,top,i} - E + Q_{root, i}

        .. math::
            Q_{ax,bottom,i} = \\frac{k_i \\: A_{ax,i} \\: \\rho_w}{\\eta_i \\: l_i}(P_{i+1} - P_{i} - P_h)

            Q_{ax,top,i} = \\frac{k_i \\: A_{ax,i+1} \\: \\rho_w}{\\eta_i \\: l_i}(P_{i-1} - P_{i} + P_h)

        where

        * :math:`E_i`: transpiration rate of the ith element (:math:`\\frac{kg}{s}`)
        * :math:`k_i`: axial permeability of the ith element (:math:`m^2`)
        * :math:`A_{ax,i}`: base surface area of xylem or phloem (:math:`m^2`)
        * :math:`\\rho_w`: liquid phase density of water (:math:`\\frac{kg}{m^3}`)
        * :math:`\\eta`: viscosity of the sap in the ith element (:math:`Pa \\: s`)
        * :math:`l_i`: length (height) of the ith element (:math:`m`)
        * :math:`P_{i}`: Pressure in the ith element (:math:`Pa`)
        * :math:`P_h`: Hydrostatic pressure (:math:`Pa`) :math:`P_h = \\rho_w a_{gravitation} l_i`

        Returns:
            numpy.ndarray (dtype=float, ndim=2)[self.tree.num_elements, 2]: The axial fluxes in units kg/s

        """

        pressures: np.ndarray = self.tree.pressure
        k: np.ndarray = self.tree.axial_permeability
        eta: np.ndarray = self.tree.viscosity
        length: np.ndarray = self.tree.element_height
        l: np.ndarray = np.concatenate(([0], length.reshape(len(length),)))
        cumulative_sum: np.ndarray = np.cumsum(l).reshape(len(l), 1)
        l_midpoints = length/2 + cumulative_sum[:-1]
        dl_midpoints = np.diff(l_midpoints, axis=0)
        E: np.ndarray = self.tree.transpiration_rate
        C: np.ndarray = self.tree.sugar_concentration_as_numpy_array()
        C = np.concatenate([np.zeros((self.tree.num_elements, 1)), C], axis=1)
        RWU: np.ndarray = self.root_fluxes()
        A = np.concatenate([self.tree.element_area([], 0),
                            self.tree.element_area([], 1)], axis=1)
        # calculate transport coefficients
        # TODO: add calculation for phloem sap density
        transport_ax: np.ndarray = k[1:, :]/eta[1:, :]/np.repeat(dl_midpoints, 2, axis=1)*RHO_WATER * A[1:, :]

        # calculate downward and upward fluxes separately
        Q_ax_down: np.ndarray = np.zeros((self.tree.num_elements, pressures.shape[1]))
        Q_ax_down[0:-1, :] = (np.diff(pressures, axis=0)
                              - RHO_WATER
                              * GRAVITATIONAL_ACCELERATION
                              * np.repeat(dl_midpoints, 2, axis=1)
                              ) * transport_ax
        Q_ax_down[-1, 0] = 0  # flux between xylem and soil is handled in the RWU
        Q_ax_down[-1, 1] = 0  # no flux from phloem to soil

        Q_ax_up: np.ndarray = np.zeros((self.tree.num_elements, pressures.shape[1]))

        Q_ax_up[1:, :] = -1.0*(np.diff(pressures, axis=0) - RHO_WATER*GRAVITATIONAL_ACCELERATION
                               * np.repeat(dl_midpoints, 2, axis=1))*transport_ax
        Q_ax_up[0, :] = 0  # no upward flux in the highest element
        Q_ax: np.ndarray = (Q_ax_up + Q_ax_down)

        Q_ax[:, 0] = Q_ax[:, 0]-E.reshape(self.tree.num_elements,)+RWU.reshape(self.tree.num_elements,)

        return Q_ax, Q_ax_up, Q_ax_down

    def radial_fluxes(self) -> np.ndarray:
        """ Calculates radial sap mass flux for every element.

        The radial flux for the phloem of the ith axial is calculated similar to
        [Hölttä et. al. 2006](https://link.springer.com/article/10.1007/s00468-005-0014-6)

        .. math::
            Q_{radial,phloem} = L_r A_{rad,i}\\rho_{w}
            [P_{i,xylem} - P_{i,phloem} - \\sigma(C_{i,xylem} - C_{i,phloem})RT)]

        where

        * :math:`L_r`: radial hydraulic conductivity (:math:`\\frac{m}{Pa \\: s}`)
        * :math:`A_{rad,i}`: lateral surface area of the xylem (:math:`m^2`)
        * :math:`\\rho_w`: liquid phase density of water (:math:`\\frac{kg}{m^3}`)
        * :math:`P_{i}`: Pressure in the ith element (:math:`Pa`)
        * :math:`\\sigma`: Reflection coefficient (Van't hoff factor) (unitless)
        * :math:`C_{i}`: Sucrose concentration in the ith element (:math:`\\frac{mol}{m^3}`)
        * :math:`R`: Universal gas constant (:math:`\\frac{J}{K \\: mol}`)
        * :math:`T`: Ambient temperature (:math:`K`)

        The radial flux for the xylem is equal to the additive inverse of the phloem flux

        .. math::
            Q_{radial,xylem} = -Q_{radial,phloem}


        Returns:
            numpy.ndarray (dtype=float, ndim=2)[self.tree.num_elements, 2]: The radial fluxes in units kg/s
        """

        pressures: np.ndarray = self.tree.pressure
        Lr: np.ndarray = self.tree.radial_hydraulic_conductivity
        C: np.ndarray = self.tree.sugar_concentration_as_numpy_array()
        Q_rad_phloem: np.ndarray = Lr.reshape(
            (self.tree.num_elements, 1)) * self.tree.cross_sectional_area() * RHO_WATER * (
            np.diff(np.flip(pressures, axis=1), axis=1) + C.reshape((self.tree.num_elements, 1))
            * MOLAR_GAS_CONSTANT*self.tree.temperature)
        Q_rad_xylem: np.ndarray = -(Q_rad_phloem.copy())
        return np.concatenate((Q_rad_xylem, Q_rad_phloem), axis=1)

    def root_fluxes(self):
        """ Calculates root water uptake for every element. If the part of the tree is not part
        of the root system, the flux is set to zero.

        Conductance between soil and root xylem is calculated similar to
        [Volpe et. al. 2013](https://doi.org/10.1016/j.advwatres.2013.07.008)

        .. math::
            Q_{root,i} = \\frac{G_i}{g} (P_{soil,i} - P_{root,xylem,i})

        where

        * :math:`G_i`: Total conductance from soil to root xylem (:math:`\\frac{m^2}{s}`). See
            [Roots class](index.html#src.roots.Roots.conductivity) for details
        * :math:`g`: gravitational acceleration (:math:`\\frac{m}{s^2}`)
        * :math:`P`: Pressure in either the soil element or root xylem element

        Returns:
            numpy.ndarray (dtype=float, ndim=2)[self.tree.num_elements, 1]: The root water uptake in units kg/s

        References:
            Volpe, V. et. al., "Root controls on water redistribution and carbon uptake in the soil-plant
            system under current and future climate", Advances in Water Resources, 60, 110-120, 2013.
        """

        ind = self.tree.root_elements
        soil_ind = self.tree.roots.soil_elements(self.soil)
        gi: np.ndarray = self.tree.roots.conductivity(self.soil)
        length: np.ndarray = self.tree.element_height
        P_root = self.tree.pressure[ind, 0].reshape(len(ind), 1)
        P_soil = self.soil.pressure[soil_ind].reshape(len(ind), 1)
        Q_root = np.zeros((self.tree.num_elements, 1))
        Q_root[ind, 0] = (gi/GRAVITATIONAL_ACCELERATION*(P_soil-P_root)).reshape(len(ind),)

        return Q_root

    def run(self, time_start: float = 1e-3, time_end: float = 120.0,
            dt: float = 0.01, output_interval: float = 60) -> None:
        """ Propagates the tree in time using explicit Euler method (very slow).

        NB! This function needs to be updated. Use run_scipy instead!

        Args:
            time_start (float): Time in seconds where to start the simulation.
            time_ned (float): Time in seconds where to end the simulation.
            dt (float): time step in seconds
            output_interval: Time interval in seconds when to save the tree stage

        """
        # TODO: This function needs to be updated so that dr/dt is calculated like in odefun.py
        dmdt_ax: np.ndarray = np.zeros((self.tree.num_elements, 2))
        dmdt_rad: np.ndarray = np.zeros((self.tree.num_elements, 2))
        for (ind, time) in enumerate(np.linspace(time_start, time_end, int((time_end-time_start)/dt))):
            # get the change in every elements mass

            dmdt_ax, _, _ = self.axial_fluxes()
            dmdt_rad = self.radial_fluxes()
            if np.abs(time % output_interval) < 1e-2:
                print(datetime.datetime.now(), "\t", time)
                results = self.properties_to_dict()
                results['simulation_time'] = time
                results['model_index'] = ind
                write_netcdf(self.ncf, results)

            self.tree.pressure += dt*self.tree.elastic_modulus/(np.transpose(
                np.array([self.tree.element_volume_water([], 0),
                          self.tree.element_volume_water([], 1)])) * RHO_WATER)\
                * (dmdt_ax + dmdt_rad)

            for i in range(self.tree.num_elements):
                self.tree.solutes[i, 1].concentration += dt / self.tree.element_volume_water([i], 1)\
                    * (dmdt_ax[i, 1]
                       * self.tree.solutes[i, 1].concentration/RHO_WATER
                       + self.tree.sugar_loading_rate[i]
                       - self.tree.sugar_unloading_rate[i])
                if i > self.tree.num_elements - 6:
                    self.tree.sugar_unloading_rate[i] = (self.tree.solutes[i, 1].concentration
                                                         - self.tree.sugar_target_concentration)*1e-5
            self.tree.element_radius += (dmdt_ax + dmdt_rad)\
                / (np.pi*RHO_WATER * np.repeat(self.tree.element_height, 2, axis=1) * self.tree.element_radius)

            # update sap viscosity
            self.tree.update_sap_viscosity()
            return results

    def run_scipy(self, time_start: float = 1e-3, time_end: float = 120.0, ind: int = 0) -> None:
        """ Propagates the tree in time using the solve_ivp function in the SciPy package.

        The stage of the tree is saved only at the start of the simulation if time_start < 1e-3 and at time_end.
        If the tree stage is desired on multiple time points the function needs to be called recurrently by splitting
        the time interval into multiple sub intervals.

        Args:
            time_start (float): Time in seconds where to start the simulation.
            time_ned (float): Time in seconds where to end the simulation.
            ind (int): index which refers to the index in self.outputfile.
                The last stage of the tree is saved to self.outputfile[ind].

        """
        # If time < 0 save the first stage of the tree
        # print(datetime.datetime.now(), "\t", time_end)
        if(time_start < 1e-3 and len(self.outputfile) != 0):
            results = self.properties_to_dict()
            results['simulation_time'] = time_start
            results['model_index'] = ind
            write_netcdf(self.ncf, results)

        # Initial values from self.tree
        initial_values = [self.tree.pressure,
                          self.tree.sugar_concentration_as_numpy_array().reshape(self.tree.num_elements, 1),
                          self.tree.element_radius]
        # broadcast initial values into 1D array
        yinit = np.concatenate([initial_values[0].reshape(self.tree.num_elements*2, order='F'),
                                initial_values[1].reshape(self.tree.num_elements, order='F'),
                                initial_values[2].reshape(self.tree.num_elements*3, order='F')])

        sol = solve_ivp(lambda t, y: odefun_tree(t, y, self), (time_start, time_end), yinit, method='BDF',
                        rtol=1e-6, atol=1e-3)
        last_tree_stage = sol.y[:, -1]
        pressures = last_tree_stage[0:self.tree.num_elements*2].reshape(
            self.tree.num_elements, MAX_ELEMENT_COLUMNS, order='F')
        sugar_concentration = last_tree_stage[self.tree.num_elements*2:self.tree.num_elements*3].reshape(
            self.tree.num_elements, 1, order='F')
        element_radius = last_tree_stage[self.tree.num_elements*3:].reshape(self.tree.num_elements,
                                                                            MAX_ELEMENT_COLUMNS+1, order='F')
        self.tree.pressure = pressures
        self.tree.update_sugar_concentration(sugar_concentration)
        self.tree.element_radius = element_radius
        self.tree.update_sap_viscosity()
        # save the tree status
        if len(self.outputfile) != 0:
            print(datetime.datetime.now(), "\t", time_end/86400)
            results = self.properties_to_dict()
            results['simulation_time'] = time_end
            results['model_index'] = ind

            write_netcdf(self.ncf, results)

    def properties_to_dict(self) -> Dict:
        """ Transfers tree properties into a dictionary.

            Args:
                tree (Tree): Instance of the tree class.

            Returns:
                (Dict): Dictionary of the tree properties.

            """
        properties = {}
        properties['height'] = self.tree.height
        properties['element_height'] = self.tree.element_height
        properties['num_elements'] = self.tree.num_elements
        properties['transpiration_rate'] = self.tree.transpiration_rate
        properties['photosynthesis_rate'] = self.tree.photosynthesis_rate
        sugar_conc = self.tree.sugar_concentration_as_numpy_array().reshape(self.tree.num_elements, 1)
        zeros = np.zeros((self.tree.num_elements, 1))
        properties['sugar_concentration'] = np.concatenate((zeros, sugar_conc), axis=1)
        properties['sugar_loading_rate'] = self.tree.sugar_loading_rate
        properties['sugar_unloading_rate'] = self.tree.sugar_unloading_rate
        properties['axial_permeability'] = self.tree.axial_permeability
        properties['radial_hydraulic_conductivity'] = np.repeat(
            self.tree.radial_hydraulic_conductivity.reshape(self.tree.num_elements, 1), 2, axis=1)
        properties['viscosity'] = self.tree.viscosity
        properties['elastic_modulus'] = self.tree.elastic_modulus
        properties['pressure'] = self.tree.pressure
        properties['radius'] = self.tree.element_radius[:, 1:]
        properties['area'] = np.concatenate([self.tree.element_area([], 0), self.tree.element_area([], 1)], axis=1)
        properties['volume_water'] = np.concatenate([self.tree.element_volume_water([], 0),
                                                     self.tree.element_volume_water([], 1)], axis=1)
        properties['sapflow'] = np.concatenate(convert_tree_flux_to_velocity(self), axis=1)
        properties['dqroot'] = self.root_fluxes()
        properties['dqrad'] = self.radial_fluxes()
        properties['dqax'], properties['dqax_up'], properties['dqax_down'] = self.axial_fluxes()
        properties['dqroot'] = self.root_fluxes()
        properties['soil_dz'] = self.soil.layer_thickness
        properties['soil_kh'] = self.soil.hydraulic_conductivity
        properties['soil_pressure'] = self.soil.pressure
        properties['soil_root_k'] = self.tree.roots.conductivity(self.soil)
        properties['rooting_depth'] = self.tree.roots.rooting_depth
        properties['root_area_density'] = self.tree.roots.area_density

        return properties
