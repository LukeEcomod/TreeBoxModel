from typing import Dict
import numpy as np
import math
import datetime
from scipy.integrate import solve_ivp
from .tree import Tree
from .constants import GRAVITATIONAL_ACCELERATION, HEARTWOOD_RADIUS, RHO_WATER, MOLAR_GAS_CONSTANT, TEMPERATURE
from .tools.iotools import initialize_netcdf, write_netcdf, tree_properties_to_dict
from .model_variables import all_variables
from .odefun import odefun
from netCDF4 import Dataset


class Model:
    def __init__(self, tree: Tree, outputfile: str = "a.nc"):
        self.tree: Tree = tree
        self.outputfile: str = outputfile
        self.ncf: Dataset = initialize_netcdf(self, all_variables)

    def axial_fluxes(self) -> np.ndarray:
        """ calculates axial change in sap mass of every element

        variable names in this function's scope are the same as in Hölttä et al. (2006)
        shorter variable names are used to make errror prone calculations more readable
        """

        pressures: np.ndarray = self.tree.pressure
        k: np.ndarray = self.tree.axial_permeability
        eta: np.ndarray = self.tree.viscosity
        length: np.ndarray = self.tree.element_height
        E: np.ndarray = self.tree.transpiration_rate
        # calculate transport coefficients
        # TODO: add calculation for phloem sap density
        transport_ax: np.ndarray = k/eta/length*RHO_WATER * np.concatenate([self.tree.element_area([], 0),
                                                                           self.tree.element_area([], 1)], axis=1)

        # calculate upward fluxes
        Q_ax_down: np.ndarray = np.zeros((self.tree.num_elements, pressures.shape[1]))
        Q_ax_down[0:-1, :] = (np.diff(pressures, axis=0)
                              - RHO_WATER*GRAVITATIONAL_ACCELERATION*self.tree.element_height[0:-1].repeat(2, axis=1)
                              ) * transport_ax[0:-1, :]
        Q_ax_down[-1, 0] = (self.tree.ground_water_potential-pressures[-1, 0])*transport_ax[-1, 0]
        Q_ax_down[-1, 1] = 0  # no flux from phloem to soil

        Q_ax_up: np.ndarray = np.zeros((self.tree.num_elements, pressures.shape[1]))
        # TODO: Think if there is a better way to achieve this without flipping twice
        Q_ax_up[1:, :] = (np.flip(
            np.diff(
                np.flip(
                    pressures, axis=0), axis=0), axis=0)
            + RHO_WATER*GRAVITATIONAL_ACCELERATION*self.tree.element_height[1:].repeat(2, axis=1)
        ) * transport_ax[1:, :]
        Q_ax_up[0, :] = 0  # the upward flux is handled in transpiration rate for the highest element

        Q_ax: np.ndarray = Q_ax_up + Q_ax_down
        Q_ax[:, 0] = Q_ax[:, 0] - E.reshape(40,)  # subtract transpiration
        return Q_ax

    def radial_fluxes(self) -> np.ndarray:
        """ calculates radial change in sap mass of every element

        variable names in this function's scope are the same as in Hölttä et al. (2006)
        shorter variable names are used to make errror prone calculations more readable
        """

        pressures: np.ndarray = self.tree.pressure
        Lr: np.ndarray = self.tree.radial_hydraulic_conductivity
        C: np.ndarray = self.tree.sugar_concentration_as_numpy_array()
        Q_rad_phloem: np.ndarray = Lr.reshape((40, 1))*self.tree.element_height.reshape((40, 1))\
            * 2.0*math.pi*(self.tree.element_radius[:, 0]+HEARTWOOD_RADIUS).reshape(40, 1)\
            * RHO_WATER * (
            np.diff(np.flip(pressures, axis=1), axis=1) + C.reshape((40, 1))*MOLAR_GAS_CONSTANT*TEMPERATURE)

        Q_rad_xylem: np.ndarray = -Q_rad_phloem

        return np.concatenate((Q_rad_xylem, Q_rad_phloem), axis=1)

    def run(self, time_start: float = 1e-3, time_end: float = 120.0, dt: float = 0.01, output_interval: float = 60):
        # TODO: do not use explicit euler method
        for (ind, time) in enumerate(np.linspace(time_start, time_end, int((time_end-time_start)/dt))):
            # get the change in every elements mass
            dmdt_ax: np.ndarray = self.axial_fluxes()
            dmdt_rad: np.ndarray = self.radial_fluxes()
            if(np.abs(time % output_interval) < 1e-2):
                print(datetime.datetime.now(), "\t", time)
                results = tree_properties_to_dict(self.tree)
                results['dqrad'] = dmdt_rad
                results['dqax'] = dmdt_ax
                results['simulation_time'] = time
                results['model_index'] = ind
                write_netcdf(self.ncf, results)

            self.tree.pressure += dt*self.tree.elastic_modulus/(np.transpose(
                np.array([self.tree.element_volume([], 0),
                          self.tree.element_volume([], 1)])) * RHO_WATER)\
                * (dmdt_ax + dmdt_rad)

            for i in range(self.tree.num_elements):
                self.tree.solutes[i, 1].concentration += dt / self.tree.element_volume([i], 1)\
                    * (dmdt_ax[i, 1]
                       * self.tree.solutes[i, 1].concentration/RHO_WATER
                       + self.tree.sugar_loading_rate[i]
                       - self.tree.sugar_unloading_rate[i])
                if(i > self.tree.num_elements - 6):
                    self.tree.sugar_unloading_rate[i] = (self.tree.solutes[i, 1].concentration
                                                         - self.tree.sugar_target_concentration)*1e-5
            self.tree.element_radius += (dmdt_ax + dmdt_rad)\
                / (math.pi*RHO_WATER * np.repeat(self.tree.element_height, 2, axis=1) * self.tree.element_radius)

            # update sap viscosity
            self.tree.update_sap_viscosity()

    def run_scipy(self, time_start: float = 1e-3, time_end: float = 120.0, ind: int = 0):
        # Initial values from model.tree
        initial_values = [self.tree.pressure,
                          self.tree.sugar_concentration_as_numpy_array().reshape(40, 1),
                          self.tree.element_radius]
        # broadcast initial values into 1D array
        yinit = np.concatenate([initial_values[0].reshape(self.tree.num_elements*2, order='F'),
                                initial_values[1].reshape(self.tree.num_elements, order='F'),
                                initial_values[2].reshape(self.tree.num_elements*2, order='F')])

        sol = solve_ivp(lambda t, y: odefun(t, y, self), (time_start, time_end), yinit, method='BDF',
                        rtol=1e-3, atol=1e-1)
        # save the tree status
        print(datetime.datetime.now(), "\t", time_end)
        results = tree_properties_to_dict(self.tree)
        results['simulation_time'] = time_start
        results['model_index'] = ind
        results['dqrad'] = self.radial_fluxes()
        results['dqax'] = self.axial_fluxes()
        write_netcdf(self.ncf, results)
