from typing import Dict
import numpy as np
import math
import datetime
from .tree import Tree
from .constants import RHO_WATER, MOLAR_GAS_CONSTANT, TEMPERATURE
from .tools.iotools import initialize_netcdf, write_netcdf, tree_properties_to_dict
from .model_variables import all_variables


class Model:
    def __init__(self, tree: Tree, outputfile: str = "a.nc"):
        self.tree: Tree = tree
        self.outputfile: str = outputfile

    def axial_fluxes(self) -> np.ndarray:
        """ calculates axial change in sap mass of every element

        variable names in this function's scope are the same as in Hölttä et al. (2006)
        shorter variable names are used to make errror prone calculations more readable
        """
        # TODO:figure out why ndarray type hint does not raise error in pyright
        # FIXME: refactor into axial and radial fluxes
        # TODO: think is this the smartest way to do the calculations
        pressures: np.ndarray = self.tree.element_property_as_numpy_array('pressure')
        k: np.ndarray = self.tree.element_property_as_numpy_array('permeability')
        eta: np.ndarray = self.tree.element_property_as_numpy_array('viscosity')
        length: np.ndarray = self.tree.element_property_as_numpy_array('height')
        E: np.ndarray = self.tree.element_property_as_numpy_array('transpiration_rate')

        # calculate transport coefficients
        # TODO: add calculation for phloem sap density
        transport_ax: np.ndarray = k/eta/length*RHO_WATER*np.transpose(
            np.array([self.tree.element_area([], 0),
                      self.tree.element_area([], 1)]))
        # calculate upward fluxes
        Q_ax_down: np.ndarray = np.zeros((self.tree.num_elements, pressures.shape[1]))
        Q_ax_down[0:-1, :] = np.diff(pressures, axis=0)*transport_ax[0:-1, :]
        Q_ax_down[-1, 0] = (self.tree.ground_water_potential-pressures[-1, 0])*transport_ax[-1, 0]
        Q_ax_down[-1, 1] = 0  # no flux from phloem to soil

        Q_ax_up: np.ndarray = np.zeros((self.tree.num_elements, pressures.shape[1]))
        # TODO: Think if there is a better way to achieve this without flipping twice
        Q_ax_up[1:len(Q_ax_up), :] = np.flip(
            np.diff(
                np.flip(
                    pressures, axis=0), axis=0), axis=0) * transport_ax[1:len(Q_ax_up), :]

        Q_ax_up[0, :] = 0  # the upward flux is handled in transpiration rate for the highest element

        Q_ax: np.ndarray = Q_ax_up + Q_ax_down - E
        return Q_ax

    def radial_fluxes(self) -> np.ndarray:
        """ calculates radial change in sap mass of every element

        variable names in this function's scope are the same as in Hölttä et al. (2006)
        shorter variable names are used to make errror prone calculations more readable
        """

        pressures: np.ndarray = self.tree.element_property_as_numpy_array('pressure')
        Lr: np.ndarray = self.tree.element_property_as_numpy_array('hydraulic_conductivity')
        C: np.ndarray = self.tree.sugar_concentration_as_numpy_array()
        Q_rad_phloem: np.ndarray = Lr[:, 1].reshape((40, 1))*self.tree.cross_sectional_area().reshape((40, 1))\
            * RHO_WATER * (
            np.diff(np.flip(pressures, axis=1), axis=1) + C.reshape((40, 1))*MOLAR_GAS_CONSTANT*TEMPERATURE)

        Q_rad_xylem: np.ndarray = -Q_rad_phloem
        return np.concatenate((Q_rad_xylem, Q_rad_phloem), axis=1)

    def run(self, time_start: float = 1e-3, time_end: float = 120.0, dt: float = 0.01, output_interval: float = 60):
        # TODO: do not use explicit euler method
        # set up the netcdf writer
        ncf = initialize_netcdf(self, all_variables)
        for (ind, time) in enumerate(np.linspace(time_start, time_end, int(time_end/dt))):
            # get the change in every elements mass
            dmdt_ax: np.ndarray = self.axial_fluxes()
            dmdt_rad: np.ndarray = self.radial_fluxes()
            # print(dmdt_rad.shape)
            if(np.abs(time % output_interval) < 1e-2):
                print(datetime.datetime.now(), "\t", time)
                results = tree_properties_to_dict(self.tree)
                results['dqrad'] = dmdt_rad
                results['dqax'] = dmdt_ax
                results['simulation_time'] = time
                results['model_index'] = ind
                write_netcdf(ncf, results)
            for i in range(self.tree.num_elements):
                # TODO: do this without for loop
                # update pressures
                self.tree.elements[i][0].pressure += dt * self.tree.elements[i][0].elastic_modulus\
                    / float((self.tree.element_volume([i], 0)) * RHO_WATER) * (dmdt_ax[i, 0] + dmdt_rad[i, 0])
                self.tree.elements[i][1].pressure += dt * self.tree.elements[i][1].elastic_modulus\
                    / float((self.tree.element_volume([i], 1)) * RHO_WATER) * (dmdt_ax[i, 1] + dmdt_rad[i, 1])
                # update sugar concentration
                self.tree.elements[i][1].solutes[0].concentration += dt * dmdt_ax[i, 1]\
                    * self.tree.elements[i][1].solutes[0].concentration/RHO_WATER\
                    + self.tree.elements[i][1].sugar_loading_rate\
                    + self.tree.elements[i][1].sugar_unloading_rate

                # update radius of each element

                self.tree.elements[i][0].radius += dt * (dmdt_ax[i, 0] + dmdt_rad[i, 0])\
                    / (math.pi * RHO_WATER
                       * self.tree.elements[i][0].height
                       * self.tree.elements[i][0].radius)
                self.tree.elements[i][1].radius += dt * (dmdt_ax[i, 1] + dmdt_rad[i, 1])\
                    / (math.pi * RHO_WATER
                       * self.tree.elements[i][1].height
                       * self.tree.elements[i][1].radius)

            # update sap viscosity
            self.tree.update_sap_viscosity()
