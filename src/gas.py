from typing import List, Dict
import numpy as np
from scipy.integrate import solve_ivp

from .tools.iotools import gas_properties_to_dict
from .odefun_gas import odefun_gas
from .tools.iotools import initialize_netcdf, write_netcdf
from .model_variables import gas_variables


class Gas:
    ''' Defines a Gas class for calculating axial advection and radial diffusion of a self in a tree.'''

    def __init__(self, num_radial_elements: int,
                 num_axial_elements: int,
                 element_radius: List[List[float]],
                 element_height: List[List[float]],
                 diffusion_coefficients: List[List[float]],
                 equilibration_rate: float,
                 velocity: List[List[float]],
                 space_division: List[List[List[float]]],
                 concentration: List[List[List[float]]],
                 henrys_law_coefficient: float,
                 temperature: float,
                 ambient_concentration: List[float],
                 sources_and_sinks: List[float],
                 soil_gas_aq_concentration: List[List[float]],
                 axial_water_uptake: List[float],
                 axial_diffusion_coefficient: float = 0.0,
                 outputfile: str = '',
                 max_output_lines: int = 100
                 ):
        self.nr: float = num_radial_elements
        self.na: float = num_axial_elements
        self.element_radius: np.ndarray = np.asarray(element_radius).reshape(self.na, self.nr)
        self.element_height: np.ndarray = np.asarray(element_height).reshape(self.na, self.nr)
        self.diffusion_coefficients: np.ndarray = np.asarray(diffusion_coefficients).reshape(self.na, self.nr)
        self.equilibration_rate: float = equilibration_rate
        self.velocity: np.ndarray = np.asarray(velocity).reshape(self.na, self.nr)
        # space division marks the fraction between air, water and cell wall in this order
        self.space_division: np.ndarray = np.asarray(space_division).reshape(self.na, self.nr, 3)
        self.concentration: np.ndarray = np.asarray(concentration).reshape(self.na, self.nr, 2)
        self.kh: float = henrys_law_coefficient  # kh = cwater/cair [kh] = (m^3 water / m^3 air)
        self.temperature: float = temperature
        self.ambient_concentration: np.ndarray = np.asarray(ambient_concentration)
        self.n_out: np.ndarray = np.zeros((self.na, 1))
        self.sources_and_sinks: np.ndarray = sources_and_sinks.reshape(self.na, self.nr, 2)
        # soil gas concentration includes gas and aq phase soil concentration in every axial element
        self.soil_gas_aq_concentration = np.asarray(soil_gas_aq_concentration).reshape(self.na, 1)
        self.axial_water_uptake = np.asarray(axial_water_uptake).reshape(self.na, 1)
        self.outputfile: str = outputfile
        self.max_output_lines: int = max_output_lines

        # set head A and volume for the tree. Should speed long computations
        self.radius_from_pith: np.ndarray = np.cumsum(self.element_radius, axis=1)
        r_to_end: np.ndarray = np.concatenate((np.zeros((self.na, 1)), self.element_radius), axis=1)
        cumulative_sum: np.ndarray = np.cumsum(r_to_end, axis=1)
        self.radius_mid_point = np.diff(cumulative_sum/2) + cumulative_sum[:, :-1]
        self.head_area: np.ndarray = self.radius_from_pith**2.0
        self.head_area[:, 1:] = self.head_area[:, 1:] - self.head_area[:, :-1]
        self.head_area = self.head_area * np.pi
        self.head_area_water = self.head_area * space_division[:, :, 1]
        self.element_volume: np.ndarray = self.head_area * self.element_height
        self.element_volume_air = self.space_division[:, :, 0]*self.element_volume
        self.element_volume_water = self.space_division[:, :, 1]*self.element_volume
        self.element_volume_cell = self.space_division[:, :, 2]*self.element_volume
        self.D_ax = axial_diffusion_coefficient

        if len(outputfile) != 0:
            dims: Dict = {"axial_layers": self.na, "radial_layers": self.nr, "space_layers": 2}
            self.ncf = initialize_netcdf(outputfile, dims, gas_variables)

    def radial_surface_area(self):
        ''' Returns the radial surface A through which molecules diffuse in radial direction '''
        r = self.radius_from_pith
        h = self.element_height
        return 2*np.pi*r*h

    def radial_fluxes(self):
        ''' Calculates radial diffusion of the self.'''
        # transport_coefficient: np.ndarray = 2.0*np.pi*self.element_height*self.diffusion_coefficients
        A = self.radial_surface_area()
        dr = np.diff(self.radius_mid_point)
        D = self.diffusion_coefficients
        dr_atm = (self.radius_from_pith[:, -1] - self.radius_mid_point[:, -1])
        # initialize flux matrices
        Q_rad = np.zeros((self.na, self.nr))

        Q_out = np.zeros((self.na, self.nr))

        Q_in = np.zeros((self.na, self.nr))

        Q_out[:, :-1] = D[:, :-1] * A[:, :-1] / dr * np.diff(self.concentration[:, :, 0], axis=1)
        Q_out[:, -1] = D[:, -1] * A[:, -1] / dr_atm \
            * (self.ambient_concentration.reshape(self.na,) - self.concentration[:, -1, 0])

        # print('new line')
        # print(dr_atm[[88, 118]])
        # print(Q_out[88, -1])
        # print(Q_out[118, -1])
        # print(self.concentration[88, -1, 0])
        # print(self.concentration[118, -1, 0])
        #Q_out[-24: -1] = 0.0

        Q_in[:, 1:] = -1.0*D[:, 1:] * A[:, :-1] / dr * np.diff(self.concentration[:, :, 0], axis=1)
        Q_rad = Q_out+Q_in

        return Q_rad, Q_out, Q_in

    def axial_fluxes(self):
        ''' Calculates axial advection of the gas'''
        Q_ax_top_pos = np.zeros((self.na, self.nr))
        Q_ax_bottom_pos = np.zeros((self.na, self.nr))
        Q_ax_top_neg = np.zeros((self.na, self.nr))
        Q_ax_bottom_neg = np.zeros((self.na, self.nr))
        v = self.velocity
        A = self.head_area
        c_aq = self.concentration[:, :, 1]
        Q_ax_top_pos[1:, :] = -v[1:, :]*A[1:, :]*c_aq[1:, :]
        Q_ax_bottom_pos[:-1, :] = v[1:, :]*A[1:, :]*c_aq[1:, :]
        Q_ax_top_neg[1:, :] = -v[1:, :]*A[1:, :]*c_aq[:-1, :]
        Q_ax_bottom_neg[:-1, :] = v[1:, :]*A[1:, :]*c_aq[:-1, :]

        # Change amounts that are not real to zero

        # Q_ax_top_pos is always negative when velocity > 0
        Q_ax_top_pos[Q_ax_top_pos > 0.0] = 0.0

        # Q_ax_bottom_pos is alays positive when velocity > 0
        Q_ax_bottom_pos[Q_ax_bottom_pos < 0.0] = 0.0

        # Q_ax_top_neg is always positive when velocity < 0
        Q_ax_top_neg[Q_ax_top_neg < 0.0] = 0.0

        # Q_ax_bottom_neg is always negative when velocity < 0
        Q_ax_bottom_neg[Q_ax_bottom_neg > 0.0] = 0.0

        Q_ax = Q_ax_top_pos + Q_ax_bottom_pos + Q_ax_top_neg + Q_ax_bottom_neg

        return Q_ax

    def axial_diffusion(self):
        Q_ax_top = np.zeros((self.na, self.nr))
        Q_ax_bottom = np.zeros((self.na, self.nr))
        Q_ax = np.zeros((self.na, self.nr))
        A = self.head_area
        c_gas = self.concentration[:, :, 0]
        length: np.ndarray = self.element_height[:, 0]
        l: np.ndarray = np.concatenate(([0], length.reshape(len(length),)))
        cumulative_sum: np.ndarray = np.cumsum(l)
        l_midpoints = ((length/2) + cumulative_sum[:-1]).reshape(len(length), 1)
        l_midpoints = np.repeat(l_midpoints, repeats=self.nr, axis=1)
        dl_midpoints = np.diff(l_midpoints, axis=0)
        Q_ax_top[1:] = -1.0*self.D_ax*A[1:]*np.diff(c_gas, axis=0)/dl_midpoints
        Q_ax_bottom[:-1] = self.D_ax*A[1:]*np.diff(c_gas, axis=0)/dl_midpoints

        # calculate the bottom and top fluxes using the radial diffusion coefficient
        Q_ax_top[0, :] = A[0, :] * self.diffusion_coefficients[0, :] * (self.ambient_concentration[0] - c_gas[0, :])
        Q_ax_bottom[-1, :] = A[-1, :] * self.diffusion_coefficients[-1, :] * (
            self.ambient_concentration[-1] - c_gas[-1, :])
        Q_ax = Q_ax_top+Q_ax_bottom
        return Q_ax

    def air_water_fluxes(self):
        """ calculates the equilibration between air and water phase self concentration at every lev.
        """
        Q_air_water = np.zeros((self.na, self.nr, 2))
        c_aq = self.concentration.copy()

        c_air = c_aq[:, :, 0]
        c_water = c_aq[:, :, 1]

        Q_air_water[:, :, 0] = (c_water - c_air*self.kh)*self.element_volume_water*self.equilibration_rate
        Q_air_water[:, :, 1] = (c_air*self.kh - c_water)*self.element_volume_water*self.equilibration_rate

        return Q_air_water

    def axial_water_uptake_source(self):
        '''Calculates gas uptake or loss from water uptake or loss along each axial element.'''
        R_root = np.zeros((self.na, self.nr, 2))
        Q_gas_aq_root = self.axial_water_uptake * self.soil_gas_aq_concentration
        neg_mask, _ = np.where(Q_gas_aq_root < 0)
        Q_gas_aq_root[neg_mask] = self.axial_water_uptake[neg_mask]\
            * self.concentration[neg_mask, -1, 1].reshape(len(neg_mask), 1)
        R_root[:, -1, 1] = Q_gas_aq_root.reshape(self.na,)

        return R_root

    def run(self, time_start: float, time_end: float):
        """ function to run the self module independently. """

        yinit = np.concatenate((self.concentration.reshape(self.na * self.nr * 2, order='F'),
                                self.n_out.reshape(self.na, order='F')))
        if time_start == 0:
            time_start = 1e-10
        sol = solve_ivp(lambda t, y: odefun_gas(t, y, self), (time_start, time_end), yinit, method='RK45',
                        rtol=1e-6, atol=1e-35)
        c_aq = sol.y[:, -1]
        self.concentration = c_aq[:self.na*self.nr*2].reshape(self.na, self.nr, 2, order='F')

        self.n_out = c_aq[self.na*self.nr*2:].reshape(self.na, 1, order='F')
        res = sol.y
        time = sol.t
        if len(self.outputfile) != 0:

            results = gas_properties_to_dict(self)
            # trim the solution for saving if the length of solution points exceeds maximum
            if self.max_output_lines == 1:
                res = res[:, -1]
                time = time[-1]
                results['gas_concentration'] = res[:self.na*self.nr*2].reshape(self.na, self.nr, 2, order='F')
                results['gas_moles_out'] = res[self.na*self.nr*2:].reshape(self.na, 1, order='F')
                results['gas_simulation_time'] = time
                results['gas_radial_flux'], _, _ = self.radial_fluxes()
                results['gas_axial_flux'] = self.axial_fluxes()
                results['gas_eq_flux'] = self.air_water_fluxes()
                write_netcdf(self.ncf, results)
            elif res.shape[1] > self.max_output_lines:
                ind = np.concatenate(([0],
                                      np.arange(1, res.shape[1]-2,
                                                int(np.round((res.shape[1]-2)/(self.max_output_lines-2)))),
                                      [res.shape[1]-1]))
                res = res[:, ind]
                time = time[ind]

                for (ind, line) in enumerate(np.rollaxis(res, 1)):

                    results['gas_concentration'] = line[:self.na*self.nr*2].reshape(self.na, self.nr, 2, order='F')
                    results['gas_moles_out'] = line[self.na*self.nr*2:].reshape(self.na, 1, order='F')
                    results['gas_simulation_time'] = time[ind]
                    results['gas_radial_flux'] = self.radial_fluxes()
                    results['gas_axial_flux'] = self.axial_fluxes()
                    results['gas_eq_flux'] = self.air_water_fluxes()
                    write_netcdf(self.ncf, results)

        return sol

    def run_euler(self, time_start: float, time_end: float, dt: float):

        time = time_start
        while time < time_end:
            Q_ax = self.axial_fluxes()
            Q_rad, Q_out, _ = self.radial_fluxes()
            Q_air_water = self.air_water_fluxes()
            R = self.sources_and_sinks()

            dcdt_air = (Q_rad + R[:, :, 0] + Q_air_water[:, :, 0])/self.element_volume_air
            dcdt_water = (Q_ax + R[:, :, 1] + Q_air_water[:, :, 1])/self.element_volume_water
            dndt_out = -1.0*Q_out[:, -1]  # -1.0*np.minimum(0, Q_rad[:, -1])

            self.n_out = self.n_out + dndt_out*dt
            self.concentration[:, :, 0] = self.concentration[:, :, 0] + dcdt_air*dt
            self.concentration[:, :, 1] = self.concentration[:, :, 1] + dcdt_water*dt
            time = time + dt
