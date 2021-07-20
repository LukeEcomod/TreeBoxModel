from typing import List, Callable, Dict
import numpy as np
from scipy.integrate import solve_ivp
from .tools.iotools import gas_properties_to_dict
from .odefun_gas import odefun_gas
from .tools.iotools import initialize_netcdf, write_netcdf
from .model_variables import gas_variables


class Gas:

    def __init__(self, num_radial_elements: int,
                 num_axial_elements: int,
                 element_radius: List[List[float]],
                 element_height: List[List[float]],
                 diffusion_coefficients: List[List[float]],
                 equilibration_rate: float,
                 velocity: List[List[float]],
                 space_division: List[List[List[float]]],
                 concentration: List[List[float]],
                 henrys_law_coefficient: float,
                 temperature: float,
                 ambient_concentration: float,
                 outputfile: str = '',
                 sources_and_sinks_func: Callable = None,
                 max_output_lines: int = 100):
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
        self.kh: float = henrys_law_coefficient  # unitless form kh = cair/cwater
        self.temperature: float = temperature
        self.ambient_concentration: float = ambient_concentration
        self.n_out: np.ndarray = np.zeros((self.na, 1))
        self.sources_and_sinks_func: Callable = sources_and_sinks_func
        self.outputfile: str = outputfile
        self.max_output_lines: int = max_output_lines

        # set head area and volume for the tree. Should speed long computations
        self.radius_from_pith: np.ndarray = np.cumsum(self.element_radius, axis=1)
        r_to_end: np.ndarray = np.concatenate((np.zeros((self.na, 1)), self.element_radius), axis=1)
        cumulative_sum: np.ndarray = np.cumsum(r_to_end, axis=1)
        self.radius_mid_point = np.diff(cumulative_sum/2) + cumulative_sum[:, :-1]
        self.head_area: np.ndarray = self.radius_from_pith**2.0
        self.head_area[:, 1:] = self.head_area[:, 1:] - self.head_area[:, :-1]
        self.head_area = self.head_area * np.pi
        self.element_volume: np.ndarray = self.head_area * self.element_height

        self.element_volume_air = self.space_division[:, :, 0]*self.element_volume
        self.element_volume_water_cell = np.sum(self.space_division[:, :, 1:], axis=2)*self.element_volume
        if len(outputfile) != 0:
            dims: Dict = {"axial_layers": self.na, "radial_layers": self.nr, "space_layers": 2}
            self.ncf = initialize_netcdf(outputfile, dims, gas_variables)

    def radial_fluxes(self):
        transport_coefficient: np.ndarray = 2.0*np.pi*self.element_height*self.diffusion_coefficients
        r = self.radius_mid_point
        # initialize flux matrices
        Q_rad = np.zeros((self.na, self.nr))

        Q_out = np.zeros((self.na, self.nr))

        Q_in = np.zeros((self.na, self.nr))

        Q_out[:, :-1] = transport_coefficient[:, :-1]*np.diff(self.concentration[:, :, 0], axis=1)\
            / np.log(np.true_divide(r[:, 1:], r[:, :-1]))
        Q_out[:, -1] = transport_coefficient[:, -1]*(self.ambient_concentration - self.concentration[:, -1, 0])\
            / np.log((r[:, -1] + 0.5 * (r[:, -1] - r[:, -2]))/r[:, -1])

        Q_in[:, 1:] = -1.0*transport_coefficient[:, 1:]*np.diff(self.concentration[:, :, 0], axis=1)\
            / np.log(np.true_divide(r[:, 1:], r[:, :-1]))

        Q_rad = Q_out+Q_in

        return Q_rad

    def axial_fluxes(self):
        # A = self.head_area()
        # initialize flux matrices
        Q_ax = np.zeros((self.na, self.nr))
        Q_ax_top = np.zeros((self.na, self.nr))
        Q_ax_bottom = np.zeros((self.na, self.nr))

        Q_ax_bottom[:-1, :] = self.velocity[:-1, :]*self.head_area[:-1, :]*np.diff(self.concentration[:, :, 1], axis=0)
        Q_ax_top[1:, :] = -1.0*self.velocity[1:, :] * self.head_area[1:, :]*np.diff(self.concentration[:, :, 1], axis=0)
        Q_ax = Q_ax_top + Q_ax_bottom
        return Q_ax

    def air_water_fluxes(self):
        """ calculates the equilibration between air and water phase gas concentration at every level.
        """
        # volume_air = self.space_division[:, :, 0]*self.element_volume()
        Q_air_water = np.zeros((self.na, self.nr, 2))
        conc = self.concentration.copy()

        conc[:, :, 1] = conc[:, :, 1]/self.kh
        equilibrium_flux = np.diff(conc, axis=2).reshape(self.na, self.nr)\
            * self.equilibration_rate*self.element_volume_air
        Q_air_water[:, :, 0] = equilibrium_flux
        Q_air_water[:, :, 1] = -equilibrium_flux

        return Q_air_water

    def sources_and_sinks(self):
        R = np.zeros((self.na, self.nr, 2))
        if self.sources_and_sinks_func is not None:
            R = self.sources_and_sinks_func
        return R

    # def radius_from_pith(self):
    #     return np.cumsum(self.element_radius, axis=1)

    # def head_area(self):
    #     # r = self.radius_from_pith()
    #     distance_sq = self.radius_area**2.0
    #     distance_sq[:, 1:] = distance_sq[:, 1:] - distance_sq[:, :-1]
    #     return np.pi*distance_sq

    # def element_volume(self):
    #     return self.head_area()*self.element_height

    def run(self, time_start: float, time_end: float):
        """ function to run the gas module independently. """
        yinit = np.concatenate((self.concentration.reshape(self.na * self.nr * 2, order='F'),
                                self.n_out.reshape(self.na, order='F')))
        if time_start == 0:
            time_start = 1e-10
        sol = solve_ivp(lambda t, y: odefun_gas(t, y, self), (time_start, time_end), yinit, method='DOP853',
                        rtol=1e-9, atol=1e-9)
        conc = sol.y[:, -1]
        self.concentration = conc[:self.na*self.nr*2].reshape(self.na, self.nr, 2, order='F')
        self.n_out = conc[self.na*self.nr*2:].reshape(self.na, 1, order='F')
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
                results['gas_radial_flux'] = self.radial_fluxes()
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
