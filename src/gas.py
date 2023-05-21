from typing import List, Dict
import numpy as np
from scipy.integrate import solve_ivp

from .tools.iotools import compound_properties_to_dict, initialize_netcdf, write_netcdf
from .odefun_gas import odefun_gas
from .model_variables import cmp_variables


class Gas:
    ''' Defines a Gas class for calculating axial advection, diffusion and radial diffusion of a single compound in a cylindrical tree.
        List[float], List[List[float]] and List[List[List[float]]] argumnets are passed through numpy.asarray() function i.e., they can be also given as appropriate
        numpy arrays. List[float] argumetns need to have shape (num_axial_elements,), List[List[float]] shape (num_axial_elements, num_radial_elements)
        and List[List[List[float]]] either shape (num_axial_elements, num_radial_elements,3) for space_division or (num_axial_elements, num_radial_elements,2) for concentration
        where the last dimensions index 0 marks gaseous phase and index 1 aqueous phase and index 2 cell (in space division)
        Args:
            num_radial_elements (int): Number of radial elements.
            num_axial_elements (int): Number of axial elements.
            element_length (List[List[float]]): Length of each radial element (:math:`m`). The element length is determined by
                :math:`\\text{element length} = \\frac{\\text{tree radius}}{\\text{number of radial elements}}`.
            element_height (List[List[float]]): Height of each vertical element (:math:`m`). The element height is determined by
                :math:`\\text{element height} = \\frac{\\text{tree height}}{\\text{number of axial elements}}`.
            radial_diffusion_coefficient (List[List[float]]): Gas phase radial diffusion coefficient in each element (:math:`\\frac{m^2}{s}`).
            equilibration_rate (float): Rate which gaseous and aqueous phase equlibrate according to Henry's law (:math:`s^{-1}`).
            velocity (List[List[float]]): Sap flow velocity in each element (:math:`\\frac{m}{s}`)
            space_division (List[List[List[float]]]): Volume fraction of gas, aqeous phases and cell in each element.
            concentration (List[List[List[float]]]): Concentration of the compound (:math:`\\frac{mol}{m^3}`)
            henrys_law_coefficient (float): Coefficient of Henry's law (:math:`\frac{mol_{\\text{aq}}{mol_{\\{gas}}}`)
            temperature (float): Temperature in each element (:math:`K`). The compound is assumed to be in thermal equilibrium with the element.
            ambient_concentration (List[float]): Concentration of the compound in the gas phase outside the tree (:math:`\\frac{mol}{m^3}`).
            sources_and_sinks (List[float]): Source and/or sink strenghts in each layer (:math:`\\frac{mol}{s}`).
            soil_compound_aq_concentration (List[List[float]]): Concentration of the compound in the soil in the aqueous phase (:math:`\\frac{mol}{m^3}`).
            axial_water_uptake (List[float]): Uptake rate of water (:math:`\\frac{m^3}{s}`). This is ment to be set nonzero for the root elements.
            axial_diffusion_coefficient (List[List[float]]): Gas phase diffusion coefficient in the axial direction for each element (:math:`\\frac{m^2}{s}`).
            outputfile (str, optional): Name of the outputfile. Defaults to None. If None no file will be saved.
            max_output_lines (int, optional): How many lines to write for every call of run-method. Defaults to 100.

        Attributes:
            self.nr (float): Number of radial elements.
            self.na (float): Number of axial elements.
            self.element_length (np.ndarray): Length of each radial element (:math:`m`).
            self.element_height (np.ndarray): Height of each vertical element (:math:`m`).
            self.radial_diffusion_coefficient (np.ndarray): Gas phase radial diffusion coefficient in each element (:math:`\\frac{m^2}{s}`).
            self.axial_diffusion_coefficinet (np.ndarray): Gas phase diffusion coefficient in the axial direction for each element (:math:`\\frac{m^2}{s}`).
            self.equilibration_rate (float): Rate which gaseous and aqueous phase equlibrate according to Henry's law (:math:`s^{-1}`).
            self.velocity (np.ndarray): Sap flow velocity in each element (:math:`\\frac{m}{s}`).
            self.space_division (np.ndarray): Volume fraction of gas, aqeous phases and cell in each element.
            self.concentration (np.ndarray): Concentration of the compound (:math:`\\frac{mol}{m^3}`).
            self.kh (float): Coefficient of Henry's law (:math:`\frac{mol_{\\text{aq}}{mol_{\\{gas}}}`).
            self.temperature (float): Temperature in each element (:math:`K`).
            self.ambient_concentration (np.ndarray): Concentration of the compound in the gas phase outside the tree (:math:`\\frac{mol}{m^3}`).
            self.n_out (np.ndarray): Amount of the compound that has diffused out of the system in each vertical layer (:math:`mol`)
            self.n_in  (np.ndarray): Amount of the compound in each element (:math:`mol`)
            self.sources_and_sinks (np.ndarray): Source and/or sink strenghts in each layer (:math:`\\frac{mol}{s}`).
            self.radius_from_pith (np.ndarray): Radius from the pith to the end of the element (:math:`m`)
            self.radius_mid_point (np.ndarray): Radius from the pit to the middle of the leement (:math:`m`)
            self.head_area (np.ndarray): Top surface area of the element (:math:`m^2`)
            self.head_area_water (np.ndarray): Top surface area of the aqueous phase of the element (:math:`m^2`)
            self.element_volume (np.ndarray): Volume of each element (:math:`m^3`)
            self.element_volume_air (np.ndarray): Volume of the gaseous phase in each element (:math:`m^3`)
            self.element_volume_water (np.ndarray): Volume of the aqueous phase in each element (:math:`m^3`)
            self.element_volume_cell (np.ndarray): Cell volume in each element (:math:`m^3`)
            self.soil_compound_aq_concentration (np.ndarray): Concentration of the compound in the soil in the aqueous phase (:math:`\\frac{mol}{m^3}`).
                This is ment to be nonzero for the root elements.
            self.axial_water_uptake (np.ndarray): Uptake rate of water (:math:`\\frac{m^3}{s}`). This is ment to be nonzero for the root elements.
            self.outputfile (str): Name of the outputfile. Defaults to None. If None no file will be saved.
            self.max_output_lines (int): How many lines to write for every call of run-method. Defaults to 100.
    '''

    def __init__(self, num_radial_elements: int,
                 num_axial_elements: int,
                 element_length: List[List[float]],
                 element_height: List[List[float]],
                 radial_diffusion_coefficient: List[List[float]],
                 equilibration_rate: float,
                 velocity: List[List[float]],
                 space_division: List[List[List[float]]],
                 concentration: List[List[List[float]]],
                 henrys_law_coefficient: float,
                 temperature: float,
                 ambient_concentration: List[float],
                 sources_and_sinks: List[float],
                 soil_compound_aq_concentration: List[List[float]],
                 axial_water_uptake: List[float],
                 axial_diffusion_coefficient: List[List[float]] = None,
                 outputfile: str = None,
                 max_output_lines: int = 100
                 ):

        self.nr: float = num_radial_elements
        self.na: float = num_axial_elements
        self.element_length: np.ndarray = np.asarray(element_length).reshape(self.na, self.nr)
        self.element_height: np.ndarray = np.asarray(element_height).reshape(self.na, self.nr)
        self.radial_diffusion_coefficient: np.ndarray = np.asarray(radial_diffusion_coefficient).reshape(self.na, self.nr)
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
        # soil compound concentration includes gas and aq phase soil concentration in every axial element
        self.soil_compound_aq_concentration = np.asarray(soil_compound_aq_concentration).reshape(self.na, 1)
        self.axial_water_uptake = np.asarray(axial_water_uptake).reshape(self.na, 1)
        self.outputfile: str = outputfile
        self.max_output_lines: int = max_output_lines

        # set radii and midpoint lengths, head A and volume for the tree. Should speed long computations
        # TODO: define methods to update these if tree dimensions change
        self.radius_from_pith: np.ndarray = np.cumsum(self.element_length, axis=1)
        r_to_end: np.ndarray = np.concatenate((np.zeros((self.na, 1)), self.element_length), axis=1)
        cumulative_sum: np.ndarray = np.cumsum(r_to_end, axis=1)
        self.radius_mid_point = np.diff(cumulative_sum)/2.0 + cumulative_sum[:, :-1]

        self.length_from_bottom = np.flip(np.cumsum(self.element_height, axis=0))
        self.length_to_end: np.ndarray = np.concatenate((self.length_from_bottom, np.zeros((1, self.nr))), axis=0)

        self.length_mid_point = np.diff(self.length_to_end, axis=0)/2 + self.length_to_end[:-1, :]

        self.head_area: np.ndarray = self.radius_from_pith**2.0
        self.head_area[:, 1:] = self.head_area[:, 1:] - self.head_area[:, :-1]
        self.head_area = self.head_area * np.pi
        self.head_area_water = self.head_area * space_division[:, :, 1]
        self.element_volume: np.ndarray = self.head_area * self.element_height
        self.element_volume_air = self.space_division[:, :, 0]*self.element_volume
        self.element_volume_water = self.space_division[:, :, 1]*self.element_volume
        self.element_volume_cell = self.space_division[:, :, 2]*self.element_volume

        if axial_diffusion_coefficient is not None:
            self.axial_diffusion_coefficient = axial_diffusion_coefficient
        else:
            self.axial_diffusion_coefficient = np.zeros((self.na, self.nr))

        # Set number of moles of compound in gas and aq phases
        self.n_in = np.zeros_like(self.concentration)
        self.n_in[:, :, 0] = self.concentration[:, :, 0]*self.element_volume_air
        self.n_in[:, :, 1] = self.concentration[:, :, 1]*self.element_volume_water

        # Initialize a netcdf file for saving results
        if self.outputfile is not None:
            dims: Dict = {"axial_layers": self.na, "radial_layers": self.nr, "space_layers": 2}
            self.ncf = initialize_netcdf(outputfile, dims, cmp_variables)

    def radial_surface_area(self):
        """ Returns the radial surface A through which molecules diffuse in radial direction

        Returns:
            np.ndarray: Radial surface area of each element (:math:`m^2`)
        """
        r = self.radius_from_pith
        h = self.element_height
        return 2*np.pi*r*h

    def gas_radial_diffusion(self):
        """ Calculates radial diffusion flux of the compound in the gas phase

        Returns:
            np.ndarray: Radial diffusion flux in each element (:math:`\\frac{mol}{s}`)
        """
        A = self.radial_surface_area()
        dr = np.diff(self.radius_mid_point)
        D = self.radial_diffusion_coefficient
        dr_atm = (self.radius_from_pith[:, -1] - self.radius_mid_point[:, -1])

        # calculate concentration

        c_gas = self.n_in[:, :, 0]/self.element_volume_air

        # initialize flux matrices
        Q_rad = np.zeros((self.na, self.nr))

        Q_out = np.zeros((self.na, self.nr))

        Q_in = np.zeros((self.na, self.nr))

        # Calculate flux
        Q_out[:, :-1] = D[:, :-1] * A[:, :-1] / dr * np.diff(c_gas, axis=1)
        Q_out[:, -1] = D[:, -1] * A[:, -1] / dr_atm \
            * (self.ambient_concentration.reshape(self.na,) - c_gas[:, -1])

        Q_in[:, 1:] = -1.0*D[:, 1:] * A[:, :-1] / dr * np.diff(c_gas, axis=1)
        Q_rad = Q_out+Q_in

        return Q_rad, Q_out, Q_in

    def aq_axial_advection(self):
        """Calculates axial advection flux of the compound in the aqueous phase

        Returns:
            np.ndarray: Axial advection flux in the aqueous phase in each element (:math:`\\frac{mol}{s}`)
        """

        # Initialize flux matrices
        Q_ax_top_pos = np.zeros((self.na, self.nr))
        Q_ax_bottom_pos = np.zeros((self.na, self.nr))
        Q_ax_top_neg = np.zeros((self.na, self.nr))
        Q_ax_bottom_neg = np.zeros((self.na, self.nr))

        # Convenience determination of variables. Can be removed for speed up
        v = self.velocity
        A = self.head_area

        # Define concentration
        c_aq = self.n_in[:, :, 1] / self.element_volume_water

        # Calculate fluxs
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

    def gas_axial_diffusion(self):
        """ Calculates axial diffusion of the compound in the gas phase.

        Returns:
            np.ndarray: Axial diffusion flux in the gas phase in each element (:math:`\\frac{mol}{s}`)
        """
        # Initialize flux matrices
        Q_ax_top = np.zeros((self.na, self.nr))
        Q_ax_bottom = np.zeros((self.na, self.nr))
        Q_ax = np.zeros((self.na, self.nr))

        # Convienence variables. Can be removed for speed up
        A = self.head_area


        # gas phase concentration
        c_gas = self.n_in[:, :, 0] / self.element_volume_air

        Q_ax_top[1:, :] = -1.0*self.axial_diffusion_coefficient[1:, :]*A[1:,:]*np.diff(c_gas, axis=0)/(-1.0*np.diff(self.length_mid_point, axis=0))
        Q_ax_bottom[:-1, :] = self.axial_diffusion_coefficient[1:, :]*A[1:, :]*np.diff(c_gas, axis=0)/(-1.0*np.diff(self.length_mid_point, axis=0))

        dl_top = (self.length_from_bottom[0, :] - self.length_mid_point[0, :])
        dl_bottom = (self.length_from_bottom[-1, :] - self.length_mid_point[-1, :])
        Q_ax_top[0, :] = A[0, :] * self.axial_diffusion_coefficient[0, :] * (self.ambient_concentration[0] - c_gas[0, :])/dl_top
        Q_ax_bottom[-1, :] = A[-1, :] * self.axial_diffusion_coefficient[-1, :] * (
            self.ambient_concentration[-1] - c_gas[-1, :]) /dl_bottom
        Q_ax = Q_ax_top+Q_ax_bottom
        return Q_ax, Q_ax_top, Q_ax_bottom

    def air_water_fluxes(self):
        """ Calculates the equilibration flux between gas and aqeous phases.

        Returns:
            np.ndarray: Equilibration flux between gas and aqueous phase (:math:`\\frac{mol}{s}`)
        """
        # Initialize flux matrix. In the last dimension index 0 is for gas phase and index 1 for aqueous phase
        Q_air_water = np.zeros((self.na, self.nr, 2))

        # Concentrations
        c_gas = self.n_in[:, :, 0] / self.element_volume_air
        c_aq = self.n_in[:, :, 1] / self.element_volume_water

        # Calculate fluxes
        Q_air_water[:, :, 0] = (c_aq - c_gas*self.kh)*self.element_volume_water*self.equilibration_rate
        Q_air_water[:, :, 1] = -1.0*Q_air_water[:, :, 0]  # What leaves from gas phase must enter to aqueous phase and vice versa

        return Q_air_water

    def axial_water_uptake_source(self):
        """ Calculates compound uptake or loss from water uptake in each axial element. This flux meant to be zero if the element is not a root element

        Returns:
             Axial diffusion flux in the gas phase in each element (:math:`\\frac{mol}{s}`)
        """
        R_root = np.zeros((self.na, self.nr, 2))
        c_aq = self.n_in[:, :, 1] / self.element_volume_water

        Q_compound_aq_concentration = self.axial_water_uptake * self.soil_compound_aq_concentration

        # Elements where compound travel from roots to soil need to be handeld differently.
        # In this case the correct concentration is not the soil_compound_aq_concentration but
        # the concentration in that root element
        neg_mask = np.where(self.axial_water_uptake < 0)[0]
        if neg_mask.shape[0] > 0:
            Q_compound_aq_concentration[neg_mask] = self.axial_water_uptake[neg_mask]\
                * c_aq[neg_mask, -1].reshape(len(neg_mask), 1)

        # Reshape and set to aq phase of outermost element of R_root
        R_root[:, -1, 1] = Q_compound_aq_concentration.reshape(self.na,)

        return R_root

    def run(self, time_start: float, time_end: float):
        """ Runs the gas transport object from time_start to time_end and save results starting from the current state of the instance.

        Args:
            time_start (float): Time where the ODE solution starts
            time_end (float): Time where the ODE solution ends

        Returns:
            Object: Solution object of scipy.interpolate.solve_ivp
        """

        yinit = np.concatenate((self.n_in.reshape(self.na * self.nr * 2, order='F'),
                                self.n_out.reshape(self.na, order='F')))
        if time_start == 0:
            time_start = 1e-10
        sol = solve_ivp(lambda t, y: odefun_gas(t, y, self), (time_start, time_end), yinit, method='RK45',
                        rtol=1e-4, atol=1e-11)
        n_in = sol.y[:, -1]
        self.n_in = n_in[:self.na*self.nr*2].reshape(self.na, self.nr, 2, order='F')

        self.n_out = n_in[self.na*self.nr*2:].reshape(self.na, 1, order='F')
        res = sol.y
        time = sol.t
        if self.outputfile is not None:

            results = compound_properties_to_dict(self)
            # trim the solution for saving if the length of solution points exceeds maximum
            if self.max_output_lines == 1 or res.shape[1] <= self.max_output_lines:
                res = res[:, -1]
                time = time[-1]
                results['cmp_moles_in'] = res[:self.na*self.nr*2].reshape(self.na, self.nr, 2, order='F')
                results['cmp_moles_out'] = res[self.na*self.nr*2:].reshape(self.na, 1, order='F')
                results['cmp_simulation_time'] = time
                results['cmp_gas_radial_diffusion_flux'], _, _ = self.gas_radial_diffusion()
                results['cmp_gas_axial_diffusion_flux'], _, _ = self.gas_axial_diffusion()
                results['cmp_aq_axial_advection_flux'] = self.aq_axial_advection()
                results['cmp_eq_flux'] = self.air_water_fluxes()

                c_gas = results['cmp_moles_in'][:, :, 0] / self.element_volume_air
                c_aq = results['cmp_moles_in'][:, :, 1] / self.element_volume_water

                results['cmp_concentration'] = np.stack((c_gas, c_aq), axis=2)
                write_netcdf(self.ncf, results)
            elif res.shape[1] > self.max_output_lines:
                ind = np.concatenate(([0],
                                      np.arange(1, res.shape[1]-2,
                                                int(np.round((res.shape[1]-2)/(self.max_output_lines-2)))),
                                      [res.shape[1]-1]))
                res = res[:, ind]
                time = time[ind]

                for (ind, line) in enumerate(np.rollaxis(res, 1)):
                    results['cmp_moles_in'] = line[:self.na*self.nr*2].reshape(self.na, self.nr, 2, order='F')
                    results['cmp_moles_out'] = line[self.na*self.nr*2:].reshape(self.na, 1, order='F')
                    results['cmp_simulation_time'] = time[ind]
                    results['cmp_gas_radial_diffusion_flux'], _, _ = self.gas_radial_diffusion()
                    results['cmp_gas_axial_diffusion_flux'] = self.gas_axial_diffusion()
                    results['cmp_aq_axial_advection_flux'] = self.aq_axial_advection()
                    results['cmp_eq_flux'] = self.air_water_fluxes()
                    c_gas = results['cmp_moles_in'][:, :, 0] / self.element_volume_air
                    c_aq = results['cmp_moles_in'][:, :, 1] / self.element_volume_water

                    results['cmp_concentration'] = np.stack((c_gas, c_aq), axis=2)
                    write_netcdf(self.ncf, results)

        return sol

    def run_euler(self, time_start: float, time_end: float, dt: float):
        """ Runs the gas transport object from time_start to time_end.
        NB! It is better to use the run method than run_euler method. The run_euler method is meant for debugging.

        Args:
            time_start (float): Time where the ODE solution starts (:math:`s`)
            time_end (float): Time where the ODE solution ends (:math:`s`)
            dt (float): time step (:math:`s`)
        """
        time = time_start
        while time < time_end:
            Q_ax = self.aq_axial_advection()
            Q_rad, Q_out, _ = self.gas_radial_diffusion()
            Q_air_water = self.air_water_fluxes()
            R = self.sources_and_sinks()

            dndt_air = (Q_rad + R[:, :, 0] + Q_air_water[:, :, 0])
            dndt_water = (Q_ax + R[:, :, 1] + Q_air_water[:, :, 1])
            dndt_out = -1.0*Q_out[:, -1]

            self.n_out = self.n_out + dndt_out*dt
            self.n_in[:, :, 0] = self.n_in[:, :, 0] + dndt_air*dt
            self.n_in[:, :, 1] = self.n_in[:, :, 1] + dndt_water*dt
            time = time + dt
