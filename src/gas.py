import numpy as np
from typing import List


class Gas:
    """
        Description of Gas

        Attributes:
            nr (type):
            na (type):
            element_radius (type):
            na, (type):
            nr) (type):
            element_height (type):
            na, (type):
            nr) (type):
            diffusion_coefficients (type):
            na, (type):
            nr) (type):
            velocity (type):
            na, (type):
            nr) (type):
            space_division (type):
            na, (type):
            nr, (type):
            concentration (type):
            na, (type):
            nr) (type):
            kh (type):
            temperature (type):
            ambient_concentration (type):

        Args:
            num_radial_elements (int):
            num_axial_elements (int):
            elemenent_radius (List[List[float]],element_height:List[List[float]]):
            diffusion_coefficients (List[List[float]],velocity:List[List[float]]):
            space_division (List[List[List[float]]],concentration:List[List[float]]):
            henrys_law_coefficient (float):
            temperature (float):
            ambient_concentration (float):

        """

    def __init__(self, num_radial_elements: int,
                 num_axial_elements: int,
                 elemenent_radius: List[List[float]],
                 element_height: List[List[float]],
                 diffusion_coefficients: List[List[float]],
                 equilibration_rate: float,
                 velocity: List[List[float]],
                 space_division: List[List[List[float]]],
                 concentration: List[List[float]],
                 henrys_law_coefficient: float,
                 temperature: float,
                 ambient_concentration: float):
        self.nr: float = num_radial_elements
        self.na: float = num_axial_elements
        self.element_radius: np.ndarray = np.asarray(elemenent_radius).reshape(self.na, self.nr)
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
        self.flux_out: np.ndarray = np.zeros((self.nr, self.na))

    def radial_fluxes(self):
        transport_coefficient: np.ndarray = 2.0*np.pi*self.element_height*self.diffusion_coefficients
        r = self.radius_from_pith()

        # initialize flux matrices
        Q_rad = np.zeros((self.na, self.nr))

        Q_out = np.zeros((self.na, self.nr))

        Q_in = np.zeros((self.na, self.nr))

        Q_out[:, :-1] = -1.0*transport_coefficient[:, :-1]*np.diff(self.concentration[:, :, 1], axis=1)\
            / np.log(np.true_divide(r[:, 1:], r[:, :-1]))
        Q_out[:, -1] = -1.0*transport_coefficient[:, -1]*(self.concentration[:, -1, 1]-self.ambient_concentration)\
            / np.log((self.r[:, -1] + 0.5 * (r[:, -1] - r[:, -2]))/r[:, -1])

        Q_in[:, 1:] = transport_coefficient[:, 1:]*np.diff(self.concentration[:, :, 1], axis=1)\
            / np.log(np.true_divide(r[:, 1:], r[:, :-1]))

        Q_rad = Q_out+Q_in

        return Q_rad

    def axial_fluxes(self):
        A = self.head_area()
        # initialize flux matrices
        Q_ax = np.zeros((self.na, self.nr))
        Q_ax_top = np.zeros((self.na, self.nr))
        Q_ax_bottom = np.zeros((self.na, self.nr))

        Q_ax_bottom[1:, :] = -1.0*self.velocity[1:, :]*A[0:-1, :]*np.diff(self.concentration[:, :, 0])
        Q_ax_top[:-1, :] = self.velocity[:-1, :] * A[1:, :]*np.diff(self.concentration[:, :, 0])
        Q_ax = Q_ax_top + Q_ax_bottom
        return Q_ax

    def air_water_fluxes(self):
        """ calculates the equilibration between air and water phase gas concentration at every level.
        """
        volume_air = self.space_division[:, :, 0]*self.element_volume()
        Q_air_water = np.zeros((self.na, self.nr, 2))
        equilibrium_flux = np.diff(self.concentration, axis=2)*self.equilibration_rate*volume_air()*volume_air
        Q_air_water[:, :, 0] = -equilibrium_flux
        Q_air_water[:, :, 1] = equilibrium_flux

        return Q_air_water

    def sources(self):
        """ calculates the sources in the tree """
        return 0

    def sinks(self):
        """ calculates the sinks in the tree """
        return 0

    def radius_from_pith(self):
        return np.cumsum(self.element_radius, axis=1)

    def head_area(self):
        r = self.radius_from_pith()
        distance_sq = r**2
        distance_sq[:, 1:] = distance_sq[:, 1:]-distance_sq[:, :-1]

        return np.pi*distance_sq

    def element_volume(self):
        return self.head_area()*self.element_height
