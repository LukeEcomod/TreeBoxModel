import numpy as np
from typing import List


class Gas:
    """ Model for gas transport calculations in the tree model.

    Args:

    Attributes:

    """

    def __init__(self, num_radial_elements: int,
                 num_axial_elements: int,
                 elemenent_radius: List[List[float]],
                 element_height: List[List[float]],
                 diffusion_coefficients: List[List[float]],
                 velocity: List[List[float]],
                 space_division: List[List[List[float]]],
                 concentration: List[List[float]],
                 henrys_law_coefficient: float,
                 temperature: float,
                 ambient_concentration: float):

        self.nr = num_radial_elements
        self.na = num_axial_elements
        self.element_radius = np.asarray(elemenent_radius).reshape(self.na, self.nr)
        self.element_height = np.asarray(element_height).reshape(self.na, self.nr)
        self.diffusion_coefficients = np.asarray(diffusion_coefficients).reshape(self.na, self.nr)
        self.velocity = np.asarray(velocity).reshape(self.na, self.nr)
        # space division marks the fraction between air, water and cell wall in this order
        self.space_division = np.asarray(space_division).reshape(self.na, self.nr, 3)
        self.concentration = np.asarray(concentration).reshape(self.na, self.nr)
        self.kh = henrys_law_coefficient
        self.temperature = temperature
        self.ambient_concentration = ambient_concentration

    def radial_fluxes(self):
        transport_coefficient: np.ndarray = 2.0*np.pi*self.element_height*self.diffusion_coefficients
        r = self.radius_from_pith()

        # initialize flux matrices
        Q_rad = np.zeros((self.na, self.nr))

        Q_out = np.zeros((self.na, self.nr))

        Q_in = np.zeros((self.na, self.nr))

        Q_out[:, :-1] = -1.0*transport_coefficient[:, :-1]*np.diff(self.concentration, axis=1)\
            / np.log(np.true_divide(r[:, 1:], r[:, :-1]))
        Q_out[:, -1] = -1.0*transport_coefficient[:, -1]*(self.concentration[:, -1]-self.ambient_concentration)\
            / np.log((self.r[:, -1] + 0.5 * (r[:, -1] - r[:, -2]))/r[:, -1])

        Q_in[:, 1:] = transport_coefficient[:, 1:]*np.diff(self.concentration, axis=1)\
            / np.log(np.true_divide(r[:, 1:], r[:, :-1]))

        Q_rad = Q_out+Q_in

        return Q_rad

    def axial_fluxes(self):
        A = self.head_area()
        # initialize flux matrices
        Q_ax = np.zeros((self.na, self.nr))
        Q_ax_top = np.zeros((self.na, self.nr))
        Q_ax_bottom = np.zeros((self.na, self.nr))

        Q_ax_bottom[:, 1:] = -1.0*self.velocity[:, 1:]*A[:, 0:-1]*np.diff(self.concentration)
        Q_ax_top[:, :-1] = self.velocity[:, :-1] * A[:, 1:]*np.diff(self.concentration)
        Q_ax = Q_ax_top + Q_ax_bottom
        return Q_ax

    def radius_from_pith(self):
        return np.cumsum(self.element_radius, axis=1)

    def head_area(self):
        r = self.radius_from_pith()
        distance_sq = r**2
        distance_sq[:, 1:] = distance_sq[:, 1:]-distance_sq[:, :-1]

        return np.pi*distance_sq
