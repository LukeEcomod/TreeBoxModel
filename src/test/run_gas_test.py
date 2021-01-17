import numpy as np
from numpy.core.arrayprint import _get_format_function
from ..gas import Gas
import matplotlib.pyplot as plt

nr = 25
na = 25


def init_gas():
    element_radius = np.repeat(0.1/nr, na*nr).reshape(na, nr)
    element_height = np.repeat(2/na, na*nr).reshape(na, nr)
    D = np.repeat(1e-11, na*nr).reshape(na, nr)
    eq_rate = 0.005
    velocity = np.repeat(1e-4, na*nr).reshape(na, nr)
    velocity[:, :12] = 0.0
    air_fraction = np.repeat(0.29, na*nr).reshape(na, nr, 1)
    water_fraction = np.repeat(0.52, na*nr).reshape(na, nr, 1)
    cell_fraction = np.repeat(0.19, na*nr).reshape(na, nr, 1)
    space_division = np.concatenate((air_fraction, water_fraction, cell_fraction), axis=2)
    temperature = 298.15
    kh = 1.4e-5*temperature*8.3145  # Sander: Compilation of Henry’s law constants (version 4.0) forwater as solvent
    ambient_concentration = 8.8e-5
    air_concentration = np.repeat(ambient_concentration, repeats=na*nr).reshape(na, nr, 1)
    water_concentration = (kh*air_concentration).reshape(na, nr, 1)
    concentration = np.concatenate((air_concentration, water_concentration), axis=2)
    sources_and_sinks_func = source(1e-9)

    return Gas(num_radial_elements=nr,
               num_axial_elements=na,
               element_radius=element_radius,
               element_height=element_height,
               diffusion_coefficients=D,
               equilibration_rate=eq_rate,
               velocity=velocity,
               space_division=space_division,
               concentration=concentration,
               henrys_law_coefficient=kh,
               temperature=temperature,
               ambient_concentration=ambient_concentration,
               sources_and_sinks_func=sources_and_sinks_func)


def source(source_strength):
    R = np.zeros((na, nr, 2))
    R[-1, 11:, 1] = source_strength
    return R


if __name__ == '__main__':
    gas = init_gas()
    sol = gas.run(1e-12, 60)
    gg = sol.y[:na*nr*2, -1].reshape(na, nr, 2, order='F')
    nout = sol.y[na*nr*2:, -1].reshape(na, 1, order='F')
    r = gas.radius_from_pith()
    h = gas.element_height
    flux_out = nout/(2*np.pi*r[:, -1].reshape(na, 1)*h[:, -1].reshape(na, 1))*16.04  # g / m2 / h
    flux_out[flux_out < 1e-11] = 0.0
    print(flux_out)
    gg = sol.y[:na*nr*2, -1].reshape(na, nr, 2, order='F')
    print((gg > 0.0).all())
    print(np.where(gg < 0))
    print(gg[24, :, :])
    plt.plot(sol.t, sol.y[624, :]/8.8e-5)
    plt.show()
