import numpy as np


def odefun_gas(t: float, y: np.ndarray, gas) -> np.ndarray:  # pylint: disable=unused-argument
    ''' The right hand sides of ODEs used in the gas module. '''
    gas.concentration = y[:gas.nr*gas.na*2].reshape(gas.na, gas.nr, 2, order='F')
    Q_ax = gas.axial_fluxes()
    Q_rad, Q_out, _ = gas.radial_fluxes()
    Q_air_water = gas.air_water_fluxes()
    R = gas.sources_and_sinks
    Q_ax_diffusion = gas.axial_diffusion()
    # Source from root water uptake
    R_root = gas.axial_water_uptake_source()
    dcdt_air = (Q_rad + R[:, :, 0] + R_root[:, :, 0] + Q_air_water[:, :, 0] + Q_ax_diffusion)/gas.element_volume_air
    dcdt_water = (Q_ax + R[:, :, 1] + R_root[:, :, 1] + Q_air_water[:, :, 1])/gas.element_volume_water
    dndt_out = -1.0 * np.minimum(0, Q_out[:, -1])
    dydt = np.concatenate([dcdt_air.reshape(gas.na*gas.nr, order='F'),
                           dcdt_water.reshape(gas.na*gas.nr, order='F'),
                           dndt_out.reshape(gas.na, order='F')])
    return dydt
