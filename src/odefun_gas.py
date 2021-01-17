import numpy as np


def odefun_gas(t: float, y: np.ndarray, gas) -> np.ndarray:
    # update gas parameters
    gas.concentration = y[:gas.nr*gas.na*2].reshape(gas.na, gas.nr, 2, order='F')
    Q_ax = gas.axial_fluxes()
    Q_rad = gas.radial_fluxes()
    Q_air_water = gas.air_water_fluxes()
    R = gas.sources_and_sinks()
    # print(R[:, :, 0])
    volume_air = gas.space_division[:, :, 0]*gas.element_volume()
    volume_water_cell = np.sum(gas.space_division[:, :, 1:], axis=2)*gas.element_volume()

    dcdt_air = (Q_rad + R[:, :, 0] + Q_air_water[:, :, 0])/volume_air
    dcdt_water = (Q_ax + R[:, :, 1] + Q_air_water[:, :, 1])/volume_water_cell
    dndt_out = -1.0*np.minimum(0, Q_rad[:, -1])
    # print(Q_rad[:, -1])
    dydt = np.concatenate([dcdt_air.reshape(gas.na*gas.nr, order='F'),
                           dcdt_water.reshape(gas.na*gas.nr, order='F'),
                           dndt_out.reshape(gas.na, order='F')])
    return dydt
