import numpy as np
from src.gas import Gas


def odefun(t: float, y: np.ndarray, gas: Gas) -> np.ndarray:
    # update gas parameters
    gas.concentration = y.reshape(gas.na, gas.nr, 2)

    Q_ax = gas.axial_fluxes()
    Q_rad = gas.radial_fluxes()
    Q_air_water = gas.air_water_fluxes()
    P = gas.sources()
    S = gas.sinks()

    volume_air = gas.space_division[:, :, 0]*gas.element_volume()
    volume_water_cell = np.sum(gas.space_division[:, :, 1:], axis=2)*gas.element_volume()

    dcdt_air = (Q_rad + P[:, :, 0] + S[:, :, 0] + Q_air_water[:, :, 0])/volume_air
    dcdt_water = (Q_ax + P[:, :, 1] + S[:, :, 1] + Q_air_water[:, :, 1])/volume_water_cell

    dydt = np.concatenate([dcdt_air.reshape(gas.na*gas.nr, order='F'),
                           dcdt_water.reshape(gas.na*gas.nr, order='F')])
    return dydt
