from src.constants import PHLOEM_RADIUS, XYLEM_RADIUS
from src.model import Model
from src.tree import Tree
from typing import List
import numpy as np
import math


if __name__ == "__main__":

    height: float = 12.0

    num_elements: int = 40

    radii: List[float] = [XYLEM_RADIUS, PHLOEM_RADIUS]
    transpiration_profile: List[float] = [0.0 for i in range(num_elements)]
    # transpiration_profile[0:10] = [0.9*1e-3/10.0]*10  # kg/s
    transpiration_max = 0.9e-3/10.0  # kg/s
    daytime = np.linspace(0, 12, 12*60+1)  # hours starting from 06:00 AM, 5 minute interval
    transpiration = np.sin(daytime*math.pi/12.0)*transpiration_max
    photosynth_profile: List[float] = [0 for i in range(num_elements)]
    photosynth_profile[0:5] = [4e-6, 4e-6, 4e-6, 4e-6, 4e-6]

    sugar_profile: List[float] = [10.0 for i in range(num_elements)]
    sugar_loading_profile = photosynth_profile

    sugar_unloading_profile: List[float] = [0.0 for i in range(num_elements)]

    axial_permeability_profile: List[List[float]] = [[1.5e-12, 6.0e-12]]*num_elements

    radial_hydr_conductivity: List[float] = [1e-13] * num_elements

    elastic_modulus_profile: List[List[float]] = [[1000e6, 30e6]] * num_elements

    ground_water_potential: float = 0.0

    tree = Tree(height=height,
                num_elements=num_elements,
                initial_radius=radii,
                transpiration_profile=transpiration_profile,
                photosynthesis_profile=photosynth_profile,
                sugar_profile=sugar_profile,
                sugar_loading_profile=sugar_loading_profile,
                sugar_unloading_profile=sugar_unloading_profile,
                axial_permeability_profile=axial_permeability_profile,
                radial_hydraulic_conductivity_profile=radial_hydr_conductivity,
                elastic_modulus_profile=elastic_modulus_profile,
                ground_water_potential=ground_water_potential)

    model = Model(tree, outputfile="test_qrad_24h.nc")

    for (ind, time) in enumerate(daytime):
        # set new transpiration rate
        transpiration_rate = transpiration_profile
        transpiration_rate[0:10] = [transpiration[ind]]*10
        model.tree.transpiration_rate = np.asarray(transpiration_rate).reshape(40, 1)
        # print('New transpiration rate')
        # print(model.tree.transpiration_rate.reshape(40,))
        # print('Running from ' + str(round(time, 2)) + 'h to ' + str(round(daytime[ind+1], 2)) + 'h')
        model.run(time_start=time*60*60, time_end=daytime[ind+1]*60*60, dt=0.01, output_interval=60)

    # set transpiration to zero
    model.tree.transpiration_rate = np.asarray(transpiration_profile).reshape(40, 1)
    model.run(time_start=12.00001*60*60, time_end=24*60*60, dt=0.01, output_interval=300)

    print('Model simulation finished')
