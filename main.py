from src.constants import PHLOEM_RADIUS, XYLEM_RADIUS
from src.model import Model
from src.tree import Tree
from typing import List
import numpy as np
import math
from datetime import datetime

if __name__ == "__main__":

    height: float = 12.0

    num_elements: int = 40

    radii: List[float] = [XYLEM_RADIUS, PHLOEM_RADIUS]
    transpiration_profile: List[float] = [0.0 for i in range(num_elements)]
    # transpiration_profile[0:10] = [0.9*1e-3/10.0]*10  # kg/s
    transpiration_max = 0.9e-3/10.0  # kg/s
    daytime = np.linspace(0, 12, 12*12+1)  # hours starting from 06:00 AM, 5 minute interval
    transpiration = np.sin(daytime*math.pi/12.0)*transpiration_max
    photosynth_profile: List[float] = [0 for i in range(num_elements)]
    # photosynth_profile[0:5] = [2e-5, 2e-5, 2e-5, 2e-5, 2e-5]
    photosynthesis_max = 2.5e-5/20.0
    photosynthesis = np.sin(daytime*math.pi/12.0)*photosynthesis_max

    sugar_profile: np.ndarray = np.zeros((num_elements, 1))
    sugar_profile[:, 0] += 1500

    sugar_loading_profile = photosynth_profile

    sugar_unloading_profile: List[float] = [0.0 for i in range(num_elements)]

    sugar_target_concentration: float = 1400

    sugar_unloading_slope = 3.5e-7

    axial_permeability_profile: np.ndarray = np.zeros((num_elements, 2))
    axial_permeability_profile[:, 0] += 1.5e-12*0.5
    axial_permeability_profile[:, 1] += 6.0e-12*0.5

    for (row, permeability) in enumerate(axial_permeability_profile):
        axial_permeability_profile[row, 0] = np.min([permeability[0], permeability[0]*1/math.sqrt(0.3*(row+1))])
        axial_permeability_profile[row, 1] = np.min([permeability[1], permeability[1]*1/math.sqrt(0.3*(row+1))])

    axial_permeability_profile = np.flip(axial_permeability_profile, axis=0)

    radial_hydr_conductivity: List[float] = [1e-13] * num_elements
    # radial_hydr_conductivity = [0]*num_elements

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
                sugar_target_concentration=sugar_target_concentration,
                sugar_unloading_slope=sugar_unloading_slope,
                axial_permeability_profile=axial_permeability_profile,
                radial_hydraulic_conductivity_profile=radial_hydr_conductivity,
                elastic_modulus_profile=elastic_modulus_profile,
                ground_water_potential=ground_water_potential)

    outputfname = 'test_'
    outputfname = outputfname + datetime.now().strftime("%y-%m-%dT%H:%M:%S") + ".nc"
    model = Model(tree, outputfile=outputfname)
    for day in range(0, 4):
        for (ind, time) in enumerate(daytime[0:-1]):
            # set new transpiration rate
            transpiration_rate = transpiration_profile
            transpiration_rate[0:10] = [transpiration[ind]]*10
            model.tree.transpiration_rate = np.asarray(transpiration_rate).reshape(40, 1)

            photosynthesis_rate = photosynth_profile
            photosynthesis_rate[0:20] = [photosynthesis[ind]]*20
            model.tree.photosynthesis_rate = np.asarray(photosynthesis_rate).reshape(40, 1)
            # print('New transpiration rate')
            # print(model.tree.transpiration_rate.reshape(40,))
            # print('Running from ' + str(round(time, 2)) + 'h to ' + str(round(daytime[ind+1], 2)) + 'h')
            model.run_scipy(time_start=(day*60*60*24)+time*60*60, time_end=(day*60*60*24)+daytime[ind+1]*60*60, ind=ind)

        # set transpiration to zero
        model.tree.transpiration_rate = np.asarray([0.0 for i in range(num_elements)]).reshape(40, 1)
        model.tree.photosynthesis_rate = np.asarray([0.0 for i in range(num_elements)]).reshape(40, 1)
        model.tree.sugar_loading_rate = model.tree.photosynthesis_rate
        nighttime = np.linspace(12.0001, 24, 12*30+1)

        for (ind, time) in enumerate(nighttime[0:-1]):
            model.run_scipy(time_start=(day*60*60*24)+time*60*60,
                            time_end=(day*60*60*24)+nighttime[ind+1]*60*60, ind=ind)

    print('Model simulation finished')
