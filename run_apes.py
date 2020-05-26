import xarray as xr
from src.model import Model
from src.tree import Tree
from typing import List
import numpy as np
import math
from datetime import datetime

# read pyAPES data
data = xr.open_dataset('pyapes_res.nc')

height: float = 15.0

num_elements: int = 60

radii: List[float] = [2.96e-2, 1e-3]
daytime = np.linspace(0, 24*9, 433)

transpiration_all = np.asarray(data.pt_transpiration[:, 0, 1, 1:61]*0.25/2.1*20*0.0180152)

photosynthesis_all = np.asarray(data.pt_net_co2[:, 0, 1, 1:61]*0.25/2.1*20*1e-6/12)

sugar_profile: np.ndarray = np.zeros((num_elements, 1))
sugar_profile[0:30, 0] = 1400
sugar_profile[30:40, 0] = 1000
sugar_profile[40:60, 0] = 800

sugar_loading_profile = [0.0 for i in range(num_elements)]
transpiration_profile = [0.0 for i in range(num_elements)]
photosynth_profile = [0.0 for i in range(num_elements)]

sugar_unloading_profile: List[float] = [0.0 for i in range(num_elements)]

sugar_target_concentration: float = 1200

sugar_unloading_slope = 3.5e-7

axial_permeability_profile: np.ndarray = np.zeros((num_elements, 2))
axial_permeability_profile[:, 0] += 1.5e-12
axial_permeability_profile[:, 1] += 6.0e-12

for (row, permeability) in enumerate(axial_permeability_profile):
    axial_permeability_profile[row, 0] = np.min([permeability[0],
                                                 permeability[0]*1/math.sqrt(height/num_elements*(row+1))])
    axial_permeability_profile[row, 1] = np.min([permeability[1],
                                                 permeability[1]*1/math.sqrt(height/num_elements*(row+1))])

axial_permeability_profile = np.flip(axial_permeability_profile, axis=0)

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
            sugar_target_concentration=sugar_target_concentration,
            sugar_unloading_slope=sugar_unloading_slope,
            axial_permeability_profile=axial_permeability_profile,
            radial_hydraulic_conductivity_profile=radial_hydr_conductivity,
            elastic_modulus_profile=elastic_modulus_profile,
            ground_water_potential=ground_water_potential)

outputfname = 'pyapes_'
outputfname = outputfname + datetime.now().strftime("%y-%m-%dT%H:%M:%S") + ".nc"
model = Model(tree, outputfile=outputfname)

for (ind, time) in enumerate(daytime[0:-1]):
    # set new transpiration rate

    model.tree.transpiration_rate = np.flip(transpiration_all[ind, :]).reshape(num_elements, 1)

    model.tree.photosynthesis_rate = np.flip(photosynthesis_all[ind, :]).reshape(num_elements, 1)
    model.tree.sugar_loading_rate = model.tree.photosynthesis_rate.copy()
    model.run_scipy(time_start=daytime[ind]*60*60, time_end=daytime[ind+1]*60*60, ind=ind)
