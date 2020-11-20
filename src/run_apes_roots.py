import xarray as xr
from datetime import datetime
from src.model import Model
from src.tree import Tree
from src.soil import Soil
from src.roots import Roots
from src.constants import GRAVITATIONAL_ACCELERATION, RHO_WATER
import numpy as np
import math

# read pyAPES data
data = xr.open_dataset('pyapes_results/202010151003_pyAPES_results.nc')
# choose tree i=0 is tree1 ... etc
i = 1
# create the soil object
soil_layer_thickness = np.array(np.diff(-data.soil_z))
pressure = np.zeros((len(soil_layer_thickness), 1))
hydraulic_conductivity = np.ones((len(soil_layer_thickness), 1))
soil = Soil(layer_thickness=soil_layer_thickness,
            hydraulic_conductivity=hydraulic_conductivity,
            pressure=pressure)
# create the roots object
# set area density until 1 meter so that RAI=30 using a linear function
root_elements = 24
area_density = np.zeros((root_elements, 1))
LAI_tree = np.asarray([8., 28., 49.])
RAI = 2*LAI_tree[i]
area_density = np.ones((root_elements, 1))
effective_radius = 0.5e-3*np.ones((root_elements, 1))
rooting_depth = 0.5
roots = Roots(rooting_depth=rooting_depth, area_density=area_density,
              effective_radius=effective_radius, soil_conductance_scale=3e8, num_elements=root_elements)

dz = roots.layer_thickness(soil)
length = roots.layer_depth(soil)
a = 200
b = RAI/(np.sum(np.exp(-a*length)*dz))

roots.area_density = b*np.exp(-a*length).reshape(roots.num_elements, 1)
# roots.area_density = -a**(-length)*np.log(a)

daytime = np.linspace(0, 24*9, 433)

LAI = [0.2, 1.4, 0.5]
tree_height = [13., 18., 20.]
diameter = [10., 18., 25.]
radii_all = [[0.02130682, 0.02840909, 0.00028409], [0.03835227, 0.05113636, 0.00051136],
             [0.05326705, 0.07102273, 0.00071023]]


height = tree_height[i]
tree_elements = int(height/0.25)
num_elements = tree_elements + roots.num_elements
element_height = np.concatenate((np.repeat(height/tree_elements, repeats=tree_elements),
                                 roots.layer_thickness(soil).reshape(roots.num_elements,)))
print(num_elements)
radii = radii_all[i]
transpiration_all = np.asarray(data.pt_transpiration[:, 0, i+1, 1: tree_elements+1]*0.25/LAI[i]*LAI_tree[i]*0.0180152)

photosynthesis_all = np.asarray(data.pt_net_co2[:, 0, i+1, 1: tree_elements+1]*0.25/LAI[i]*LAI_tree[i]*1e-6/12)

# soil_water_potential_all = np.asarray(data.soil_water_potential[:, 0, 10]*1e4)

sugar_profile: np.ndarray = np.zeros((num_elements, 1))
upper_third = int(num_elements/3)
middle_third = int(num_elements/3*2)
sugar_profile[0: upper_third, 0] = 1400
sugar_profile[upper_third: middle_third, 0] = 1000
sugar_profile[middle_third: num_elements, 0] = 800

sugar_loading_profile = [0.0 for i in range(num_elements)]
transpiration_profile = [0.0 for i in range(num_elements)]
photosynth_profile = [0.0 for i in range(num_elements)]

sugar_unloading_profile = [0.0 for i in range(num_elements)]

sugar_target_concentration: float = 1200

sugar_unloading_slope = 3.5e-7

axial_permeability_profile: np.ndarray = np.zeros((num_elements, 2))
axial_permeability_profile[:, 0] += 1.5e-12
axial_permeability_profile[:, 1] += 6.0e-12
tree_ind = np.arange(0, tree_elements)
root_ind = np.arange(tree_elements, num_elements)
for (row, permeability) in enumerate(axial_permeability_profile[tree_ind]):
    axial_permeability_profile[row, 0] = np.min([permeability[0],
                                                 permeability[0]*1/math.sqrt(height/tree_elements*(row+1))])
    axial_permeability_profile[row, 1] = np.min([permeability[1],
                                                 permeability[1]*1/math.sqrt(height/tree_elements*(row+1))])


axial_permeability_profile[tree_ind] = np.flip(axial_permeability_profile[tree_ind], axis=0)

axial_permeability_profile[root_ind] = [1.5e-14, 6.0e-14]

radial_hydr_conductivity = [1e-13] * num_elements

elastic_modulus_profile = [[1000e6, 30e6]] * num_elements

tree = Tree(height=height,
            num_elements=num_elements,
            element_height=element_height,
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
            roots=roots)

outputfname = 'pyapes_'
# outputfname = outputfname + 'tree' + str(i+1) + '_' + datetime.now().strftime("%y-%m-%dT%H:%M:%S") + ".nc"
outputfname = outputfname + 'tree' + str(i+1) + '_' + datetime.now().strftime("%y-%m-%dT%H:%M:%S") + 'roots.nc'

model = Model(tree, soil, outputfile=outputfname)
model.tree.pressure[model.tree.root_elements, 0] = np.asarray(data.soil_water_potential[0, 0, :24]
                                                              * RHO_WATER*GRAVITATIONAL_ACCELERATION)
base_permeability = np.array(axial_permeability_profile).reshape(num_elements, 2)
max_pressure = 2e6
scale_exponent = 2
for (ind, time) in enumerate(daytime[0:-1]):
    # set new transpiration, photosynthesis rate and soil water potential
    print(model.tree.pressure[model.tree.tree_elements]/1e6)
    model.tree.transpiration_rate[model.tree.tree_elements, 0] = np.flip(transpiration_all[ind, :])\
        .reshape(tree_elements,)

    model.soil.pressure = (np.asarray(data.soil_water_potential[ind, 0, :-1])*RHO_WATER*GRAVITATIONAL_ACCELERATION)\
         .reshape(model.soil.num_elements, 1)
    model.soil.hydraulic_conductivity = np.asarray(data.soil_hydraulic_conductivity[ind, 0, :-1])\
         .reshape(model.soil.num_elements, 1)
    # model.soil.hydraulic_conductivity = model.soil.hydraulic_conductivity * 1e9
    model.tree.photosynthesis_rate[model.tree.tree_elements, 0] = np.flip(photosynthesis_all[ind, :])\
        .reshape(tree_elements,)
    model.tree.sugar_loading_rate[model.tree.tree_elements] = model.tree.photosynthesis_rate[model.tree.tree_elements]\
        .copy()
    model.run_scipy(time_start=daytime[ind]*60*60, time_end=daytime[ind+1]*60*60, ind=ind)
