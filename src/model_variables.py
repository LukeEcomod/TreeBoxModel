all_variables = {'height': ['height of the tree', 'm', ("index"), 'f4'],
                 'element_height': ['Height of each element', 'm', ("index", "axial_layers"), 'f4'],
                 'radius': ['radius of an element', 'm', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'num_elements':
                 ['number of tree+root elements', 'pcs.', ("index"), 'i4'],
                 'transpiration_rate':
                 ['Transpiration rate in a tree element', 'kg/s', ("index", "axial_layers"), 'f4'],
                 'photosynthesis_rate':
                 ['Photosynthesis rate in a tree element', 'mol/s', ("index", "axial_layers"), 'f4'],
                 'sugar_concentration':
                 ['Sugar concentration in a tree element', 'mol/m3', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'sugar_loading_rate':
                 ['Sugar loading rate in a tree element', 'mol/s', ("index", "axial_layers"), 'f4'],
                 'sugar_unloading_rate':
                 ['Sugar unloading rate in a tree element', 'mol/s', ("index", "axial_layers"), 'f4'],
                 'axial_permeability':
                 ['axial_permeability (k) in a tree element', 'm^2', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'viscosity':
                 ['viscosity of solution in a tree element', 'Pa s', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'radial_hydraulic_conductivity':
                 ['radial hydraulic conductivity (l) in a tree element', 'm/Pa/s',
                  ("index", "axial_layers", "radial_layers"), 'f4'],
                 'elastic_modulus':
                 ['Elastic modulus of a tree element (E)', 'Pa', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'dqrad':
                 ['radial solution flux', 'kg/s', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'dqax':
                 ['axial solution flux', 'kg/s', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'dqax_up':
                 ['upward axial solution flux', 'kg/s', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'dqax_down':
                 ['downward axial solution flux', 'kg/s', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'pressure':
                 ['pressure of an element', 'Pa', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'simulation_time':
                 ['Time in simulation', 's', ("index"), 'f4'],
                 'model_index':
                 ['number of index in the model loop', '#', ("index"), 'i4'],
                 'area':
                 ['Element base area', 'm^2', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'volume_water':
                 ['Element volume', 'm^3', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'soil_dz':
                 ['Soil layer thickness', 'm', ("index", "soil_elements"), 'f4'],
                 'soil_pressure':
                 ['Soil water potential / pressure', 'Pa', ("index", "soil_elements"), 'f4'],
                 'soil_kh':
                 ['Soil hydraulic conductivity', 'm/s', ("index", "soil_elements"), 'f4'],
                 'soil_root_k':
                 ['Soil-root system total conductivity', '1/s', ("index", "root_elements"), 'f4'],
                 'rooting_depth':
                 ['Maximum depth of root system', 'm', ("index"), 'f4'],
                 'root_area_density':
                 ['Area density of fine roots in element', 'm^2/m^3', ("index", "root_elements"), 'f4'],
                 'area_per_tree':
                 ['Ground area of the tree', 'm^2', ("index"), 'f4'],
                 'dqroot':
                 ['Root water uptake', 'kg/s', ("index", "axial_layers"), 'f4'],
                 'sapflow':
                 ['Sapflow velocity in the tree', 'm/s', ("index", "axial_layers", "radial_layers"), 'f4']}

cmp_variables = {
    'cmp_element_height':
    ['height of the tree elements', 'm', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_concentration':
    ['Concentration of the gas', 'mol/m^3', ("index", "axial_layers", "radial_layers", "space_layers"), 'f4'],
    'cmp_space_division':
    ['Volume fraction of air / water + cell space in an element', '1',
     ("index", "axial_layers", "radial_layers", "space_layers"), 'f4'],
    'cmp_element_radii':
    ['Radii of each element', 'm', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_head_area':
    ['Head (top) area of each element', 'm^2', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_element_volume':
    ['Volume of each element', 'm^3', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_element_volume_air':
    ['Air volume of each element', 'm^3', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_element_volume_water':
    ['Water volume of each element', 'm^3', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_element_volume_cell':
    ['Cell volume of each element', 'm^3', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_gas_radial_diffusion_coefficient':
    ['Diffusion coefficient of CH4 in each element', 'm^2/s', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_gas_axial_diffusion_coefficient':
    ['Diffusion coefficient of CH4 in each element', 'm^2/s', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_eq_rate':
    ['Equilibration rate of CH4 between water and air', '1/s', ("index"), 'f4'],
    'cmp_velocity':
    ['Sap flow velocity in the tree', 'm/s', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_henry_coef':
    ['Henrys law coefficient Kh = c_water/c_air', '1', ("index"), 'f4'],
    'cmp_temperature':
    ['Temperature of the element', 'K', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_ambient_concentration':
    ['Ambient concentration of the gas', 'mol/m^3', ("index", "axial_layers"), 'f4'],
    'cmp_moles_out':
    ['Amount of gas that has diffused out of the tree', 'mol', ("index", "axial_layers"), 'f4'],
    'cmp_simulation_time':
    ['Simulation time in the gas simulation', 's', ("index"), 'f4'],
    'cmp_gas_radial_diffusion_flux':
    ['Radial Diffusion flux of the compound in the gas phase', 'mol/s', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_gas_axial_diffusion_flux':
    ['Axial diffusion flux of the compound in the gas phase', 'mol/s', ("index", "axial_layers", "radial_layers"), 'f4'],
    'cmp_eq_flux':
    ['Equilibrium flux between aqueous and gas phase of an element', 'mol/s',
     ("index", "axial_layers", "radial_layers", "space_layers"), 'f4'],
    'cmp_aq_axial_advection_flux':
    ['Aqueous phase axial advection flux of the compound', 'mol/s',
     ("index", "axial_layers", "radial_layers",), 'f4'],
    'cmp_moles_in':
    ['Molar amount of the compound inside the tree', 'mol',
     ("index", "axial_layers", "radial_layers", "space_layers"), 'f4']
}

index_dim_vars = [key for key, val in all_variables.items() if val[2] == ("index")]
soil_dim_vars = [key for key, val in all_variables.items() if val[2] == ("index", "soil_elements")]
root_dim_vars = [key for key, val in all_variables.items() if val[2] == ("index", "root_elements")]
axial_layer_dim_vars = [key for key, val in all_variables.items() if val[2] == ("index", "axial_layers")]
gas_three_dim_vars = [key for key, val in cmp_variables.items() if val[2] == ("index", "axial_layers",
                                                                              "radial_layers", "space_layers")]
gas_index_dim_vars = [key for key, val in cmp_variables.items() if val[2] == ("index")]
gas_axial_dim_vars = [key for key, val in cmp_variables.items() if val[2] == ("index", "axial_layers")]
