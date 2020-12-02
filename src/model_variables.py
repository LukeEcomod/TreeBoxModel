all_variables = {'height': ['height of the tree', 'm', ("index"), 'f4'],
                 'radius': ['radius of an element', 'm', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'num_elements': ['number of tree elements', 'pcs.', ("index"), 'i4'],
                 'transpiration_rate':
                 ['Transpiration rate in a tree element', 'kg/s', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'photosynthesis_rate':
                 ['Photosynthesis rate in a tree element', 'mol/s', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'sugar_concentration':
                 ['Sugar concentration in a tree element', 'mol/m3', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'sugar_loading_rate':
                 ['Sugar loading rate in a tree element', 'mol/s', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'sugar_unloading_rate':
                 ['Sugar unloading rate in a tree element', 'mol/s', ("index", "axial_layers", "radial_layers"), 'f4'],
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
                 'ground_water_potential':
                 ['Ground water potential', 'Pa', ("index"), 'f4'],
                 'pressure':
                 ['pressure of an element', 'Pa', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'num_elements':
                 ['number of tree elements', 'pcs.', ("index"), 'i4'],
                 'simulation_time':
                 ['Time in simulation', 's', ("index"), 'f4'],
                 'model_index':
                 ['number of index in the model loop', '#', ("index"), 'i4'],
                 'area':
                 ['Element base area', 'm^2', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'volume':
                 ['Element volume', 'm^3', ("index", "axial_layers", "radial_layers"), 'f4'],
                 'soil_dz':
                 ['Soil layer thickness', 'm', ("index", "soil_elements"), 'f4'],
                 'soil_pressure':
                 ['Soil water potential / pressure', 'Pa', ("index", "soil_elements"), 'f4'],
                 'soil_kh':
                 ['Soil hydraulic conductivity', 'm/s', ("index", "soil_elements"), 'f4'],
                 'soil_root_k':
                 ['Soil-root system total conductivity', '1/s', ("index", "root_elements"), 'f4'],
                 'dqroot':
                 ['Root water uptake', 'kg/s', ("index", "axial_layers"), 'f4']}

index_dim_vars = [key for key in all_variables.keys() if all_variables[key][2] == ("index")]
soil_dim_vars = [key for key in all_variables.keys() if all_variables[key][2] == ("index", "soil_elements")]
root_dim_vars = [key for key in all_variables.keys() if all_variables[key][2] == ("index", "root_elements")]
axial_layer_dim_vars = [key for key in all_variables.keys() if all_variables[key][2] == ("index", "axial_layers")]
