class TreeElement:
    def __init__(self, pressure, solutes, viscosity, permeability, height,
                 hydraulic_conductivity, transpiration_rate,
                 photosynthesis_rate, sugar_loading_rate,
                 sugar_unloading_rate):
        self.pressure = pressure  # unit: Pa
        self.solutes = solutes  # Solute array
        self.viscosity = viscosity  # unit: Pa s
        self.permeability = permeability  # unit: m2
        self.height = height  # unit: m
        self.hydraulic_conductivity = hydraulic_conductivity  # unit: m/(Pa s)
        self.transpiration_rate = transpiration_rate  # unit: m3/s
        self.photosynthesis_rate = photosynthesis_rate  # unit: mol/s
        self.sugar_loading_rate = sugar_loading_rate  # unit: mol/s
        self.sugar_unloading_rate = sugar_unloading_rate  # unit: mol/s
