M_WATER: float = 0.0182  # molar mass of water unit: kg/mol
RHO_WATER: float = 1000  # density of liquid water unit: kg/m3
VISCOSITY_WATER: float = 1e-3  # dynamic viscosity of water, unit: Pa s

M_SUCROSE: float = 0.3423  # molar mass of sucrose unit: kg/mol
RHO_SUCROSE: float = 1590.0  # density of sucrose unit: kg/m3

GRAVITATIONAL_ACCELERATION: float = 9.81  # acceleration due to Earth's gravity unit: m/s2
AVOGADROS_CONSTANT: float = 6.022e23  # avogadro's constant unit: 1/mol

# These are needed for water volume calculation
# TODO: add these to treeElement parameters
# From Hölttä et al., (2016)
HEARTWOOD_RADIUS: float = 60e-3
XYLEM_RADIUS: float = 80e-3
PHLOEM_RADIUS: float = 8e-3

MAX_ELEMENT_COLUMNS: int = 2  # Max number of columns in tree.elements
