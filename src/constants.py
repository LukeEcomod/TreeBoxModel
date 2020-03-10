M_WATER = 0.0182  # molar mass of water unit: kg/mol
RHO_WATER = 1000  # density of liquid water unit: kg/m3

M_SUCROSE = 0.3423  # molar mass of sucrose unit: kg/mol
RHO_SUCROSE = 1590  # density of sucrose unit: kg/m3

GRAVITATIONAL_ACCELERATION = 9.81  # acceleration due to Earth's gravity unit: m/s2
AVOGADROS_CONSTANT: float = 6.022e23  # avogadro's constant unit: 1/mol

# These are needed for water volume calculation
# TODO: add these to treeElement parameters
# From Hölttä et al., (2016)
HEARTWOOD_RADIUS = 60e-3
XYLEM_RADIUS = 80e-3
PHLOEM_RADIUS = 8e-3

MAX_ELEMENT_COLUMNS = 2  # Max number of columns in tree.elements
