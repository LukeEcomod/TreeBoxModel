"""

"""

import math

M_WATER: float = 0.0182  # molar mass of water unit: kg/mol
RHO_WATER: float = 1000  # density of liquid water unit: kg/m3
VISCOSITY_WATER: float = 1e-3  # dynamic viscosity of water, unit: Pa s

M_SUCROSE: float = 0.3423  # molar mass of sucrose unit: kg/mol
RHO_SUCROSE: float = 1590.0  # density of sucrose unit: kg/m3

GRAVITATIONAL_ACCELERATION: float = 9.81  # acceleration due to Earth's gravity unit: m/s2
AVOGADROS_CONSTANT: float = 6.022e23  # avogadro's constant unit: 1/mol
MOLAR_GAS_CONSTANT: float = 8.3145  # molar gas constant unit: J/K/mol

# These are needed for water volume calculation
# TODO: add these to treeElement parameters
# From Hölttä et al., (2016)
HEARTWOOD_RADIUS: float = 0.55e-2
XYLEM_RADIUS: float = 0.52e-2
PHLOEM_RADIUS: float = 1e-3

MAX_ELEMENT_COLUMNS: int = 2  # Max number of columns in tree.elements

TEMPERATURE: float = 298.0  # Temperature of the tree unit: K
