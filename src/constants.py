import math


M_WATER: float = 0.0182  #: Molar mass of water :math:`\left(\frac{kg}{mol}\right)`

RHO_WATER: float = 1000  #: density of liquid water :math:`\left(\frac{kg}{m^3}\right)`

VISCOSITY_WATER: float = 1e-3  #: dynamic viscosity of water :math:`\left( Pa \cdot s \right)`

M_SUCROSE: float = 0.3423  #: molar mass of sucrose :math:`\left(\frac{kg}{mol}\right)`

RHO_SUCROSE: float = 1590.0  #: density of sucrose :math:`\left(\frac{kg}{m^3}\right)`

GRAVITATIONAL_ACCELERATION: float = 9.81  #: acceleration due to Earth's gravity :math:`\left(\frac{m}{s^2}\right)`

AVOGADROS_CONSTANT: float = 6.022e23  #: avogadro's constant :math:`\left(\frac{1}{mol}\right)`

MOLAR_GAS_CONSTANT: float = 8.3145  #: molar gas constant :math:`\left(\frac{J}{K \cdot mol}\right)`

# TODO: add these to tree parameters
# From Hölttä et al., (2016)
# HEARTWOOD_RADIUS: float = 0.55e-2  # for nikinmaa
HEARTWOOD_RADIUS = 4.44e-2  #: The radius of the heartwood :math:`\left( m \right)`
XYLEM_RADIUS: float = 0.52e-2  # for nikinmaa
PHLOEM_RADIUS: float = 1e-3

MAX_ELEMENT_COLUMNS: int = 2  #: Max number of columns in the [tree](index.html#src.tree.Tree) class

TEMPERATURE: float = 298.0  #: Temperature of the tree :math:`\left( K \right)`
