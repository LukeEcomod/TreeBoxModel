from dataclasses import dataclass


@dataclass
class Solute:

    molar_mass: float  # unit: kg/mol
    density: float  # unit: kg/m3
    concentration: float  # unit: mol/m3
