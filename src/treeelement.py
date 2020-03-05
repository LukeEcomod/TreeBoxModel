from typing import List
from .solute import Solute
from dataclasses import dataclass


@dataclass
class TreeElement:
    pressure: float  # unit: Pa

    solutes: List[Solute]  # Solute array

    viscosity: float  # unit: Pa s

    permeability: float  # unit: m2

    elastic_modulus: float  # unit: Pa

    hydraulic_conductivity: float  # unit: m/(Pa s)

    height: float  # unit: m

    transpiration_rate: float  # unit: m3/s

    photosynthesis_rate: float  # unit: mol/s

    sugar_loading_rate: float  # unit: mol/s

    sugar_unloading_rate: float  # unit: mol/s
