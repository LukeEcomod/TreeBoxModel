from dataclasses import dataclass


@dataclass
class Solute:
    """ Contains the variables to model a solute compound

        Args:
            name (str): Name of the compound
            molar_mass (float): Molar mass of the compound (:math:`\\frac{kg}{mol}`)
            density: Liquid phase density of the compound (:math:`\\frac{kg}{m^3}`)
            concentration: Concentration of the compound in the sap solution (:math:`\\frac{mol}{m^3}`)
    """
    name: str
    molar_mass: float
    density: float
    concentration: float
