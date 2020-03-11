from src.model import Model
from src.tree import Tree
from typing import List


if __name__ == "__main__":

    height: float = 12.0

    num_elements: int = 40

    radii: List[float] = [1.0, 2.0]
    transpiration_profile: List[float] = [0 for i in range(num_elements)]
    transpiration_profile[0] = 0.9*1e-6  # m3/s

    photosynth_profile: List[float] = [0 for i in range(num_elements)]
    photosynth_profile[0:3] = [1e-6, 1e-6, 1e-6]

    sugar_profile: List[float] = [10.0 for i in range(num_elements)]
    sugar_loading_profile = photosynth_profile

    sugar_unloading_profile: List[float] = [0.0 for i in range(num_elements)]
    sugar_unloading_profile[-1] = 2

    axial_permeability_profile: List[List[float]] = [[1.5e-12, 6.0e-12]]*num_elements

    radial_hydr_conductivity: List[float] = [1e-13] * num_elements

    elastic_modulus_profile: List[List[float]] = [[1000e6, 30e6]] * num_elements

    ground_water_potential: float = 0.0

    tree = Tree(height=height,
                num_elements=num_elements,
                initial_radius=radii,
                transpiration_profile=transpiration_profile,
                photosynthesis_profile=photosynth_profile,
                sugar_profile=sugar_profile,
                sugar_loading_profile=sugar_loading_profile,
                sugar_unloading_profile=sugar_unloading_profile,
                axial_permeability_profile=axial_permeability_profile,
                radial_hydraulic_conductivity_profile=radial_hydr_conductivity,
                elastic_modulus_profile=elastic_modulus_profile,
                ground_water_potential=ground_water_potential)

    model = Model(tree, outputfile="test.out")

    model.run(time_start=1e-5, time_end=5.0, dt=1.0)
    print('Model simualtion finished')
