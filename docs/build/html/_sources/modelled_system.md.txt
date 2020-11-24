Description of the modelled system
==================================

Introduction
------------

![Modelled system](../../source/_static/modelled_system.png "modelled system")


Modelled system is a 2-dimensional array representation of a tree. Each property of the tree is saved in a NumPy array. 
The first column in each array represents values of the xylem and the second column values of the phloem. 
The rows of the arrays represent vertical elements in the tree. The height of an element (l<sub>i</sub>) is calculated 
from the input total tree height (H) and number of elements (N)

$$l_i = \\frac{H}{N}$$

The height of an element does not change during the simulation.

The row numbering of the arrays starts from the top of the tree (i=0) and ends in the bottom of the tree (i=N).

The number of elements that transpire, produce or unload sugar can be changed by the user. Currently only the bottommost
element takes up water from the soil. The model is build with the idea, that the sap in both xylem and phloem can contain
any number of soluts. However, currently functionality to set, retrieve and update solutes is limited only to sugar (sucrose)
in the phloem.

Input and output of the model
-----------------------------
![Model input and output](../../source/_static/model_input_output.png "model input and output")

See [instructions how to run the model](instructions_to_run.md) how to give the input to the model.

The fluxes inside the tree
-------------------------------------------------------------------------------

![Fluxes](../../source/_static/fluxes.png "calculated fluxes in the model")

The fluxes shown in the figure above are calculated as follows (the index i refers to rows in any 2D property array in the model). The fluxes are in the units (kg/s)

$$ Q_{ax,i} = Q_{ax,bottom,i} + Q_{ax,top,i} - E_i$$

$$ Q_{ax,bottom,i} = \\frac{k_i \\: A_{ax,i} \\: \\rho_w}{\eta_i \\: l_i}(P_{i+1} - P_{i} - P_h) $$

$$ Q_{ax,top,i} = \\frac{k_i \\: A_{ax,i+1} \\: \\rho_w}{\eta_i \\: l_i}(P_{i-1} - P_{i} + P_h) $$

$$ Q_{radial,phloem} = L_r A_{rad,i} \\rho_{w} [P_{i,xylem} - P_{i,phloem} - \\sigma(C_{i,xylem} - C_{i,phloem}) R T)] $$

$$ Q_{radial,xylem} = -Q_{radial,phloem} $$

where

* E<sub>i</sub>: transpiration rate of the ith element (kg/s)
* k<sub>i</sub>: axial permeability of the ith element (m<sup>2</sup>)
* A<sub>ax,i</sub>: base surface area of xylem or phloem (m<sup>2</sup>)
* &rho;<sub>w</sub>: liquid phase density of water (kg/m<sup>3</sup>)
* &eta;: viscosity of the sap in the ith element (Pa s)
* l<sub>i</sub>: length (height) of the ith element (m)
* P<sub>i</sub>: Pressure in the ith element (Pa)
* P<sub>h</sub>: Hydrostatic pressure (Pa) - P<sub>h</sub> = &rho;<sub>w</sub> a<sub>gravitation</sub> l<sub>i</sub>
* L<sub>r</sub>: radial hydraulic conductivity (m/Pa/s)
* A<sub>rad,i</sub>: lateral surface area of the xylem (m<sup>2</sup>)
* &sigma;: reflection coefficient (Van't hoff factor) (unitless)
* C<sub>i</sub>: sugar concentration in the ith element (mol/m<sup>3</sup>)
* R: universal gas constant (J/K/mol)
* T: ambient temperature (K)

The change of pressure, sugar concentration and element radius due to sap flux
------------------------------------------------------------------------------
The equations are calcualted based on [Hölttä et al., (2006)](https://link.springer.com/article/10.1007/s00468-005-0014-6)

The change of pressure in the xylem and the phloem due to sap flux is calculated separately for both horizontal elements according to

$$ \frac{\text{d}P_i}{\text{d}t} = \frac{\varepsilon_i}{V_i \rho_w}(Q_{ax,i}+Q_{rad,i}) $$

The change in sugar concentration in the phloem is calculated according to

$$ \frac{\text{d}C_i}{\text{d}t} = \frac{Q_{ax,i}}{V_i} \left(\frac{C_i}{\rho_w} + L_i + U_i \right) $$

The change in element radius in the xylem is calculated according to

$$ \frac{\text{d}R_{xylem,i}}{\text{d}t} = \frac{Q_{ax,i} + Q_{rad,i}}{2 \pi \rho_w l_i (R_{xylem,i}+R_{heartwood})} $$

The change in element radius in the phloem is calculated according to

$$ \frac{\text{d}R_{i,phloem}}{\text{d}t} = \frac{Q_{ax,i} + Q_{rad,i}}{2 \pi \rho_w l_i} \left( \frac{1}{\Sigma_{j=1}^3 R_{j,i}} - \frac{R_{phloem,i}}{(\Sigma_{j=1}^2 R_{j,i})(\Sigma_{j=1}^3 R_{j,i}} \right) $$

where the sum goes from j=1 (the heartwood) to j=3 (the phloem)

* &varepsilon;<sub>i</sub>: elastic modulus of the ith element (Pa)
* V<sub>i</sub>: volume of the ith element (m<sup>3</sup>)
* L<sub>i</sub>: sugar loading rate of the ith element (mol/s)
* U<sub>i</sub>: sugar unloading rate of the ith element (mol/s)
* &rho;<sub>w</sub>: liquid phase density of water (kg/m<sup>3</sup>)
* l<sub>i</sub>: length (height) of the ith element (m)
* C<sub>i</sub>: Sucrose concentration in the ith element (mol/m<sup>3</sup>)
* Q<sub>ax/rad,i</sub>: axial or radial sap flux of the ith element (kg/s)
* R<sub>i</sub>: radius of the ith element (m)
* P<sub>i</sub>: pressure in the ith element (Pa)

Sugar loading and unloading rates
---------------------------------

Currently in the model the sugar loading rate should be set equal to the photosynthesis rate when the tree is created.

The initial unloading rate needs to be set but it is updated according to
[Nikinmaa et. al., (2014)](https://academic.oup.com/aob/article/114/4/653/2769025)

$$ U_i = A_{rad,i} \max{ [0, u_s (C_i - C_0)]}$$

where

* A<sub>rad,i</sub>: radial surface area between the xylem and the phloem (i.e., lateral area of the xylem)
* u<sub>s</sub>: sugar unloading slope set
* C<sub>i</sub>: sugar concentration in the ith element
* C<sub>0</sub>: sugar target concentration


How the element radii are modelled
-------------------------------------------------
![Tree dimension](../../source/_static/tree.png "definition of tree radii")

The radius of each horizontal element (heartwood, xylem and phloem) are given as the width of the element.
For example, **given total tree diameter of 10 cm, heartwood radius of 2 cm and 90/10% split between the phloem
and the xylem**. The traditional radii of each element (distance from the pith) would be
* Heartwood: 2 cm
* Xylem: 4.7 cm
* Phloem 5 cm

**In the model the radiis are given as the width of the horizontal element, i.e.:**
* Heartwood: 2 cm
* Xylem: 2.7 cm
* Phloem: 0.3 cm