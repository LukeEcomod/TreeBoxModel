.. Tree Box Model documentation master file, created by
   sphinx-quickstart on Tue May 19 13:49:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Tree Box Model documentation
============================

.. toctree::
   :maxdepth: 2
   :caption: Details of the model

   modelled_system
   instructions_to_run

Installation
------------

- Download the source

>>> git clone git@github.com:LukeEcomod/TreeBoxModel.git

or

download the source https://github.com/LukeEcomod/TreeBoxModel

- Install the required packages

Ideally you have created a new virtual environment for this project.

To install all the packages required for the model to run use 

>>> pip install -r requirements_pip.txt

Quick start
-----------
run main.py.

>>> python main.py

See instructions to run the model for detailed instructions.

main.py
-------
.. automodule:: main

Modules, Classes \& functions
------------------------------

.. automodule:: src.model
   :members:
   :undoc-members:

.. automodule:: src.tree
   :members:
   :undoc-members:

.. automodule:: src.solute
   :members:
   :undoc-members:

.. automodule:: src.soil
   :members:
   :undoc-members:

.. automodule:: src.roots
   :members:
   :undoc-members:

.. autofunction:: src.odefun.odefun

.. autodata:: src.constants.MAX_ELEMENT_COLUMNS

.. autodata:: src.constants.TEMPERATURE

.. autodata:: src.constants.M_WATER

.. autodata:: src.constants.RHO_WATER

.. autodata:: src.constants.VISCOSITY_WATER

.. autodata:: src.constants.M_SUCROSE

.. autodata:: src.constants.RHO_SUCROSE

.. autodata:: src.constants.GRAVITATIONAL_ACCELERATION

.. autodata:: src.constants.AVOGADROS_CONSTANT

.. autodata:: src.constants.MOLAR_GAS_CONSTANT

.. autodata:: src.model_variables
   :annotation: = Name, descriptions, unit, dimension and precision of each variable that is saved to the netcdf dataset

.. automodule:: src.tools.iotools
   :members:
   :undoc-members:

.. automodule:: src.tools.plotting
   :members:
   :undoc-members:

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

