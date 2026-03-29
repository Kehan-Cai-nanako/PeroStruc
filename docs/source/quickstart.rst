Quickstart
==========

Install
-------

For quick install:

.. code-block:: bash

   python -m pip install .

For development (editable install):

.. code-block:: bash

   python -m pip install -e ".[docs]"

Minimal workflow
----------------

A typical workflow is:

1. Build an initial perovskite structure (e.g., ABO\ :sub:`3`).
2. Specify a force-field model for the perovskite. This step is not handled by this package. Please use your favorite force-field model. Note that this package can work with any force-field model as long as it can be used as a ASE calculator. For example, you can use the `Deep Potential (DP) <https://github.com/deepmodeling/deepmd-kit>`_ model, which already has a ASE calculator interface.
3. Use the ``FIRESwapOptimizer`` class to initialize the optimizer. You should specify the initial configuration, the force-field model, and the simulation parameters.
4. Optimize the structure by calling ``optimize()``, which conducts Monte Carlo swaps and records accepted structures / energies.
5. (Optional) Restart the optimization from a previous run by using the ``restart()`` method.

Example
-------

The following example demonstrates optimizing the structure of :math:`\mathrm{Pb\left(Mg_{1/3}Nb_{2/3}\right)O_3}` (PMN), with a :math:`6\times 6\times 6` simulation cell.

.. code-block:: python

   import numpy as np
   import ase
   import ase.io
   from perostruc import FIRESwapOptimizer, update_element
   from deepmd.calculator import DP    ## you need to install deepmd-kit to use the DP model. deepmd-kit is not included in this package.

   init_config = ase.io.read('initial_configs/disordered_L6X6X6.lmp', format='lammps-data', atom_style='atomic')

   ## This is the order of the elements in the lammps data file. We need this because lammps dump file does not have the element names. 
   ## Instead, it uses implcit numbering (1, 2, 3, ...). So we need to specify how the implicit numbering corresponds to the element names.
   element_specorder = ['Mg', 'Nb', 'O', 'Pb']   
   update_element(init_config, element_specorder)

   ## Load the force-field model. You can use any other force-field model as long as it has a ASE calculator interface. 
   dpmodel = DP(model="/global/cfs/projectdirs/m5025/Ferroic/PMN_paper_repo/01.DP_Model_Training/Production_Model/model-compress.pb")
   init_config.calc = dpmodel

   simulation = FIRESwapOptimizer(
      temperature = 600.0,  ## temperature in K for the Monte Carlo swaps.
      A_site_elements = ['Pb'],
      B_site_elements = ['Mg', 'Nb'],
      element_specorder = element_specorder,
      neighbor_cutoff = 6.0,  ## atoms within this cutoff will be considered for the swaps.
      fmax = 0.02, ## this can be made larger if the convergence is hard to achieve.
      fire_max_steps = 100, ## this is the maximum number of steps for each FIRE optimization.
      nepoch = 10000,
   )

   ## Initialize from scratch
   simulation.initialize(supercell_size = [6, 6, 6], init_config = init_config)
   ## Relax the initial geometric configuration using plain FIRE.
   simulation.run_FIRE(fmax = 0.02, steps = 1000)
   ## Conduct the FIRESwap optimization, which alternates between FIRE and Monte Carlo swaps.
   ## swap_type = 'B' means swapping B-site species. You can also set it to 'A' to swap A-site species.
   ## flexible_cell = False means fixing the cell parameters. You can set it to True to allow the cell parameters to change.
   simulation.optimize(swap_type = 'B', flexible_cell = False, profile = True)

Examples with more details can be found on `GitHub <https://github.com/Kehan-Cai-nanako/PeroStruc/tree/main/examples>`_.

Tips
----

- Start at higher temperature to explore configurations, then cool down (annealing).
- Adjust the arguments ``nepoch``, ``fmax``, ``fire_max_steps`` if the convergence is bad.
