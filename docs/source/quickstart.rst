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
2. Build a force-field model for the perovskite.
3. Specify which sites are **B-sites** and the allowed B-site species.
4. Optimize the structure by calling ``optimize.optimize()``, which conducts Monte Carlo swaps and records accepted structures / energies.
5. (Optional) Restart the optimization from a previous run by calling ``optimize.restart()``.

Example
---------------------

This example demonstrates optimizing the structure of :math:`\mathrm{Pb\left(Mg_{1/3}Nb_{2/3}\right)O_3}` (PMN), with a :math:`12\times 12\times 12` simulation cell.
The **B-sites** are occupied by Mg/Nb atoms. (Replace the configuration and model with your own files.)

.. code-block:: python

   from perostruc import optimize

   optimize.optimize(filename = './PMN_conf.lmp',
                     model = "./PMN_model.pb",
                     temperature = 300.0, 
                     cell_size = np.array([12, 12, 12]), 
                     element_list = ['Mg', 'Nb','O','Pb'], 
                     B_site_elements = ['Mg', 'Nb'], 
                     neighbor_cutoff = 6,    # allow nearest and next-nearest neighbor swap.
                     traj_dir = 'trajs', 
                     logging_dir = 'logs', 
                     nepoch = 50, 
                     fmax = 0.03,            # force threshold for the FIRE optimization
                     fire_max_steps = 50,    # maximum number of steps for the FIRE optimization
                     )
    
    optimize.restart(model = "./PMN_model.pb",
                     temperature = 300.0, 
                     cell_size = np.array([12, 12, 12]), 
                     element_list = ['Mg', 'Nb','O','Pb'], 
                     B_site_elements = ['Mg', 'Nb'], 
                     neighbor_cutoff = 6,    # allow nearest and next-nearest neighbor swap.
                     traj_dir = 'trajs', 
                     logging_dir = 'logs', 
                     nepoch = 100, 
                     fmax = 0.03,          # force threshold for the FIRE optimization
                     fire_max_steps = 50,    # maximum number of steps for the FIRE optimization
                     )


Tips
----

- Start at higher temperature to explore configurations, then cool down (annealing).
- Adjust the arguments ``nepoch``, ``fmax``, ``fire_max_steps`` if the convergence is bad.
