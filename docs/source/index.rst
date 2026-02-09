PeroStruc: Disordered Perovskite Monte Carlo Optimizer
======================================================

A Python package to sample optimized atomic configurations of **disordered perovskites**
where **B-sites** are occupied by multiple chemical species (e.g., Mg/Nb/Ti).

The package uses a **Monte Carlo (MC)** approach:

- Propose a trial move by **swapping B-site species**;
- Optimize the new strucutre;
- Compute the **energy difference** between the old and new structures;
- Accept/reject using the **Metropolis criterion**;
- Repeat to obtain low-energy configurations and statistics.

.. toctree::
   :maxdepth: 2
   :caption: Documentation

   quickstart
   theory
   api
