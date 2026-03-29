PeroStruc: Disordered Perovskite Monte Carlo Optimizer
======================================================

A Python package for sampling compositional and geometric structures of **complex perovskites** (e.g., :math:`\mathrm{Pb\left(Mg_{1/3}Nb_{2/3}\right)O_3}`) 
and **perovskite solid solutions** (e.g., :math:`\mathrm{Pb_{x}Sr_{1-x}TiO_3}`).

The package uses a hybrid **Molecular Dynamics (MD)** and **Monte Carlo (MC)** approach called **FIRESwap**:

- Propose a trial move by **swapping A-site (or B-site) species**;
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
