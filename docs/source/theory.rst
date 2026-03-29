Method Overview
===============

Disordered perovskites
----------------------

We consider complex perovskites or perovskite solid solutions of the form ABO\ :sub:`3`, where the **A (or B) sublattice** is
occupied by a mixture of chemical species (e.g., Pb/Sr or Mg/Nb/Ti).

Monte Carlo sampling
--------------------

At each MC step:

1. Select two A-sites (or B-sites) and propose a **swap attempt** of their species.
2. Evaluate the energy difference:

   .. math::

      \Delta E = E_{\mathrm{new}} - E_{\mathrm{old}}

3. Accept the move with Metropolis probability:

   .. math::

      P(\mathrm{accept}) =
      \begin{cases}
      1, & \Delta E \le 0 \\
      \exp(-\Delta E / k_B T), & \Delta E > 0
      \end{cases}

where :math:`T` is the temperature and :math:`k_B` is Boltzmann's constant.

Practical notes
---------------

- The "energy" may come from an empirical potential, cluster expansion,
  machine-learning potential, or any callable energy model in principle.
- Swaps conserve overall composition.
