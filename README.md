# PeroStruc: Disordered Perovskite Monte Carlo Optimizer

A Python package aimed to sample chemically disordered atomic configurations of **complex perovskites** (e.g., Pb(Mg<sub>1/3</sub>Nb<sub>2/3</sub>)O<sub>3</sub>) and perovskite solid solutions (e.g., Pb<sub>x</sub>Sr<sub>1-x</sub>TiO<sub>3</sub>)

## Methodology

The package uses a hybrid **Molecular Dynamics (MD)** and **Monte Carlo (MC)** approach called **FIRESwap**

- Propose a trial move by **swapping A-site (or B-site) species**;
- Optimize the new strucutre;
- Compute the **energy difference** between the old and new structures;
- Accept/reject using the **Metropolis criterion**;
- Repeat to obtain low-energy configurations and statistics.

## Installation

For quick install:

```bash
python -m pip install .
```

For development (editable install):

```bash
python -m pip install -e .
```

## Minimal workflow

A typical workflow is:

1. Build an initial perovskite structure (e.g., ABO<sub>3</sub>).
2. Specify a force-field model for the perovskite. This step is not handled by this package. Please use your favorite force-field model. Note that this package can work with any force-field model as long as it can be used as a ASE calculator. For example, you can use the Deep Potential (DP) (https://github.com/deepmodeling/deepmd-kit) model, which already has a ASE calculator interface.
3. Use the `FIRESwapOptimizer` class to initialize the optimizer. You should specify the initial configuration, the force-field model, and the simulation parameters.
4. Optimize the structure by calling `optimize()`, which conducts Monte Carlo swaps and records accepted structures / energies.
5. (Optional) Restart the optimization from a previous run by using the `restart()` method.

## Example

The following example demonstrates optimizing the structure of Pb(Mg<sub>1/3</sub>Nb<sub>2/3</sub>)O<sub>3</sub> (PMN), with a 6x6x6 simulation cell.

```python
import numpy as np
import ase
import ase.io
from perostruc import FIRESwapOptimizer, update_element
from deepmd.calculator import DP    ## you need to install deepmd-kit to use the DP model. deepmd-kit is not included in this package.

init_config = ase.io.read('initial_configs/disordered_L6X6X6.lmp', format='lammps-data', atom_style='atomic')

## this is the order of the elements in the lammps data file. We need this because lammps dump file does not have the element names. Instead, it uses implcit numbering (1, 2, 3, ...). So we need to specify how the implicit numbering corresponds to the element names.
element_specorder = ['Mg','Nb','O','Pb']   
update_element(init_config,element_specorder)

## load the force-field model. You can use any other force-field model as long as it has a ASE calculator interface. 
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

## initialize from scratch
simulation.initialize(supercell_size = [6, 6, 6], init_config = init_config)
## relax the initial geometric configuration using plain FIRE.
simulation.run_FIRE(fmax = 0.02, steps = 1000)
## conduct the FIRESwap optimization, which alternates between FIRE and Monte Carlo swaps.
## swap_type = 'B' means swapping B-site species. You can also set it to 'A' to swap A-site species.
## flexible_cell = False means fixing the cell parameters. You can set it to True to allow the cell parameters to change.
simulation.optimize(swap_type = 'B', flexible_cell = False, profile = True)
```

Examples with more details can be found in the directory [`examples/`](examples/).

## Documentation

See the full instruction on [GitHub Pages](https://kehan-cai-nanako.github.io/PeroStruc/).

## TODO

1. Make Breadth-first search (BFS) a module. 
2. Logo.
3. Documentation. 
4. PMN analyze Pb order parameter.
