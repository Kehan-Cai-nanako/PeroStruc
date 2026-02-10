# PeroStruc: Disordered Perovskite Monte Carlo Optimizer

A Python package aimed to sample optimized atomic configurations of **disordered perovskites**
where **B-sites** are occupied by multiple chemical species (e.g., Mg/Nb/Ti).

## Methodology

The package uses a **Monte Carlo (MC)** approach:

- Propose a trial move by **swapping B-site species**;
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
2. Build a force-field model for the perovskite.
3. Specify which sites are **B-sites** and the allowed B-site species.
4. Optimize the structure by calling `optimize.optimize()`, which conducts Monte Carlo swaps and records accepted structures / energies.
5. (Optional) Restart the optimization from a previous run by calling `optimize.restart()`.

## Example

The example demonstrates optimizing the structure of Pb(Mg<sub>1/3</sub>Nb<sub>2/3</sub>)O<sub>3</sub> (PMN), with a 12x12x12 simulation cell.
The **B-sites** are occupied by Mg/Nb atoms. (Replace the configuration and model with your own files.)

```python
from perostruc import optimize

optimize.optimize(filename = './PMN_conf.lmp',
                  model = "./PMN_model.pb",
                  temperature = 300.0, 
                  cell_size = np.array([12, 12, 12]), 
                  element_list = ['Mg', 'Nb', 'O', 'Pb'], 
                  B_site_elements = ['Mg', 'Nb'], 
                  neighbor_cutoff = 6,    # allow nearest and next-nearest neighbor swap.
                  traj_dir = 'trajs', 
                  logging_dir = 'logs', 
                  nepoch = 50, 
                  fmax = 0.03,            # force threshold for the FIRE optimization
                  fire_max_steps = 50,    # maximum number of steps for the FIRE optimization
                  )
```

Examples with more details can be found in the directory [`examples/`](examples/).

## Documentation

See the full instruction on [GitHub Pages](https://kehan-cai-nanako.github.io/PeroStruc/).