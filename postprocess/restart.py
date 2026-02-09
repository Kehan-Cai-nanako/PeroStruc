import os
from time import time
import numpy as np
import glob
import csv
from ase.optimize import FIRE
from ase.io import write
from deepmd.calculator import DP
from utility import *


"""
Simulation Parameters
"""
kb = 8.617e-5  ## Boltzmann constant in eV/K


def optimize_struc(model: str = "unipero.pb",
                   temperature: float = 300.0, 
                   cell_size = np.array([12, 12, 12]), 
                   element_list: list[str] = ['Mg', 'Nb','O','Pb'], 
                   B_site_elements: list[str] = ['Mg', 'Nb'], 
                   neighbor_cutoff: int = 6,    # allow nearest and next-nearest neighbor swap.
                   traj_dir: str = 'trajs', 
                   logging_dir: str = 'logs', 
                   nepoch: int = 300000, 
                   fmax: float = 0.03,          # force threshold for the FIRE optimization
                   fire_max_steps: int = 50,    # maximum number of steps for the FIRE optimization
                   ):
    
    kbT = kb * temperature
    ncell = np.prod(cell_size)

    """
    Preprocessing
    """
    ### Create directories for saving the figures and trajectories.
    os.makedirs(traj_dir, exist_ok=True)
    os.makedirs(logging_dir, exist_ok=True)
    mc_log_file = os.path.join(logging_dir, 'mc.log')
    assert os.path.exists(mc_log_file), 'The log file {} does not exist'.format(mc_log_file)
    mc_log = np.loadtxt(mc_log_file, skiprows=2, delimiter=',')[-1]
    last_epoch = int(mc_log[0])
    print('Last epoch: {}'.format(last_epoch))

    ### Load the last traj file.
    traj_files = glob.glob(os.path.join(traj_dir, 'f*.lmp'))
    traj_files.sort(key=lambda x: int(os.path.splitext(os.path.basename(x))[0][1:]))
    last_traj_file = traj_files[-1]
    nswap = int(os.path.splitext(os.path.basename(last_traj_file))[0][1:])
    print('Loading the last traj file: {}'.format(last_traj_file))
    perov = ase.io.read(last_traj_file, format='lammps-data', atom_style='atomic')
    update_element(perov, element_list)
    natoms = len(perov)
    print(perov)

    ### Get the filter for B-site atoms (Mg or Nb). Although Nb and Mg will be swapped, the filter for B-site won't change.
    sym = perov.get_chemical_symbols()
    B_filter = np.array([ _atm in B_site_elements for _atm in sym ])
    Bsites_indices = np.arange(natoms)[B_filter]   # indices of all B-site
    if B_filter.astype(int).sum() != ncell:
        raise ValueError('There are {} B-sites, but the number of B-sites should be {}'.format(B_filter.astype(int).sum(), ncell))

    ### For each B-site, we will find the index of its 6 nearest neighbors. 
    ### Again, after swapping, the neighboring relation won't change. Only the element type of B-site atoms will be changed.
    Bsites_neighborlist = []
    for bsite_idx in Bsites_indices:
        nb_idx, nb_sym, nb_dist = get_neighbor(perov, bsite_idx, cutoff=neighbor_cutoff)  ## sorted by distance from low to high.  
        nb_B_filter = np.array([ _atm in B_site_elements for _atm in nb_sym ])
        nb_idx = nb_idx[nb_B_filter]
        nb_sym = nb_sym[nb_B_filter]
        nb_dist = nb_dist[nb_B_filter]
        if nb_idx.size <6:
            raise ValueError('Each B-site should have equal to or more than 6 nearest neighbors. Here it has {}'.format(nb_idx.size))
        Bsites_neighborlist.append(nb_idx)

    ### Describe the interatomic interactions with DP model.
    dpmodel = DP(model=model)
    perov.calc = dpmodel
    print("initial E={}eV/atom".format(perov.get_potential_energy()/natoms))
    penergy = [ perov.get_potential_energy() ]

    """
    MC-MD simulation
    """
    ### Create list to save the energy.
    before_energies = []
    after_energies = []
    for i in range(last_epoch+1, nepoch):
        ### Randomly choose one Mg/Nb atom and get its neighborlist.
        choose_bsite_seed = np.random.choice(np.arange(ncell))
        atom1_idx = Bsites_indices[choose_bsite_seed]
        atom2_idx = np.random.choice(Bsites_neighborlist[choose_bsite_seed])
        if sym[atom1_idx] == sym[atom2_idx]:
            continue
        else:
            t0 = time()
            before_energies.append(penergy[-1])
            sym = perov.get_chemical_symbols()
            sym_new =  sym.copy()
            sym_new[atom1_idx] = sym[atom2_idx]
            sym_new[atom2_idx] = sym[atom1_idx]
            perov_new = perov.copy()
            perov_new.set_chemical_symbols(sym_new)
            perov_new.calc = dpmodel
            dyn = FIRE(perov_new, logfile='{}/fire.log'.format(logging_dir) )
            dyn.run(fmax=fmax, steps=fire_max_steps)   ## limit the number of steps for the FIRE optimization
            penergy_new = perov_new.get_potential_energy()
            after_energies.append(penergy_new)
            energy_diff = penergy_new - penergy[-1]
            if energy_diff < 0 or (np.random.rand() < np.exp(-energy_diff/kbT)):
                accept = 1
                perov = perov_new
                penergy.append(penergy_new)
                nswap += 1
                ase.io.write('./{}/f{}.lmp'.format(traj_dir, nswap), perov, format='lammps-data')
                ### logging
                t1 = time()
                print('========== epoch={},  swap-{},  time cost={:.3f}s =========='.format(i, nswap, t1-t0))
                print('attempt to swap {} and {} succeeded, dE={:.3f}eV, '.format(
                    sym[atom1_idx], sym[atom2_idx],energy_diff  ))
            else:
                accept = 0
                t1 = time()
                print('========== epoch={},  failed swap,  time cost={:.3f}s =========='.format(i, t1-t0))
                print('attempt to swap {} and {} failed, dE={:.3f}eV, '.format(
                    sym[atom1_idx], sym[atom2_idx],energy_diff  ))
            with open(mc_log_file, 'a') as f:
                f.write('{}, {}, {}, {}, {}\n'.format(i, before_energies[-1], after_energies[-1], accept, t1-t0))
    
    return 0