import os
from re import A
from time import time
import numpy as np
import glob
import ase
import pickle
from ase.optimize import FIRE
from ase.io import write, read
from ase.filters import FrechetCellFilter
from perostruc.utility import *


"""s
Simulation Parameters
"""

class FIRESwapOptimizer:
    """
    Monte Carlo optimizer for perovskite chemical structure optimization using FIRE-Swap methods
    """

    def __init__(self,
                 temperature: float,
                 A_site_elements: list[str],
                 B_site_elements: list[str],
                 element_specorder: list[str],
                 neighbor_cutoff: float = 6.0,
                 fmax: float = 0.03,
                 fire_max_steps: int = 50,
                 nepoch: int = 300000,
                 traj_dir: str = 'trajs',
                 logging_dir: str = 'logs',
                 mc_log_file: str = 'mc.log',
                 fire_log_file: str = 'fire.log',
                 ):
        ## thermodynamic parameters
        self.temperature = temperature
        self.kb = 8.617e-5  ## Boltzmann constant in eV/K
        self.kbT = self.kb * self.temperature
        ## system parameters
        self.atoms = None
        self.ncell = None
        self.A_site_elements = A_site_elements
        self.B_site_elements = B_site_elements
        self.element_specorder = element_specorder
        self.A_neighborlist = None
        self.B_neighborlist = None
        ## simulation parameters
        self.neighbor_cutoff = neighbor_cutoff
        self.fmax = fmax
        self.fire_max_steps = fire_max_steps
        self.start_epoch = 0
        self.current_swap = 0
        self.nepoch = nepoch
        ## IO parameters
        self.traj_dir = traj_dir
        self.logging_dir = logging_dir
        self.mc_log_file = os.path.join(self.logging_dir, mc_log_file)
        self.fire_log_file = os.path.join(self.logging_dir, fire_log_file)

    def _setup_dirs_and_log(self) -> None:
        os.makedirs(self.traj_dir, exist_ok=True)
        os.makedirs(self.logging_dir, exist_ok=True)
        with open(self.mc_log_file, 'w') as f:
            f.write('# Temperature = {:.3f}K, fmax={:.3f} eV/Å, fire_max_steps={}, neighbor_cutoff={}\n'.format(
                self.temperature, self.fmax, self.fire_max_steps, self.neighbor_cutoff))
            f.write('# iter, energy_before [eV], energy_after [eV], acceptance, time_cost, \n')

    def _build_swap_site_list(self, atoms: ase.Atoms):
        """Build neighbor list for A- and B-sites, respectively."""
        sym = atoms.get_chemical_symbols()
        natoms = len(sym)
        A_filter = np.array([_atm in self.A_site_elements for _atm in sym])
        B_filter = np.array([_atm in self.B_site_elements for _atm in sym])
        A_indices = np.arange(natoms)[A_filter]
        B_indices = np.arange(natoms)[B_filter]
        if len(A_indices) != self.ncell:
            raise ValueError('There are {} A-sites, but the number of A-sites should be {}'.format(
                len(A_indices), self.ncell))
        if len(B_indices) != self.ncell:
            raise ValueError('There are {} B-sites, but the number of B-sites should be {}'.format(
                len(B_indices), self.ncell))
        A_neighborlist = []
        B_neighborlist = []
        for asite_idx in A_indices:
            nb_idx, nb_sym, nb_dist = get_neighbor(atoms, asite_idx, cutoff=self.neighbor_cutoff)  ## this returns the indices ordered by distance from low to high
            nb_A_idx = []
            for _idx, _sym in zip(nb_idx, nb_sym):
                if _sym in self.A_site_elements:
                    nb_A_idx.append(_idx)
            if len(nb_A_idx) < 6:
                raise ValueError('Each A-site should have equal to or more than 6 nearest neighbors. Here it has {}'.format(len(nb_A_idx)))
            A_neighborlist.append(nb_A_idx)
        for bsite_idx in B_indices:
            nb_idx, nb_sym, nb_dist = get_neighbor(atoms, bsite_idx, cutoff=self.neighbor_cutoff)  ## this returns the indices ordered by distance from low to high
            nb_B_idx = []
            for _idx, _sym in zip(nb_idx, nb_sym):
                if _sym in self.B_site_elements:
                    nb_B_idx.append(_idx)
            if len(nb_B_idx) < 6:
                raise ValueError('Each B-site should have equal to or more than 6 nearest neighbors. Here it has {}'.format(len(nb_B_idx)))
            B_neighborlist.append(nb_B_idx)
        self.A_neighborlist = {
            'indices': A_indices,
            'neighborlist': A_neighborlist,
        }
        self.B_neighborlist = {
            'indices': B_indices,
            'neighborlist': B_neighborlist,
        }
        return None

    def save_neighborlist(self):
        """
        Save the neighborlist of the A- and B-site atoms to a pickle file.
        """
        with open('A_neighborlist.pkl', 'wb') as f:
            pickle.dump(self.A_neighborlist, f)
        with open('B_neighborlist.pkl', 'wb') as f:
            pickle.dump(self.B_neighborlist, f)
        return None
    
    def load_neighborlist(self):
        """
        Load the neighborlist of the A- and B-site atoms from a pickle file.
        """
        with open('A_neighborlist.pkl', 'rb') as f:
            self.A_neighborlist = pickle.load(f)
        with open('B_neighborlist.pkl', 'rb') as f:
            self.B_neighborlist = pickle.load(f)
        return None

    def get_lattice_map(self, type: str, latt_size: list[int]) -> np.ndarray:
        """
        Get the lattice map of the A- or B-site atoms.
        Parameters:
            type: str
                The type of the site to get the lattice map. Can be 'A' or 'B'.
            latt_size: list[int]
                The size of the lattice.  If this is a L1 X L2 X L3 supercell, then latt_size = [L1, L2, L3].
        Returns:
            lattice_map: np.ndarray of shape (l1, l2, l3) with integer value of indices of the site atoms
                The lattice map of the A- or B-site atoms.
        """
        ## sanity check
        if type == 'A':
            site_indices = self.A_neighborlist['indices']
        elif type == 'B':
            site_indices = self.B_neighborlist['indices']
        else:
            raise ValueError('The type should be either "A" or "B".')
        assert len(latt_size) == 3, 'The lattice size should be a list of 3 integers.'
        l1, l2, l3 = latt_size
        assert l1 * l2 * l3 == len(site_indices), 'The number of site indices should be equal to the product of the lattice size.'
        ## get the lattice map
        lattice_map = np.zeros((l1, l2, l3), dtype=int) - 1  ## -1 means the lattice site is empty
        scaled_pos = self.atoms.get_scaled_positions()[site_indices]
        assert scaled_pos.min() >= 0 and scaled_pos.max() <= 1, 'the scaled positions of the site atoms are not in the unit cell'
        origin = scaled_pos.min(axis=0)
        scaled_pos = scaled_pos - origin
        latt_i = np.round(scaled_pos[:,0] * l1).astype(int) % l1
        latt_j = np.round(scaled_pos[:,1] * l2).astype(int) % l2
        latt_k = np.round(scaled_pos[:,2] * l3).astype(int) % l3
        for site_idx, _i, _j, _k in zip(site_indices, latt_i, latt_j, latt_k):
            lattice_map[_i,_j,_k] = site_idx
        if lattice_map.min() == -1:
            raise ValueError('some lattice sites are empty. Please check if the initial configuration is correct.')
        return lattice_map

    def initialize(self, supercell_size: list[int], init_config=None, calculator=None, restart: bool = False, build_neighborlist: bool = True) -> None:
        """
        Initialize optimization from scratch. The given initial configuration is a ASE Atoms object. It should already has an calculator.
        Parameters:
            init_config: ase.Atoms
                The initial configuration.
            supercell_size: list[int]
                The size of the supercell in perovskite conventional.
            calculator: ase.calculators.calculator.Calculator, optional if the initial configuration already has a calculator.
                The calculator to use for the optimization.
        """
        if restart is False:
            if init_config is None:
                raise ValueError('The initial configuration is not provided. Please provide a ase.Atoms object or set restart to True.')
            self.atoms = init_config
            self._setup_dirs_and_log()
        else:
            if not os.path.exists(self.traj_dir):
                raise FileNotFoundError('The trajectory directory {} does not exist. Cannot restart from a non-existent trajectory directory.'.format(self.traj_dir))
            if not os.path.exists(self.mc_log_file):
                raise FileNotFoundError('The log file {} does not exist. Cannot restart from a non-existent log file.'.format(self.mc_log_file))
            mc_log = np.loadtxt(self.mc_log_file, skiprows=2, delimiter=',')[-1]
            last_epoch = int(mc_log[0])
            self.start_epoch = last_epoch + 1
            print('Restarting from the result of the last epoch: {}'.format(last_epoch))
            traj_files = glob.glob(os.path.join(self.traj_dir, 'f*.lmp'))
            traj_files.sort(key=lambda x: int(os.path.splitext(os.path.basename(x))[0][1:]))
            last_traj_file = traj_files[-1]
            self.current_swap = int(os.path.splitext(os.path.basename(last_traj_file))[0][1:])
            print('Loading the last traj file: {}'.format(last_traj_file))
            atoms = read(last_traj_file, format='lammps-data', atom_style='atomic')
            update_element(atoms, self.element_specorder)
            self.atoms = atoms
        ## set calculator
        if self.atoms.calc is None:
            if calculator is None:
                raise ValueError('The initial configuration does not have a calculator, and you did not provide a calculator.')
            else:
                self.atoms.calc = calculator
                print('Atoms loaded. The calculator is set')
        else:
            if calculator is not None:
                self.atoms.calc = calculator
                print('The existing calculator is overwritten by the provided calculator. Please check if you want to do this.')
        ## sanity check
        if len(supercell_size) != 3:
            raise ValueError('Supercell size should be a list of 3 integers.')
        self.ncell = int(np.prod(supercell_size))
        if len(self.atoms) != self.ncell * 5:
            raise ValueError('The number of atoms in the initial configuration is not equal to the product of the perovskite supercell size times 5.')
        ## initialize
        if build_neighborlist is True:
            self._build_swap_site_list(self.atoms)
        return None

    def run_FIRE(self, fmax=0.02, steps=1000, flexible_cell: bool = False) -> None:
        """
        Do a single FIRE optimization, used to relax the initial configuration. We will use a more strict force threshold and max steps here.
        Parameters:
            fmax: float
                The force threshold for the FIRE optimization.
            steps: int
                The maximum number of steps for the FIRE optimization.
        """
        if self.atoms is None:
            raise ValueError('The atoms are not initialized. Please initialize the atoms first.')
        if self.atoms.calc is None:
            raise ValueError('The atoms do not have a calculator. Please set a calculator first.')
        if flexible_cell is True:
            atoms_flexible = FrechetCellFilter(self.atoms)
            dyn = FIRE(atoms_flexible, logfile=self.fire_log_file)
        else:
            dyn = FIRE(self.atoms, logfile=self.fire_log_file)
        dyn.run(fmax=fmax, steps=steps)
        return None

    def optimize(self, swap_type: str, flexible_cell: bool = False, profile: bool = False) -> None:
        """
        Run MC swap loop from start_epoch to self.nepoch.
        Parameters:
            swap_type: str
                The type of the swap to be performed. Can be 'A' or 'B'.
            flexible_cell: bool
                Whether to use flexible cell optimization.
        """
        ## sanity check
        if swap_type == 'A':
            site_indices = self.A_neighborlist['indices']
            site_neighborlist = self.A_neighborlist['neighborlist']
        elif swap_type == 'B':
            site_indices = self.B_neighborlist['indices']
            site_neighborlist = self.B_neighborlist['neighborlist']
        else:
            raise ValueError('The swap type should be either "A" or "B".')
            
        ## start the MC swap loop
        penergy = self.atoms.get_potential_energy()
        for i in range(self.start_epoch, self.nepoch):
            if profile:
                t0 = time()

            ## determine which pair of atoms to swap
            choose_seed = np.random.choice(np.arange(self.ncell))
            atom1_idx = site_indices[choose_seed]
            atom2_idx = np.random.choice(site_neighborlist[choose_seed])
            sym = self.atoms.get_chemical_symbols()
            if sym[atom1_idx] == sym[atom2_idx]:
                continue

            ## make a new ase.Atoms object with the swapped atoms
            sym_new = sym.copy()
            sym_new[atom1_idx] = sym[atom2_idx]
            sym_new[atom2_idx] = sym[atom1_idx]
            atoms_new = self.atoms.copy()
            atoms_new.set_chemical_symbols(sym_new)
            atoms_new.calc = self.atoms.calc

            ## optimize the new atoms
            if flexible_cell is True:
                atoms_new_flexible = FrechetCellFilter(atoms_new)
                dyn = FIRE(atoms_new_flexible, logfile=self.fire_log_file)
            else:
                dyn = FIRE(atoms_new, logfile=self.fire_log_file)
            dyn.run(fmax=self.fmax, steps=self.fire_max_steps)
            penergy_new = atoms_new.get_potential_energy()

            ## compare the energy of the new and old atoms, and decide whether to accept the swap
            energy_diff = penergy_new - penergy
            if energy_diff < 0 or (np.random.rand() < np.exp(-energy_diff / self.kbT)):
                accept = 1
                self.atoms = atoms_new
                penergy = penergy_new
                self.current_swap += 1
                write(
                    os.path.join(self.traj_dir, 'f{}.lmp'.format(self.current_swap)), 
                    self.atoms, 
                    format='lammps-data', 
                    specorder=self.element_specorder)
            else:
                accept = 0
            
            ## logging
            if profile:
                t1 = time()
                print('========== epoch={},  swap-{}, accept?:{} time cost={:.3f}s =========='.format(i, self.current_swap, accept, t1 - t0))
                print('attempt to swap {} and {} failed, dE={:.3f}eV, '.format(
                    sym[atom1_idx], sym[atom2_idx], energy_diff))
            with open(self.mc_log_file, 'a') as f:
                f.write('{}, {}, {}, {}, {}\n'.format(i, penergy, penergy_new, accept, t1 - t0))
 
    def get_pair_product(self, atoms: ase.Atoms, score_map: dict[str, float], type: str = 'A') -> float:
        """
        Calculate the order parameter O = \sum_{i\sim j} s_i s_j / N /6, where s_i is the score of the i-th atom, and N is the number of cells. 
        The sum is over all pairs of nearest-neighbor site-A (or B, specified by `type`) atoms (there are 6 nearest neighbors for each site).
        score_map is a dictionary that maps the element on A (or B, specified by `type`) site to the score. For example, {'Pb': -1.0, 'Sr': 1.0} for PST.
        """
        if type == 'A':
            elements = self.A_site_elements
            site_indices = self.A_neighborlist['indices']
            site_neighborlist = self.A_neighborlist['neighborlist']
        elif type == 'B':
            elements = self.B_site_elements
            site_indices = self.B_neighborlist['indices']
            site_neighborlist = self.B_neighborlist['neighborlist']
        else:
            raise ValueError('The type should be either "A" or "B".')
        for element in elements:
            assert element in score_map, 'The element {} is not in the score map'.format(element)
        order_parameter = 0.0
        syms = atoms.get_chemical_symbols()
        for site_idx, nb_indices in zip(site_indices, site_neighborlist):
            sym_center = syms[site_idx]
            score_center = score_map[sym_center]
            if len(nb_indices) < 6:
                raise ValueError('Each site should have equal to or more than 6 nearest neighbors. Here it has {}'.format(len(nb_indices)))
            else:
                nnb_indices = nb_indices[:6]  ## note that the indices are already ordered by distance from low to high
            for nb_idx in nnb_indices:
                sym_nb = syms[nb_idx]
                score_nb = score_map[sym_nb]
                order_parameter += score_center * score_nb
        order_parameter = order_parameter / len(site_indices) / 6.0
        return order_parameter
 