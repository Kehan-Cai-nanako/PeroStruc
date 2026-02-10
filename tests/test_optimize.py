import numpy as np
from perostruc import optimize

if __name__ == '__main__':
    flag_test_optimize = 1
    if flag_test_optimize != 0:
        optimize.optimize(filename = './initial_configs/sublattice_L12X12X12.lmp',
                          model = "./unipero.pb",
                          temperature = 300.0, 
                          cell_size = np.array([12, 12, 12]), 
                          element_list = ['Mg', 'Nb','O','Pb'], 
                          B_site_elements = ['Mg', 'Nb'], 
                          neighbor_cutoff = 6,    # allow nearest and next-nearest neighbor swap.
                          traj_dir = 'trajs', 
                          logging_dir = 'logs', 
                          nepoch = 2, 
                          fmax = 0.03,            # force threshold for the FIRE optimization
                          fire_max_steps = 50,    # maximum number of steps for the FIRE optimization
                          )
    
    flag_test_restart = 0
    if flag_test_restart != 0:
        optimize.restart(model = "./unipero.pb",
                         temperature = 300.0, 
                         cell_size = np.array([12, 12, 12]), 
                         element_list = ['Mg', 'Nb','O','Pb'], 
                         B_site_elements = ['Mg', 'Nb'], 
                         neighbor_cutoff = 6,    # allow nearest and next-nearest neighbor swap.
                         traj_dir = 'trajs', 
                         logging_dir = 'logs', 
                         nepoch = 10, 
                         fmax = 0.03,          # force threshold for the FIRE optimization
                         fire_max_steps = 50,    # maximum number of steps for the FIRE optimization
                         )