import ase
import ase.io
from ase import units
import numpy as np


def update_element(atoms, element_list):
    sym = atoms.get_atomic_numbers()
    sym_new = []
    for s in sym:
        sym_new.append(element_list[s-1])
    atoms.set_chemical_symbols(sym_new)
    return


def printenergy(a):
    """Function to print the potential, kinetic and total energy"""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))


def get_neighbor(atoms, idx, cutoff=6):
    d = atoms.get_distances(idx, np.arange(len(atoms)), mic=True, vector=False)
    neighbor_filter = d < cutoff
    neighbor_idx = np.arange(len(atoms))[neighbor_filter]
    neighbor_sym = np.array(atoms.get_chemical_symbols())[neighbor_idx]
    neighbor_dist = d[neighbor_filter]
    neighbor_dist_sort = np.argsort(neighbor_dist)

    ## sort by distance from low to high. Then remove the first, which is itself
    neighbor_idx = neighbor_idx[neighbor_dist_sort][1:]
    neighbor_sym = neighbor_sym[neighbor_dist_sort][1:]
    neighbor_dist = neighbor_dist[neighbor_dist_sort][1:]
    return  neighbor_idx, neighbor_sym, neighbor_dist



"""
Helper functions for cluster analysis
"""


def get_neighbor_indices(i, j, k, l1, l2, l3):
    return [( (i+1) % l1, j, k), ( (i-1) % l1, j, k), (i, (j+1) % l2, k), (i, (j-1) % l2, k), (i, j, (k+1) % l3), (i, j, (k-1) % l3)]

def find_neighbor_values(lattice_values, i, j, k):
    l1, l2, l3 = lattice_values.shape
    neighbors = get_neighbor_indices(i, j, k, l1, l2, l3)
    neighbor_values = []
    for ni, nj, nk in neighbors:
        neighbor_values.append(lattice_values[ni, nj, nk])
    return neighbor_values

def BFS_cluster_analysis(cluster_flag):
    l1, l2, l3 = cluster_flag.shape
    visited = np.zeros_like(cluster_flag, dtype=bool)
    all_clusters = []

    for i in range(l1):
        for j in range(l2):
            for k in range(l3):
                if cluster_flag[i][j][k] == 1 and visited[i][j][k] == False:
                    # New cluster found, start BFS
                    current_cluster = []
                    queue = [(i, j, k)]
                    visited[i][j][k] = True

                    while len(queue) > 0:
                        ci, cj, ck = queue.pop(0)
                        current_cluster.append((ci, cj, ck))

                        # Define neighbors again for traversal
                        neighbors = get_neighbor_indices(ci, cj, ck, l1, l2, l3)

                        for ni, nj, nk in neighbors:
                            if cluster_flag[ni][nj][nk] == 1 and visited[ni][nj][nk] == False:
                                visited[ni][nj][nk] = True
                                queue.append((ni, nj, nk))
                    
                    all_clusters.append(current_cluster)

    return all_clusters


def calculate_cluster_surface_size(cluster_flag):
    '''
    Calculate the surface size of a cluster
    Parameters:
    cluster_flag: np.ndarray[int] with shape (l1, l2, l3). 1 for in the cluster, 0 for out of the cluster.  
        the flag of the cluster
    Returns:
    surface_size: int
        the surface size of the cluster
    '''
    surface_size = 0
    for direction in [0,1,2]:
        flag_shifted = np.roll(cluster_flag, shift=1, axis=direction)
        surface_size += np.sum((cluster_flag != flag_shifted).astype(int))
    return surface_size


def find_clusters(atoms, lattice_map, element:str, tolerance:int = 1):
    """
    Find the clusters of the given element in the trajectory.
    Parameters:
        traj: list[ase.Atoms]
            The trajectory of the atoms.
        lattice_map: np.ndarray[int]
            The lattice map of the atoms. The integer value is the index of the site atom.
        element: str
            The element to find the clusters.
        tolerance: int
            The tolerance for a site to be considered as part of the cluster. The number of neighboring elements that are not the given element.
    Returns:
        all_cluster_sum: int
            The sum of the sizes of all clusters of the given element. The default tolerance is 1, which means the site is considered as part of the cluster if it has at most one neighboring element that is not the given element.
        largest_cluster_size: int
            The size of the largest cluster of the given element.
        surface_volume_ratio: float
            The surface-volume ratio of the largest cluster of the given element.
    """
    l1, l2, l3 = lattice_map.shape
    ## initialize the lattice values (0 means Mg, 1 means Nb)
    lattice_values = np.zeros_like(lattice_map, dtype=int)
    syms = atoms.get_chemical_symbols()
    for i in range(l1):
        for j in range(l2):
            for k in range(l3):
                site_sym = syms[lattice_map[i,j,k]]
                lattice_values[i,j,k] = 1 if site_sym == element else 0

    ## find clusters (neighboring at most one another element)
    cluster_flag = np.zeros_like(lattice_values, dtype=int)
    for i in range(l1):
        for j in range(l2):
            for k in range(l3):
                if lattice_values[i,j,k] != 0:
                    neighbor_values = find_neighbor_values(lattice_values, i, j, k)
                    if neighbor_values.count(0) <= tolerance:
                        cluster_flag[i,j,k] = 1
    all_cluster_sum = cluster_flag.sum()
    ## get the clusters and the size of the largest cluster
    all_clusters = BFS_cluster_analysis(cluster_flag)
    ## if the list of clusters is empty, return 0 for all_cluster_sum, largest_cluster_size, and surface_volume_ratio
    if len(all_clusters) == 0:
        return 0, 0, 0
        
    largest_cluster = all_clusters[np.argmax([len(cluster) for cluster in all_clusters])]
    largest_cluster_size = len(largest_cluster)

    ## 
    largest_cluster_flag = np.zeros_like(lattice_values, dtype=int)
    for i, j, k in largest_cluster:
        largest_cluster_flag[i,j,k] = 1
    surface_size = calculate_cluster_surface_size(largest_cluster_flag)
    surface_volume_ratio = surface_size / largest_cluster_size
    
    return all_cluster_sum, largest_cluster_size, surface_volume_ratio
         