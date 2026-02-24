import os
import sys
import warnings
import numpy as np
import ase
import ase.io
from ase.geometry.analysis import Analysis
from ase.neighborlist import NeighborList, natural_cutoffs, NewPrimitiveNeighborList  ## use NewPrimitiveNeighborList <- linear scaling
from ase.geometry import get_distances



class PST_perovskite:
    def __init__(self, ABCO3=['Pb','Sr','Ti','O']):
        if len(ABCO3) == 4:
            self.type_A, self.type_B, self.type_C, self.type_O = ABCO3
        else:
            raise ValueError("Please specify the atom types of perovskites with the conventional order, e.g. ['Pb','Mg','Nb','O'] for lead magnesium niobate")
  
    def _create_disordered_domain_small(self, supercell=[3,3,3], a=4.0 ):
        '''
        Building disordered PST configurations perserving the stoichiometry.
        Args:
            supercell: the size of the supercell. Should be small.
            a: lattice constant in Angstrom
        '''
        conf = ase.Atoms( 
            symbols= self.type_A+self.type_C+3*self.type_O,
            scaled_positions=[ 
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.5],
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5]],
            cell=[a, a, a],
            pbc=True)
        conf = conf.repeat(supercell)
        ncell = supercell[0]*supercell[1]*supercell[2]
        syms = conf.get_chemical_symbols()
        Asites = [idx for idx, s in enumerate(syms) if s==self.type_A]
        ndopants = ncell*2//3
        replaced = np.random.choice(Asites, ndopants, replace=False)
        for idx in replaced:
            syms[idx] = self.type_B
        conf.set_chemical_symbols(syms)
        return conf
      

if __name__ == "__main__":
    pst_factory  = PST_perovskite(['Pb','Sr','Ti','O'])
    type_map = ['Pb','Sr', 'O', 'Ti']  #  alphabatic order

    supercell = [6,6,6]
    disordered = pst_factory._create_disordered_domain_small( supercell=[6,6,6], a=4.00 )
    ase.io.write('Asite_disordered_L{}X{}X{}.lmp'.format(supercell[0], supercell[1], supercell[2]), disordered, format='lammps-data', specorder=type_map)

    ## ordered configuration
    supercell = [12,12,12]
    sublattice = pst_factory._create_disordered_domain_small( supercell=[12,12,12], a=4.00 )
    ase.io.write('Asite_disordered_L{}X{}X{}.lmp'.format(supercell[0], supercell[1], supercell[2]), sublattice, format='lammps-data', specorder=type_map)