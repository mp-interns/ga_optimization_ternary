#!/usr/bin/env python

'''
Created on Mar 13, 2012
'''

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Mar 13, 2012"

import pymongo
import os
import json
from pymatgen.core.periodic_table import Element

from pymatgen.core.structure import Structure, Lattice

connection = pymongo.Connection('localhost', 12345)
_unc_db = connection.unc


def data_raw_coll(remove=False):
    coll = _unc_db.data_raw
    if remove:
        coll.remove()
    return coll


def data_process_coll(remove=False):
    coll = _unc_db.data_process
    if remove:
        coll.remove()
    return coll


def init_db_raw():
    raw_coll = data_raw_coll(remove=True)
    
    dirs = ["abo3", "abo2n"]
    for dir in dirs:
        for f in os.listdir(dir):
            if  '.json' in f:
                file = open(os.path.join(dir,f))
                data = json.loads(file.read())
                print type(data), f
                data['filename'] = os.path.join(dir,f)
                raw_coll.insert(data)
                file.close()

            
def init_db_processed():
    raw_coll = data_raw_coll()
    processed_coll = data_process_coll(remove=True)
    
    for entry in raw_coll.find():
        try:
            out_dict = {}
            copy_fields = ['filename', 'A', 'B', 'anion', 'CB_dir', "CB_ind", "VB_dir", "VB_ind", "FermiLevel", "FermiWidth", "XCFunctional", "gllbsc_dir-gap", "gllbsc_ind-gap", "heat_of_formation_all"]
            for field in copy_fields:
                out_dict[field] = entry[field]
            
            #structure
            latt = Lattice(entry['ase_cell'])
            atomic_species = entry['ase_chemical_symbols']
            positions = entry['ase_scaled_positions']
            s = Structure(latt, atomic_species, positions, True)
            out_dict['structure'] = s.to_dict
            
            out_dict['contains_anion'] = True if out_dict['anion'].strip() == "O2N" else False
            out_dict['is_direct'] = True if out_dict['gllbsc_dir-gap'] == out_dict['gllbsc_ind-gap'] else False
            out_dict['sum_magnetic_moments'] = sum(entry['MagneticMoments'])
             
            processed_coll.insert(out_dict)
        except:
            print entry
    
if __name__ == "__main__":
    #init_db_raw()
    #init_db_processed()
    find_atomic_range()
    """
        def __init__(self, lattice, atomicspecies, coords, validate_proximity = False, to_unit_cell = False, coords_are_cartesian = False):
        
        Create a periodic structure.
        
        Arguments:
            lattice:
                pymatgen.core.lattice Lattice object signify the lattice.
            atomicspecies:
                list of atomic species.  dict of elements and occupancies.
            fractional_coords:
                list of fractional coordinates of each species.
            validate_proximity:
                Whether to check if there are sites that are less than 1 Ang apart. Defaults to false.
            coords_are_cartesian:
                Set to True if you are providing coordinates in cartesian coordinates. Defaults to false.
        """