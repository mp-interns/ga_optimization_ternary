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
    
    dirs = ["../data/abn3", "../data/abo2f", "../data/abo2n", "../data/abo2s", "../data/abo3", "../data/abofn", "../data/abon2"]
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
            
            anion_dict = {}
            anion_dict["O3"] = 0
            anion_dict["O2N"] = 1
            anion_dict["ON2"] = 2
            anion_dict["N3"] = 3
            anion_dict["O2F"] = 4
            anion_dict["OFN"] = 5
            anion_dict["O2S"] = 6
            
            out_dict['anion_idx'] = anion_dict[out_dict['anion']]
            out_dict['is_direct'] = True if out_dict['gllbsc_dir-gap'] == out_dict['gllbsc_ind-gap'] else False
            out_dict['sum_magnetic_moments'] = sum(entry['MagneticMoments'])
             
            processed_coll.insert(out_dict)
        except:
            print entry


def fix_db_raw():
    db_raw = data_raw_coll()
    update_cnt = 0
    with open(os.path.join("../data","abo2f.txt")) as f:
        for line in f:
            f_name, energy = line.split()[0]+".json", float(line.split()[1])
            with open(os.path.join("../data/abo2f", f_name)) as m_file:
                data = json.loads(m_file.read())
                print data['A'], data['B'], data['anion']
                db_raw.update({"A": data['A'], "B": data['B'], "anion": data['anion']}, {"$set": {"heat_of_formation_all": energy}})
                m_file.close()
                update_cnt += 1
                print update_cnt
    
    print update_cnt

def fudge_gap_in_hf_o2f():
    db_raw = data_raw_coll()
    db_raw.update({"A": "In", "B": "Hf", "anion": "O2F"}, {"$set": {"gllbsc_ind-gap": 3.01}})
    
if __name__ == "__main__":
    # init_db_raw()
    # fix_db_raw()
    fudge_gap_in_hf_o2f()
    init_db_processed()
    
    #find_atomic_range()
