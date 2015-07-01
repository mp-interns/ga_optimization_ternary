 #!/usr/bin/env python

'''
Created on Mar 13, 2012
'''
from pymatgen import Structure, Lattice

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Mar 13, 2012"

import pymongo
import os
import json

LOCAL_MONGO_PORT = int(os.environ.get('LOCAL_MONGO_PORT', 27017))
connection = pymongo.Connection('localhost', LOCAL_MONGO_PORT)
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


def init_db_raw(dir=None):
    raw_coll = data_raw_coll(remove=True)
    if not dir:
        data_dirs = [dir for dir in os.listdir('.') if dir.endswith('_data')]
        if len(data_dirs) == 1:
            dir = data_dirs[0]
        else:
            raise ValueError("Multiple '_data' directories. Specify path.")
    for fname in os.listdir(dir):
        if fname.endswith('.json'):
            with open(os.path.join(dir,fname)) as f:
                data = json.loads(f.read())
                print type(data), fname
                data['filename'] = os.path.join(dir,fname)
                raw_coll.insert(data)


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
            print out_dict
        except:
            print entry


def fix_db_raw():
    db_raw = data_raw_coll()
    update_cnt = 0
    with open(os.path.join("/Users/ajain/Documents/papers_proposals/2012Paper_UNCOptimization/data/rawdata","abo2f.txt")) as f:
        for line in f:
            f_name, energy = line.split()[0]+".json", float(line.split()[1])
            with open(os.path.join("/Users/ajain/Documents/papers_proposals/2012Paper_UNCOptimization/data/rawdata/abo2f", f_name)) as m_file:
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

def export():
    processed_coll = data_process_coll()

    for i in processed_coll.find():
        del i['_id']
        f_name = os.path.join("alldata", "{}{}{}.json".format(i['A'], i['B'],i['anion']))
        with open(f_name, 'w') as f:
            f.write(json.dumps(i, indent=4))

# The data directory has already-processed data, so just need to
# init. "Processing" is just copying to other db.
if __name__ == "__main__":
    init_db_raw()
    #fix_db_raw()
    #fudge_gap_in_hf_o2f()
    init_db_processed()
    #export()
    #find_atomic_range()
