#!/usr/bin/env python

'''
Created on Oct 4, 2012
'''
from ga_optimization_ternary.database import Stats_Database
from ga_optimization_ternary.utils import get_reference_array

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Oct 4, 2012"

if __name__ == "__main__":
    sdb = Stats_Database()
    
    sdb_os = Stats_Database(extension="_OS")
    print get_reference_array()[20]
    print get_reference_array()[10]
    
    # print headers
    params = ["crossover_fnc", "popsize", "selection_overall", "mutation_rate", "elitism_num", "fitness_fnc"]  
    for param in params:
        print param,
    print 'ten', 'all', 'half (OS)', 'all (OS)'
        
    for item in sdb._stats_process.find():
        unique_key = item['unique_key']
        if "eval_fitness_complex_product" in unique_key:
            unique_key = unique_key.replace("eval_fitness_complex_product", "eval_fitness_complex_product_oxide_shield")
        else:
            unique_key = unique_key.replace("eval_fitness_complex", "eval_fitness_complex_oxide_shield")
            unique_key = unique_key.replace("eval_fitness_simple", "eval_fitness_simple_oxide_shield")
        
        item2 = sdb_os._stats_process.find_one({"unique_key":unique_key})
        for param in params:
            print item['parameters'][param],
            
        print item['ten'], item['all'], item2['half'], item2['all']