Code that powered [Performance of genetic algorithms in search for water splitting perovskites](http://dx.doi.org/10.1007/s10853-013-7448-9) by A. Jain et al. Data set may be obtained upon request from dwinston [at] lbl [dot] gov.

[`ga_optimization_ternary/simple_ga.py`](https://github.com/mp-interns/ga_optimization_ternary/blob/master/ga_optimization_ternary/simple_ga.py) is the main script.

For setup, once in a Python virtual environment and in the base directory:
```
$ pip install -r requirements.txt
$ export PYTHONPATH=`pwd`:$PYTHONPATH # or e.g. `add2virtualenv .`
```

To set up the db:
```
$ python db_setup/initialize_db.py
```
