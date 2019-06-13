=======================================
Overview of pymts src and related files
=======================================

To Do
=====



Files
=====

mwts.py
-------
Provides a command line interface to solvemwts procedure. User can specify inputs 
such as scenario, phase 1 dat file, output path, solver, mipgap, debug mode,
phase 1 model, and phase 2 model. It appears that there were plans to have
the ability to read a YAML formatted input config file - not implemented yet.




mwts_phase1.py
--------------

Main Phase 1 model written in Pyomo


Getting mwts working in 2019
============================

Okay, I'm going to get this multiweek model running, fixed, published,
and released.

Pyomo and conda
---------------------------

According to http://www.pyomo.org/about, Pyomo is no longer referred
to as Coopr software library. It's an official COIN-OR project.

According to installation instructions, Pyomo can be installed using conda:

https://pyomo.readthedocs.io/en/latest/installation.html#using-conda

    conda install -c conda-forge pyomo
    
    conda install -c conda-forge pyomo.extras
    
Looks like several open source solvers can also be installed via conda:

    conda install -c conda-forge ipopt coincbc glpk
    
    
After doing the above three installs in a conda environment I named mwts, 
I found a link (in the docs) to a gallery of
Jupyter notebooks which demo Pyomo. 

https://github.com/jckantor/ND-Pyomo-Cookbook

I downloaded the notebooks and started Jupyter lab in my mwts conda
environment. Upon opening one of the notebooks, Jupyter wanted to
build some extras and I let it. Then I tried the simple LP production
model with glpk and, voila, it solved.

Retried it with cbc and, voila, it solved even faster. :)

The job shop scheduling notebook is awesome! Uses pandas and makes
gantt charts. Super cool.



Running test problems
=====================

I installed pytest on 2019-03-27

https://docs.pytest.org/en/latest/goodpractices.html#choosing-a-test-layout-import-rules


Put test scripts inside /tests folder. Name them to start with 'test_'
so that they can be discovered by pytest if desired.

To run a specific test, run from main source directory - in this case,
from pymtws/pymwts.

    $ python -m pytest ../tests/test_mwts04_d02_t12_a01_ptub_moderate.py
    
Fixed the imports to get rid of `from pymwts import BLAH` and everything just works.


Running real problems
=====================

To run problems using the development version and without installing, need
to add the source directory to the Python search path. Example:

    import sys
    sys.path.insert(0, "/home/mark/Documents/research/MultiWeek/pymwts/pymwts")
    import solvemwts
    solvemwts.solvemwts("mwts05_d02_t1_a01_noptub_loose","./inputs/dat/mwts05_d02_t1_a01_noptub_loose.dat","./outputs/","gurobi",1200.0,0.02,results_db="mwts05_d2456.db")



TODO - Update this section after figuring out new way to package
and install mwts. The stuff below is old

Method 1 - calling solvemwts from a Python script
-------------------------------------------------

Here's a simple script::

    """
    basictest.py - Created on Wed May 21 11:12:49 2013

    @author: mark
    """
    
    import solvemwts
    solvemwts.solvemwts("mwts04_d02_t1_a01_ptub_moderate",
                    "../tests/inputs/mwts04_d02_t1_a01_ptub_moderate_nous.dat",
                    "../tests/outputs/",
                    "cbc",
                    1200.0,
                    0.02,
                    results_db="../tests/mwts04_d2456_cbc_testing.db")

You can run it like any Python script::

    (mwts)$ python basictest.py


Method 2 - using mwts command line Python script
-------------------------------------------------

I created a setup.py file that uses setuptools to install::

    from setuptools import setup

    setup(name='pymwts',
      version='1.0',
      description='Pyomo based mwts',
      author='Mark Isken',
      author_email='isken@oakland.edu',
      url='http://www.hselab.org/machinery',
      packages=['pymwts', 'pymwts.pymwtsio'],
      entry_points = {
        'console_scripts': [
            'mwts = pymwts.mwts:main']}     
     )

This ends up creating a Python script called **mwts** that gets put in
~/Tools/coopr/bin/. It looks like this: ::

    #!/home/mark/Tools/coopr/bin/python
    # EASY-INSTALL-ENTRY-SCRIPT: 'pymwts==1.0','console_scripts','mwts'
    __requires__ = 'pymwts==1.0'
    import sys
    from pkg_resources import load_entry_point

    if __name__ == '__main__':
        sys.exit(
            load_entry_point('pymwts==1.0', 'console_scripts', 'mwts')()
        )

So, **mwts** can be invoked from the command line. ::

    (coopr)$ mwts -h
    usage: mwts [-h] [--version] [-p PATH] [-s {cbc,glpk}] [-t TIMELIMIT]
                [-g MIPGAP] [-w] [-p1 PHASE1MODEL] [-p2 PHASE2MODEL] [-y YAML]
                scenario phase1dat

    Solve a multi-week tour scheduling problem.

    positional arguments:
      scenario              Short string to be used in output filenames
      phase1dat             DAT file for phase 1

    optional arguments:
      -h, --help            show this help message and exit
      --version             show program's version number and exit
      -p PATH, --path PATH  Relative path to output file directory. Terminate with
                            /
      -s {cbc,glpk}, --solver {cbc,glpk}
                            cbc or glpk for now
      -t TIMELIMIT, --timelimit TIMELIMIT
                            seconds
      -g MIPGAP, --mipGap MIPGAP
                            Can prevent really long run times.
      -w, --windebug        Write out start window debug info.
      -p1 PHASE1MODEL, --phase1model PHASE1MODEL
                            Model for phase 1 problem
      -p2 PHASE2MODEL, --phase2model PHASE2MODEL
                            Model for phase 2 problem
      -y YAML, --yaml YAML  YAML input config filename. NOT IMPLEMENTED.

    May the force be with you.





