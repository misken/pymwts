=======================================
Overview of pymts src and related files
=======================================

To Do
=====

Currently it seems that I need to run the solvemwts.py from scripts that live in 
main src directory for the models. 

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


Virtual Environment for Coopr
=============================

If I type `coopr`, a virtual environment is activated. I don't recall how I did this. 
Obviously it uses virtualenv. Found my blogger post:

While mucking around with the project structure for pymwts, I ended up reinstalling Coopr. I used method 1 which installs Coopr into dist-utils. Chaos ensued with all kinds of bizarre errors. Coopr forum hunting led to the suspicion that some of the errors were bugs that had been fixed in newer versions of Coopr. However, the only way to get these newer versions use to use the 2nd install method:

Option 2: Install a Coopr Release with the coopr_install Scrip

The nice thing about this method is that it installs a virtual Python environment. Why is this nice?

Gets rid of the annoying errors that you get when Pyomo tries to write to some system folders.
Gets rid of the crazy errors mentioned above.
Forces me to learn to use virtualenv which seems to be widely used.
Allows me to get the latest and greatest versions of Coopr.

Activating a virtual Python environment is easy. I installed Coopr in ~/Tools/coopr. To activate it:

source ~/Tools/coopr/bin/activate

I created an alias named coopr to run the above command. To deactivate the virtualenv you just type deactivate.

Running test problems
=====================

Method 1 - calling solvemwts from a Python script
-------------------------------------------------

Here's a simple script called basictest.py::

    """
    basictest.py - Created on Wed May 21 11:12:49 2013

    @author: mark
    """
    import pymwts.solvemwts as solve

    solve.solvemwts('mwts02_d01_t01_tight',
                        'mwts02_d01_t01_tight.dat',
                        './',
                        'cbc',1200,0.05)

You can run it like any Python script::

    (coopr)$ python basictest.py


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





