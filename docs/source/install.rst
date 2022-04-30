# Installation and dependencies

We assume you
are already comfortable with using Python and familiar with optimization solvers (e.g. CBC, glpk, or
Gurobi), ``git``, Github, installing Python programs using ``pip``,
Python virtual environments, and running programs from a command shell.

-  Need to have either CBC, glpk or Gurobi installed and available to
   use as the mixed-integer programming solver
-  Clone or download the source code from
   `https://github.com/misken/pymwts <https://github.com/misken/pymwts>`_
-  Open a command shell in the main project directory ``pymtws/``.

It is recommended to create a virtual environment within which to
install ``pymwts`` to avoid adding such tools to your base Python
environment. Then just use ``pip`` to install it and navigate to the
``examples/`` subfolder after installation is complete. The ``pymwts``
package depends on a few other Python packages, namely,
`pandas <https://pandas.pydata.org/>`__ and
`pyomo <http://www.pyomo.org/>`__. Both of these will get installed
automatically if they aren’t already installed.

-  ``pip install .``
-  ``cd examples``

After installing pymwts, you can run it from a command shell. Let’s
run it with the ``-h`` flag to see the help info.

    $ pymtws -h

    usage: pymwts [-h] [-p PATH] [-s {cbc,glpk,gurobi}] [-t TIMELIMIT] [-g MIPGAP]
                  [--version]
                  scenario phase1dat

    Solve a multi-week tour scheduling problem.

    positional arguments:
      scenario              Short string to be used in output filenames
      phase1dat             DAT file for phase 1

    optional arguments:
      -h, --help            show this help message and exit
      -p PATH, --path PATH  Relative path to output file directory. Terminate with
                            /
      -s {cbc,glpk,gurobi}, --solver {cbc,glpk,gurobi}
                            cbc, glpk or gurobi for now
      -t TIMELIMIT, --timelimit TIMELIMIT
                            seconds
      -g MIPGAP, --mipgap MIPGAP
                            Can prevent really long run times.
      --version             show program's version number and exit
      
An example Jupyter notebook and associated data files are in the ``examples/`` folder and
the ``examples/inputs/`` subfolder, respectively. The notebook shows how pymwts might
be used in a simple scheduling analysis problem.
