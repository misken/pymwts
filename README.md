
pymwts - A multi-week implicit tour scheduling model
====================================================

This project currently contains:


* Python code for Pyomo based scheduling models,
* Sample input and output data files,
* A Jupyter notebook containing a short demo of installing and using pymwts,
* Links to preprint version of the journal paper presenting this model along with associated computational experiment data.

In the future we will add:

* Python code for data input and data output processing and management,
* Documentation including user manual, tutorials, explanations.

The models are mixed integer programming models, represented with the
Pyomo modeling language (a Python library) and solvable with standard
solvers such as Cbc, GLPK, Gurobi, or CPLEX.


Usage
-----

See [scheduling_analysis_example.ipynb](https://github.com/misken/pymwts/blob/master/examples/scheduling_analysis_example.ipynb) for a short demo of installing and using pymwts. There are also html and pdf versions of this example notebook available in the `examples` folder of this repo.

Research paper and associated data
----------------------------------

Link to paper preprint (COMING SOON)

### Data files

* [prob_soln_summary.csv](https://drive.google.com/file/d/1ZLiPxYShPDksZIa-iAo6zyTWQoyOnbLu/view?usp=sharing) - summary of the 696 test problems
* [dat.zip](https://drive.google.com/file/d/1hRkDcz8eu3w3QVfsVQ5CFudB1MOy4J_g/view?usp=sharing) - DAT formatted input files
* [outputs.zip](https://drive.google.com/file/d/1RplO6_EET3UrnuVLT1LfIw0KCDFPT6OO/view?usp=sharing) - all of the various output files created
