# Imports
from coopr.pyomo import *
from pyutilib.misc import import_file
	
# Create the prod1_mdl object
prod1_mdl = import_file("prod_1.py").model
prod1_mdl.name = "prod1"

prod1_inst = prod1_mdl.create(filename = 'prod.dat')
	
solver = coopr.opt.SolverFactory("cbc")
solver.options.seconds = 30

print "\nSolve prod_1 and display results"
prod1_results = solver.solve(prod1_inst)
prod1_inst.load(prod1_results)  # Put results in model instance

print str(value(prod1_inst.Total_Profit))
for i in prod1_inst.P:
   print i
   msg = str(prod1_inst.X[i]) + " = " + str(prod1_inst.X[i].value)
   print msg

# Create the prod2_mdl object
prod2_mdl = import_file("prod_2.py").model
prod2_mdl.name = "prod2"

prod2_inst = prod2_mdl.create(filename = 'prod.dat')
prod2_inst.pprint()

print "\nTrying to set X params in prod_2 from solution of prod_1"
for i in prod2_inst.P:
   print i, str(prod2_inst.X[i].value)
   prod2_inst.X[i].value = prod1_inst.X[i].value + 100
   print i, str(prod2_inst.X[i].value)

prod2_inst.pprint()

print "\nSolving prod_2 and trying to display results"
prod2_results = solver.solve(prod2_inst)
prod2_inst.load(prod2_results)  # Put results in model instance

print str(value(prod2_inst.Total_Profit))
for i in prod2_inst.P:
   print i
   msg = str(prod2_inst.X[i].name) + " = " + str(prod2_inst.X[i].value)
   print msg
