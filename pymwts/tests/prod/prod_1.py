# Imports
from coopr.pyomo import *
from pyutilib.misc import import_file
	
# Create the model object
#model = AbstractModel()
model = import_file('prod_base.py').model
model.name = "prod_1"
	

	
# Variables
model.X = Var(model.P)
	
# Objective
def Objective_rule(model):
    return sum([model.c[j]*model.X[j] for j in model.P])
model.Total_Profit = Objective(rule=Objective_rule, sense=maximize)
	
# Time Constraint
def Time_rule(model):
    return summation(model.X, denom=model.a) < model.b
model.Time = Constraint(rule=Time_rule)

	
# Limit Constraint
def Limit_rule(model, j):
    return (0, model.X[j], model.u[j])
model.Limit = Constraint(model.P, rule=Limit_rule)

