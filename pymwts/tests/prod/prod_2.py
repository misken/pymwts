# Imports
from coopr.pyomo import *
from pyutilib.misc import import_file
	
# Create the model object
#model = AbstractModel()
model = import_file('prod_base.py').model
model.name = "prod_2"
	


# We will get values for following param from solution to prod.py model
X_init = {}
X_init['bands']=0
X_init['coils']=0
model.X = Param(model.P,initialize=X_init)
	
# Variables
model.dummy = Var(model.P,within=NonNegativeReals)

	
# Objective
def Objective_rule(model):
    return sum([model.c[j]*model.X[j]-model.dummy[j] for j in model.P]) 
model.Total_Profit = Objective(rule=Objective_rule, sense=maximize)
	

# Dummy Constraint
def dummy_rule(model):
    return summation(model.dummy, denom=model.a) < model.b
model.dummy_con = Constraint(rule=dummy_rule)
