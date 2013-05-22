# Imports
from coopr.pyomo import *
	
# Create the model object
model = AbstractModel()
	
# Sets
model.P = Set()
	
# Parameters
model.a = Param(model.P)
model.b = Param()
model.c = Param(model.P)
model.u = Param(model.P)
	

