from coopr.pyomo import *
from coopr.opt import SolverFactory

# create a solver object
#opt = SolverFactory('cplexamp', solver_io='nl') # Write probelm using the ASL interface
opt = SolverFactory('cbc') 
#opt = SolverFactory('cplex', solver_io='lp') # Write problem using the LP file interface
#opt = SolverFactory('cplex', solver_io='python') # Write problem using the Cplex Python API
# option to send solver output to screen
stream_solver = True

# example options
# tell cplex to solve the dual using the primal simplex algorithm
# Note: These particular cplex options may only be valid through the ASL interface
#solver_options = ['display=2','predual=1','primalopt=""']
# send example options to the solver object
#opt.set_options(" ".join(solver_options))

# define the model 
# At this point you could import the 
# module containing your model or call
# a function that returns a pyomo model
model = ConcreteModel()
model.S = Set(initialize=[1,2,3])
model.x = Var(model.S)
model.obj = Objective(expr=summation(model.x))
model.con = Constraint(model.S,expr=lambda model,i: model.x[i] > i)

# For ConcreteModel you really only need to call .preprocess()
# If you have an AbstractModel then call .create() with your data file
inst = model.create

# ****
# One other thing to note is that for a ConcreteModel, after calling
# inst = model.create()
# the 'inst' and 'model' python variables point to the same object in memory
# this is not the case for an AbstractModel

# solve the instance
results = opt.solve(inst, tee=stream_solver)

# load the solution results back into the instance
inst.load(results)

print 
print "Obective:", value(inst.obj)
for i in inst.S:
    print "x["+str(i)+"]:", value(inst.x[i])

N = sum(value(inst.x[i]) for i in inst.S)

# You can either start a new model from scratch or add to your old model
# *****
# If you started with an AbstractModel, be sure to call
inst.concrete_mode()
# Before trying to add any elements to it. At this point
# it will behave like a ConcreteModel
# *****
inst.Q = Set(initialize=range(1,int(N)+1))
inst.y = Var(inst.Q)
inst.obj = Objective(expr=summation(inst.y)+summation(inst.x))
inst.con2 = Constraint(model.Q,expr=lambda model,i: model.y[i] > i)
# optionally disable phase1 model components
inst.con.deactivate()
# fix the x variables, they now become constants in the constraint expressions
for i in inst.S:
    inst.x[i].fixed=True

# You must call .preprocess after fixing variables, adding/deactivating constraints, etc. 
inst.preprocess()

# solve the instance
results = opt.solve(inst, tee=stream_solver)

# load the solution results back into the instance
inst.load(results)

print 
print "Obective:", value(inst.obj)
for i in inst.Q:
    print "y["+str(i)+"]:", value(inst.y[i])

