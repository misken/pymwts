'''
Created on Jan 13, 2012

@author: mark

This script solves both phases of the mwts model. The solution from phase 1 is fed
to the phase 2 model to create tours.
'''

#Python imports
import sys

from coopr.pyomo import *
import coopr.opt
from pyutilib.misc import import_file

from mwts_phase_1_2_integration import phase_1_2_integrate


# Possible input parameters
isTest = False
if isTest:
    phase1_mod_file = '../tests/prod.py'
    phase1_dat_file = '../tests/prod.dat'
else:
    phase1_mod_file = 'mwts_phase1.py'
    phase2_mod_file = 'mwts_phase2_part1.py'
    phase1_dat_file = '../tests/simple.dat'

which_solver = 'cbc'  # glpk, cbc, gurobi
timelimit = 300
mipgap = 0.10
bWriteStartWinDebug = True


# Setup the solver
solver = None
if which_solver == 'cbc':
    solver = coopr.opt.SolverFactory("cbc")
    solver.options.seconds = timelimit
    solver.options.ratioGap = mipgap
    solver.options.solution = '../tests/solution.sol'
    
if which_solver == 'glpk':
    solver = coopr.opt.SolverFactory("glpk")
    solver.options.tmlim = timelimit
    solver.options.mipgap = mipgap

if which_solver == 'gurobi':
    solver = coopr.opt.SolverFactory("gurobi")
    solver.options.timelimit = timelimit
    solver.options.mipgap = mipgap

# Import phase 1 model and create phase 1 model instance
phase1_mdl = import_file(phase1_mod_file).model_phase1
phase1_inst = phase1_mdl.create(filename = phase1_dat_file)
phase1_inst.name = 'mwts_phase1_inst'

if bWriteStartWinDebug:
    old_stdout = sys.stdout
    sys.stdout = open("../tests/debugwin.txt","w")
    
    for w in phase1_inst.WEEKS:
        for j in phase1_inst.DAYS:
            for i in phase1_inst.PERIODS:
                print 'b_window_epoch[{0},{1},{2}] = {3}, e_window_epoch[{0},{1},{2}] = {4}'.format(i,j,w,phase1_inst.b_window_epoch[i,j,w].value, phase1_inst.e_window_epoch[i,j,w].value)
            
    for (i,j,w) in phase1_inst.PotentialStartWindow_index:     
        if phase1_inst.PotentialStartWindow[i,j,w]:
            print 'PotentialStartWindows[{},{},{}] = {}'.format(i, j, w, phase1_inst.PotentialStartWindow[i,j,w].value)
             
    for (i,j,w,k,t) in phase1_inst.okStartWindow_index:     
        if phase1_inst.okStartWindow[i,j,w,k,t]:
            print 'okStartWindows[{},{},{},{},{}] = {}'.format(i, j, w, k, t, phase1_inst.okStartWindow[i,j,w,k,t].value)
    
    print 'okStartWindowBeginnings_index = ' 
    for (t,k) in phase1_inst.okStartWindowBeginnings_index:
        print t, k
    
      
    for (t,k) in phase1_inst.okStartWindowBeginnings_index:
        print 'okStartWindowBeginnings[{},{}] = '.format(t, k)
        for (i, j, w) in phase1_inst.okStartWindowBeginnings[t, k]:
            print i,  j, w
    
    print 'okStartWindowBeginnings, bchain, and echain = '
    for (t,k) in phase1_inst.okStartWindowBeginnings_index:
        print '[{},{}] = '.format(t, k)
        for (i, j, w) in phase1_inst.bchain[t, k]:
            out = '({},{},{})'.format(i,j,w)
            for (x, y, z) in phase1_inst.echain[t, k, i, j, w]:
                out = out + '({},{},{}) {}\n'.format(x,y,z,phase1_inst.n_links[t, k, i, j, w].value)
            print out
    
    for (t,k,i,j,w) in phase1_inst.chain_index:
        out = ''
        print 'chain[{},{},{},{},{}]='.format(t,k,i,j,w)
        for (x, y, z) in phase1_inst.chain[t, k, i, j, w]:
            out = out + '({},{},{})*'.format(x,y,z,phase1_inst.chain[t, k, i, j, w].value)
        print out
        
    for (t,k,i,j,w,m) in phase1_inst.link_index:
        print 'link[{},{},{},{},{},{}]='.format(t,k,i,j,w,m)
        
    
    #chain[1,1,37,1,3]=
    #(41,1,3)*(39,1,3)*(37,1,3)*(43,1,3)*    
    
    #for m in range(1,phase1_inst.n_links[1,1,37,1,3]+1):
    #    print 'link[{},{},{},{},{},{}] = {}'.format(1,1,37,1,3,m,phase1_inst.link[1,1,37,1,3,m].value)
    #    out = ''
    #    for (x, y, z) in phase1_inst.linkspan[1,1,37,1,3, m]:
    #        out = out + '({},{},{})*'.format(x,y,z,phase1_inst.linkspan[1,1,37,1,3, m].value)
    #    print out
    
    
    sys.stdout = old_stdout


# Optimize phase 1
if solver is None:
    print "Could not get solver: " + which_solver
    sys.exit(1)
else:

    phase1_results = solver.solve(phase1_inst)
    phase1_inst.load(phase1_results)  # Put results in model instance
    
    old_stdout = sys.stdout
    sys.stdout = open("../tests/mwts_phase1_inst.txt","w")
    phase1_inst.pprint()
    sys.stdout = old_stdout
    

    
    print phase1_results['Solution'][0]['Objective'][1]['Value']
    
    if isTest:
        print str(value(phase1_inst.Total_Profit))
        for i in phase1_inst.P:
            print i
            msg = str(phase1_inst.X[i]) + " = " + str(phase1_inst.X[i].value)
            print msg
    else:
        old_stdout = sys.stdout
        sys.stdout = open("../tests/mwts_phase1_out.txt","w")
        
        phase1_solution_status = phase1_results['Solution'][0]['Status']
        print phase1_solution_status
        print str(value(phase1_inst.total_cost))
        tot_cons = 0
        tot_vars = 0
        print "\n\nConstraint summary \n------------------"
        for con_name in phase1_inst.active_components(Constraint):
            con = getattr(phase1_inst,con_name)
            print con_name + " --> " + str(len(con))
            tot_cons += len(con)
        print "\n\nVariable summary \n------------------"
        for var_name in phase1_inst.active_components(Var):
            var = getattr(phase1_inst,var_name)
            print var_name + " --> " + str(len(var))
            tot_vars += len(var)
            
        msg = "\ntotal cons = " + str(tot_cons)
        print msg
        msg = "total vars = " + str(tot_vars)
        print msg
        
        
        phase1_results.write()
        sys.stdout = old_stdout
        
        # If phase 1 solved, create phase 2 params from phase 1 vars and solve phase 2 problem. 
        # TODO - how to check for status other than 'optimal'?
         
        phase2_mdl = import_file(phase2_mod_file).model_phase2
        phase2_inst = phase2_mdl.create(filename = phase1_dat_file)
        phase2_inst.name='mwts_phase2_inst'
        
        # Need to convert the phase2 instance to concrete mode before we can add constraints, fix values, etc.
        phase2_inst.concrete_mode()
        
        # In general, for concrete models, instance and model object share the same memory space (unlike
        # with an abstract model. Not sure if putting an abstract model into concrete_mode does the same thing.
        # We will reference the instance object to be safe.
        
        phase_1_2_integrate(phase1_inst, phase2_inst)

#        n_tours = sum(phase1_inst.TourType[i,t].value for (i,t) in phase1_inst.okTourType)
#        
#        phase2_inst.n_tours = Param(initialize=n_tours)
#        
#        phase2_inst.TOURS = RangeSet(1,phase2_inst.n_tours)
#      
#        for (i,j,w,k,t) in phase1_inst.okShifts:
#            phase2_inst.Shift[i,j,w,k,t].value = phase1_inst.Shift[i,j,w,k,t].value
#
#        for (i,t) in phase1_inst.TourType_index:
#            phase2_inst.TourType[i,t].value = phase1_inst.TourType[i,t].value
#            print i, t, phase2_inst.TourType[i,t].value
#            
#        for (i,t,j,w) in phase1_inst.DailyTourType_index:
#            phase2_inst.DailyTourType[i,t,j,w].value = phase1_inst.DailyTourType[i,t,j,w].value
#                                                                         
#        for (i,t,k,j,w) in phase1_inst.ok_daily_shift_index:
#            phase2_inst.DailyShiftWorked[i,t,k,j,w].value = phase1_inst.DailyShiftWorked[i,t,k,j,w].value
#                                                                         
#        for (i,t,p) in phase1_inst.ok_weekenddaysworked_index:
#            phase2_inst.WeekendDaysWorked[i,t,p].value = phase1_inst.WeekendDaysWorked[i,t,p].value                                                                 
#      
#        # Initialize the phase 2 instance
#        
        print "n_tours = " + str(phase2_inst.n_tours.value)
        
phase2_inst.WIN_x = Param(phase2_inst.TOURS,default=0)
phase2_inst.TT_x = Param(phase2_inst.TOURS,default=0)
        
def index_tours_init(M):

    tour_num = 0
    for (i,t) in M.okTourType:
        if M.TourType[i,t].value > 0:
            print 'TourType[{},{}] = {}'.format(str(i),str(t),str(M.TourType[i,t].value))
            tnum_lower = tour_num + 1
            tnum_upper = tour_num + int(M.TourType[i,t].value)
            print tnum_lower, tnum_upper
   
            tour_num = tour_num + int(M.TourType[i,t].value)

            for idx in range(tnum_lower, tnum_upper + 1):
                M.WIN_x[idx].value = i
                M.TT_x[idx].value = t
                print idx, M.WIN_x[idx].value, M.TT_x[idx].value
                print tnum_lower, 1, M.WIN_x[1].value, M.TT_x[1].value
        
    for i in M.TOURS:
        print i, M.WIN_x[i].value, M.TT_x[i].value



index_tours_init(phase2_inst)

for i in phase2_inst.TOURS:
    print i, phase2_inst.WIN_x[i].value, phase2_inst.TT_x[i].value

phase2_inst.preprocess()


        
        
        
        





#    sys.stdout = old_stdout
#    print results['Solution'][0]['Objectivepython int'][1]['Value']
##
## Write the output
#
##for (i,j,w) in instance.bins:
#
#sys.stdout = old_stdout
## Update the results, to use the same labels as the model
#    transformed_results = instance.update_results(results)
#    sol = transformed_results.solution[0]
#    for var in sorted(sol.variable.keys()):
#            if 'Shift' in var or 'Tour' in var or 'Weekend' in var:
#                #print "  Variable",var
#                for key in sorted(sol.variable[var].keys()):
#                    if key == 'Value' and sol.variable[var]['Value'] > 0.0:
#                        print var ,key, sol.variable[var]['Value']

#
##print yout['Solution'][1]['Variable']['Shift(1,1,1,1,3)']['Value']
