"""
Solve both phases of the mwts model. 

The solution from phase 1 is fed
to the phase 2 model to create tours. Phase 1 solution values for variables are converted
to GMPL DAT parameter statements and a phase 2 DAT file is created from a copy of the
phase 1 DAT file and these additional parameters.
"""

#Python imports
import sys
import os
import time
import sqlite3
import datetime

import pyomo.opt
from pyomo.environ import *

# Need to edit the __init__.py for the package to do the detailed imports
#from .mwts_phase1 import model_phase1
#from .mwts_phase2 import model_phase2

from pymwts import mwts_phase1
from pymwts import mwts_phase2

from pymwts.mwts_utils import *
from pymwts.pymwtsio.mwts_process_out_tour import create_mwt

# Possible input parameters

# scenario - Required 
# phase 1 dat file - Required


# solver Required
# timelimit Required
# mipgap Required
# bWriteStartWinDebug Optional default=false
# phase 1 model Optional default='mwts_phase1.py'
# phase 2 model Optional default='mwts_phase2.py'
# input YAML config file Optional

def solvemwts(scenario,phase1_dat_file,path,
              which_solver,timelimit,mipgap,
              phase1_mod_file='mwts_phase1.py', 
              phase2_mod_file='mwts_phase2.py',
              results_db='exps/mwts00/mwts00.db',
              bWriteStartWinDebug=False,
              bWritePhase1Instance=True,
              bWritePhase2Instance=False):
    
    start_win_debug_file = path + scenario + '_debugwin.txt'
    phase2_dat_file = path + scenario + '_phase2.dat'
    
    phase1_inst_file = path + scenario + '_phase1_inst.txt'
    phase2_inst_file = path + scenario + '_phase2_inst.txt'
    
    log_file = path + scenario + '.log'
    f_log = open(log_file,"w")
    logger(f_log,'Scenario ' + scenario,time.clock())
    logger(f_log,'DAT ' + phase1_dat_file,time.clock())
    
    
    phase1_summary_file = path + scenario + '_phase1_summary.txt'
    phase2_summary_file = path + scenario + '_phase2_summary.txt'
    
    phase1_results_file = path + scenario + '_phase1_results.yml'
    phase2_results_file = path + scenario + '_phase2_results.yml'
    
    phase1_shiftsum_file = path + scenario + '_phase1_shiftsum.csv'
    phase1_tourskeleton_file = path + scenario + '_phase1_tourskeleton.csv'
    
    tour_file = path + scenario + '.tur'
    
    phase1_solution_status = 'untried'
    phase2_solution_status = 'untried'
    phase1_solution_value = 0.0
    phase2_solution_value = 0.0
    
    

    
    # Import phase 1 model and create phase 1 model instance
    #phase1_mdl = import_file(phase1_mod_file).model_phase1
    phase1_mdl = mwts_phase1.model_phase1
    # Trying to get all this working as a package of modules so will try to avoid
    # using import_file (since it's really just emulating a module import anyway - pyutilib.misc)
    #phase1_mdl = mwts_phase1.model_phase1
    phase1_inst = phase1_mdl.create_instance(filename = phase1_dat_file)
    phase1_inst.name = 'mwts_phase1_inst'
    logger(f_log,'Phase 1 instance created',time.clock())
    
    


# Activate/deactivate constraints
    
    b_weekend_subsets_5_4_con_active = True
    b_weekend_subsets_5_5_con_active = False # Feels redundant, see comments in constraint
    b_weekend_subsets_5_5lb_con_active = False
    b_weekend_subsets_5_5sun_con_active = False
    b_weekend_subsets_5_5sat_con_active = False
    b_weekend_subsets_5_4sun_con_active = False
    b_weekend_subsets_5_4sat_con_active = False

    b_weekend_subsets_4_3_con_active = True
    b_weekend_subsets_4_4_con_active = False

    b_weekend_subsets_3_2_con_active = True

    b_weekend_subsets_2_1_con_active = True
    
    b_DTT_TT_fullwkendadj_UB_active = True
    
    b_ad_hoc_weekend_subsets_ttype7_active = True

   
    if not b_weekend_subsets_5_4_con_active:
        phase1_inst.weekend_subsets_5_4_con.deactivate()
        
    if not b_weekend_subsets_5_5_con_active:
        phase1_inst.weekend_subsets_5_5_con.deactivate()
        
    if not b_weekend_subsets_5_5lb_con_active:
        phase1_inst.weekend_subsets_5_5lb_con.deactivate()
        
    if not b_weekend_subsets_5_5sun_con_active:
        phase1_inst.weekend_subsets_5_5sun_con.deactivate()
        
    if not b_weekend_subsets_5_5sat_con_active:
        phase1_inst.weekend_subsets_5_5sat_con.deactivate()
        
    if not b_weekend_subsets_5_4sun_con_active:
        phase1_inst.weekend_subsets_5_4sun_con.deactivate()
        
    if not b_weekend_subsets_5_4sat_con_active:
        phase1_inst.weekend_subsets_5_4sat_con.deactivate()

    if not b_weekend_subsets_4_3_con_active:
        phase1_inst.weekend_subsets_4_3_con.deactivate()

    if not b_weekend_subsets_4_4_con_active:
        phase1_inst.weekend_subsets_4_4_con.deactivate()

    if not b_weekend_subsets_3_2_con_active:
        phase1_inst.weekend_subsets_3_2_con.deactivate()

    if not b_weekend_subsets_2_1_con_active:
        phase1_inst.weekend_subsets_2_1_con.deactivate()

    if not b_DTT_TT_fullwkendadj_UB_active:
        phase1_inst.DTT_TT_fullwkendadj_UB.deactivate()   
        
    if not b_ad_hoc_weekend_subsets_ttype7_active:
        phase1_inst.ad_hoc_weekend_subsets_ttype7.deactivate()  
    
    # Optionally write out out phase 1 instance    
    if bWritePhase1Instance:
        try:          
            f1_inst = open(phase1_inst_file,'w')
            old_stdout = sys.stdout
            sys.stdout = f1_inst
            phase1_inst.pprint()
            logger(f_log,'Phase 1 instance written',time.clock())

        finally:
            sys.stdout = old_stdout
            f1_inst.close()
        
    try:
        f1_sum = open(phase1_summary_file,'w')
        old_stdout = sys.stdout
        sys.stdout = f1_sum
        
            
        tot_cons = 0
        tot_vars = 0
        print("\n\nConstraint summary \n------------------")
        for c in phase1_inst.component_objects(Constraint, active=True):
            #conobj = getattr(phase1_inst,str(c))
            print(c.name + " --> " + str(len(c)))
            tot_cons += len(c)
        print("\n\nVariable summary \n------------------")
        for v in phase1_inst.component_objects(Var):
            #vobj = getattr(phase1_inst,str(v))
            print(v.name + " --> " + str(len(v)))
            tot_vars += len(v)
                
        msg = "\ntotal cons = " + str(tot_cons)
        print(msg)
        msg = "total vars = " + str(tot_vars)
        print(msg)
        f1_sum.close()
        
    finally:
        sys.stdout = old_stdout
    
    if bWriteStartWinDebug:
        old_stdout = sys.stdout
        sys.stdout = open(start_win_debug_file,"w")
        
        for w in phase1_inst.WEEKS:
            for j in phase1_inst.DAYS:
                for i in phase1_inst.PERIODS:
                    print('b_window_epoch[{0},{1},{2}] = {3}, e_window_epoch[{0},{1},{2}] = {4}'.format(i,j,w,phase1_inst.b_window_epoch[i,j,w].value, phase1_inst.e_window_epoch[i,j,w].value))
                
        for (i,j,w) in phase1_inst.PotentialGlobalStartWindow_index:     
            if phase1_inst.PotentialGlobalStartWindow[i,j,w]:
                print('PotentialGlobalStartWindow[{},{},{}] = {}'.format(i, j, w, phase1_inst.PotentialGlobalStartWindow[i,j,w].value))
                 
        for (i,j,w,k,t) in phase1_inst.PotentialStartWindow_index:     
            if phase1_inst.PotentialStartWindow[i,j,w,k,t]:
                print('PotentialStartWindow[{},{},{},{},{}] = {}'.format(i, j, w, k, t, phase1_inst.PotentialStartWindow[i,j,w,k,t].value))
        
        print('okStartWindowRoots_index = ') 
        for (t,k) in phase1_inst.okStartWindowRoots_index:
            print(t, k)
        
          
        for (t,k) in phase1_inst.okStartWindowRoots_index:
            print('okStartWindowRoots[{},{}] = '.format(t, k))
            for (i, j, w) in phase1_inst.okStartWindowRoots[t, k]:
                print(i,  j, w)
                
        print('okTourType = ') 
        for (i,t) in phase1_inst.okTourType:
            print(i, t)
            
        print('okDailyTourType = ') 
        for (i,t,d) in phase1_inst.okDailyTourType:
            print(i, t, d)
        
        print('bchain echain #links = ')
        for (t,k) in phase1_inst.okStartWindowRoots_index:
            print('(t,k) = [{},{}] = '.format(t, k))
            for (i, j, w) in phase1_inst.bchain[t, k]:
                out = '({},{},{})'.format(i,j,w)
                for (x, y, z) in phase1_inst.echain[t, k, i, j, w]:
                    out = out + '({},{},{}) {}\n'.format(x,y,z,phase1_inst.n_links[t, k, i, j, w].value)
                print(out)
        
        for (t,k,i,j,w) in phase1_inst.chain_index:
            out = ''
            print('chain[{},{},{},{},{}]='.format(t,k,i,j,w))
            for (x, y, z) in phase1_inst.chain[t, k, i, j, w]:
                out = out + '({},{},{})*'.format(x,y,z,phase1_inst.chain[t, k, i, j, w].value)
            print(out)
            
        for (t,k,i,j,w,m) in phase1_inst.link_index:
            print('link[{},{},{},{},{},{}]='.format(t,k,i,j,w,m))
            
        
        #chain[1,1,37,1,3]=
        #(41,1,3)*(39,1,3)*(37,1,3)*(43,1,3)*    
        
        #for m in range(1,phase1_inst.n_links[1,1,37,1,3]+1):
        #    print 'link[{},{},{},{},{},{}] = {}'.format(1,1,37,1,3,m,phase1_inst.link[1,1,37,1,3,m].value)
        #    out = ''
        #    for (x, y, z) in phase1_inst.linkspan[1,1,37,1,3, m]:
        #        out = out + '({},{},{})*'.format(x,y,z,phase1_inst.linkspan[1,1,37,1,3, m].value)
        #    print out
        
        
        sys.stdout = old_stdout
        logger(f_log,'Windows debug info written',time.clock())

    # solver = pyomo.opt.SolverFactory('cplex')
    # results = solver.solve(self.m, tee=True, keepfiles=False,
    #                        options_string="mip_tolerances_integrality=1e-9 mip_tolerances_mipgap=0")
    #

    # TODO - add status messages during overall model generation and solution process

    # Setup the solver
    solver = None
    if which_solver == 'cbc':
        solver = pyomo.opt.SolverFactory('cbc')
        solver.options.seconds = timelimit
        solver.options.ratioGap = mipgap
    # solver.options.solution = '../tests/solution.sol' This isn't the correct way to specify this option

    if which_solver == 'glpk':
        solver = pyomo.opt.SolverFactory('glpk')
        solver.options.tmlim = timelimit
        solver.options.mipgap = mipgap

    if which_solver == 'gurobi':
        solver = pyomo.opt.SolverFactory('gurobi')
        solver.options.timelimit = timelimit
        solver.options.mipgap = mipgap

    # Optimize phase 1
    if solver is None:
        print("Could not get solver: " + which_solver)
        sys.exit(1)
    else:
    
        stream_solver = True
        phase1_results = solver.solve(phase1_inst, tee=stream_solver)
        
        # TODO - check if phase 1 solved        
        if (phase1_results.solver.status != pyomo.opt.SolverStatus.ok):
            #logging.warning('Check solver not ok?')
            logger(f_log, 'Check solver not ok?', phase1_results.solver.status)
        if (phase1_results.solver.termination_condition != pyomo.opt.TerminationCondition.optimal):
            #logging.warning('Check solver optimality?')
            logger(f_log, 'Check solver optimality?', phase1_results.solver.termination_condition)
        
        logger(f_log,'Phase 1 solution status=' + str(phase1_results.Solution.Status),time.clock())
        
        logger(f_log,'Phase 1 solved',time.clock())

        phase1_inst.solutions.load_from(phase1_results)  # Put results in model instance

        logger(f_log,'Phase 1 results loaded',time.clock())

        
        
        try:
            #instance.solutions.load_from(results)
            phase1_solution_value = phase1_results['Solution'][0]['Objective'][1]['Value']
            print(phase1_solution_value)
        except:
            print('Phase 1 problem not solved successfully.')
            phase1_solution_status = phase1_results['Solution'][0]['Status']
            print('Status: ' + str(phase1_solution_status))
            sys.exit(1)
     
        
        old_stdout = sys.stdout
        try:
            f1_sum = open(phase1_summary_file,'a')
            sys.stdout = f1_sum
        
            
            tot_cons = 0
            tot_vars = 0
            print("\n\nConstraint summary \n------------------")
            for c in phase1_inst.component_objects(Constraint, active=True):
                #conobj = getattr(phase1_inst, str(c))
                print(c.name + " --> " + str(len(c)))
                tot_cons += len(c)
            print("\n\nVariable summary \n------------------")
            for v in phase1_inst.component_objects(Var):
                #vobj = getattr(phase1_inst, str(v))
                print(v.name + " --> " + str(len(v)))
                tot_vars += len(v)
                
            msg = "\ntotal cons = " + str(tot_cons)
            print(msg)
            msg = "total vars = " + str(tot_vars)
            print(msg)
            
            phase1_solution_status = phase1_results['Solution'][0]['Status']
            print(phase1_solution_status)
            
            
            f1_sum.close()
            f1_res = open(phase1_results_file,"w")
            sys.stdout = f1_res
            phase1_results.write()
            f1_res.close()
            logger(f_log,'Phase 1 summary and results written',time.clock())
        
        finally:
            sys.stdout = old_stdout
        
        # Write shift summary
        phase1_shiftsummary = write_phase1_shiftsummary(phase1_inst)
        f1_shiftsum = open(phase1_shiftsum_file,'w')
        print(phase1_shiftsummary, file=f1_shiftsum)
        f1_shiftsum.close()
        
        # Write tour skeleton
        phase1_tourskeleton = weekenddaysworked_to_tourskeleton(phase1_inst)
        f1_tourskeleton = open(phase1_tourskeleton_file,'w')
        print(phase1_tourskeleton, file=f1_tourskeleton)
        phase1_tourskeleton = dailytourtype_to_tourskeleton(phase1_inst)
        print(phase1_tourskeleton, file=f1_tourskeleton)
        f1_tourskeleton.close()
        
        
        
        # If phase 1 solved, create phase 2 params from phase 1 vars and solve phase 2 problem. 
        # TODO - how to check for status other than 'optimal'?
         
        
        
        # Need to convert the phase2 instance to concrete mode before we can add constraints, fix values, etc.
        # phase2_inst.concrete_mode()
        
        # In general, for concrete models, instance and model object share the same memory space (unlike
        # with an abstract model. Not sure if putting an abstract model into concrete_mode does the same thing.
        # We will reference the instance object to be safe.
        
        #phase_1_2_integrate(phase1_inst, phase2_inst)
        
        tot_cap = value(sum(phase1_inst.cov[i,j,w].value for (i,j,w) in phase1_inst.bins))
        us1_cost = value(sum(phase1_inst.under1[i,j,w] * phase1_inst.cu1.value for (i,j,w) in phase1_inst.bins))
        us2_cost = value(sum(phase1_inst.under2[i,j,w] * phase1_inst.cu2.value for (i,j,w) in phase1_inst.bins))
    
        n_tours = int(round(sum(phase1_inst.TourType[i,t].value for (i,t) in phase1_inst.okTourType)))
        
        param_n_tours = scalar_to_param('n_tours',n_tours)
    
        param_Shift = shift_to_param('Shift',phase1_inst)
        param_TourType = tourtype_to_param('TourType',phase1_inst)
        param_DailyTourType = dailytourtype_to_param('DailyTourType',phase1_inst)
        param_DailyShiftWorked = dailyshiftworked_to_param('DailyShiftWorked',phase1_inst)
        param_WeekendDaysWorked = weekenddaysworked_to_param('WeekendDaysWorked',phase1_inst)
        
        param_tour_WIN_TT = tour_WIN_TT_to_param(phase1_inst)
        
        f2_out = open(phase2_dat_file,'w')
        f1_in = open(phase1_dat_file,'r')
        f1_str = f1_in.read()
        f2_out.write(f1_str)
        f1_in.close()
        f2_out.close()
        
        # Put the pieces together
        dat = StringIO.StringIO()
    
        print(param_n_tours, file=dat)
        print(param_Shift, file=dat)  
        print(param_TourType, file=dat)
        print(param_DailyTourType, file=dat)
        print(param_DailyShiftWorked, file=dat)
        print(param_WeekendDaysWorked, file=dat)
        print(param_tour_WIN_TT, file=dat)
    
        f2_dat = open(phase2_dat_file,'a')
        print(dat.getvalue(), file=f2_dat)
        f2_dat.close()
        logger(f_log,'Phase 2 dat file created',time.clock())

        
        
        # Initialize the phase 2 instance
        
        #phase2_mdl = import_file(phase2_mod_file).model_phase2
        phase2_mdl = mwts_phase2.model_phase2
        phase2_inst = phase2_mdl.create_instance(filename = phase2_dat_file)
        phase2_inst.name = 'mwts_phase2_inst'
        logger(f_log,'Phase 2 instance created',time.clock())
        
        # Activate/deactivate constraints
        
        bOnePatternPerTourShift_active = True
        bTours_Daily_active = True
        bTours_Daily_conservation_active = True
        
        bTour_WkendDof_conservation_active = True
        bTour_ShiftWkendDof_integration1_active = True
        
        bTours_Weekly_LB_active = True
        bTours_Weekly_UB_active = True
        bTours_Total_LB_active = True
        bTours_Total_UB_active = True

        bTours_Shiftlen_Weekly_LB_active = True
        bTours_Shiftlen_Weekly_UB_active = True
        bTours_Shiftlen_Total_LB_active = True
        bTours_Shiftlen_Total_UB_active = True 
            
        bTours_Weekly_Prds_LB_active = True
        bTours_Weekly_Prds_UB_active = True
        bTours_Total_Prds_LB_active = True
        bTours_Total_Prds_UB_active = True
    
        # Deactivate constraints per the above list of binaries
        
        if not bOnePatternPerTourShift_active:
            phase2_inst.OnePatternPerTourShift.deactivate()
         
        if not bTour_WkendDof_conservation_active: 
            phase2_inst.Tour_WkendDof_conservation.deactivate()
    
        if not bTour_ShiftWkendDof_integration1_active:
            phase2_inst.Tour_ShiftWkendDof_integration1.deactivate()
    
        if not bTours_Daily_active:
            phase2_inst.Tours_Daily.deactivate()
            
        if not bTours_Daily_conservation_active:
            phase2_inst.Tours_Daily_conservation.deactivate()
         
        if not bTours_Weekly_LB_active: 
            phase2_inst.Tours_Weekly_LB.deactivate()
            
        if not bTours_Weekly_UB_active: 
            phase2_inst.Tours_Weekly_UB.deactivate()
       
        if not bTours_Total_LB_active:
            phase2_inst.Tours_Total_LB.deactivate()
         
        if not bTours_Total_UB_active:         
            phase2_inst.Tours_Total_UB.deactivate()
                 
        if not bTours_Shiftlen_Weekly_LB_active: 
            phase2_inst.Tours_Shiftlen_Weekly_LB.deactivate()
                    
        if not bTours_Shiftlen_Weekly_UB_active: 
            phase2_inst.Tours_Shiftlen_Weekly_UB.deactivate()
       
        if not bTours_Shiftlen_Total_LB_active:
            phase2_inst.Tours_Shiftlen_Total_LB.deactivate()
         
        if not bTours_Shiftlen_Total_UB_active:         
            phase2_inst.Tours_Shiftlen_Total_UB.deactivate()
            
        if not bTours_Weekly_Prds_LB_active: 
            phase2_inst.Tours_Weekly_Prds_LB.deactivate()
            
        if not bTours_Weekly_Prds_UB_active: 
            phase2_inst.Tours_Weekly_Prds_UB.deactivate()
            
        if not bTours_Total_Prds_LB_active: 
            phase2_inst.Tours_Total_Prds_LB.deactivate()
            
        if not bTours_Total_Prds_UB_active: 
            phase2_inst.Tours_Total_Prds_UB.deactivate()
        
        # Optionally write out phase 2 problem instance
        if bWritePhase2Instance:
            try:            
                f2_inst = open(phase2_inst_file,'w')
                old_stdout = sys.stdout
                sys.stdout = f2_inst
                phase2_inst.pprint()
                logger(f_log,'Phase 2 instance',time.clock())
            finally:
                sys.stdout = old_stdout
                f2_inst.close()
        
        # Write out phase 2 problem size summary info
        old_stdout = sys.stdout
        try:
            f2_sum = open(phase2_summary_file,'w')
            sys.stdout = f2_sum
        
            
            tot_cons = 0
            tot_vars = 0
            print("\n\nConstraint summary \n------------------")
            for c in phase2_inst.component_objects(Constraint, active=True):
                conobj = getattr(phase2_inst, str(c))
                print(c + " --> " + str(len(conobj)))
                tot_cons += len(conobj)
            print("\n\nVariable summary \n------------------")
            for v in phase2_inst.component_objects(Var):
                vobj = getattr(phase2_inst, str(v))
                print(v + " --> " + str(len(vobj)))
                tot_vars += len(vobj)
                
            msg = "\ntotal cons = " + str(tot_cons)
            print(msg)
            msg = "total vars = " + str(tot_vars)
            print(msg)
            f2_sum.close()
        finally:
            sys.stdout = old_stdout
            
        # Solve phase 2
        stream_solver = True
        phase2_results = solver.solve(phase2_inst, tee=stream_solver)
        phase2_solution_status = str(phase2_results.Solution.Status)
        logger(f_log,'Phase 2 solution status=' + str(phase2_results.Solution.Status),time.clock())
        
        if str(phase2_results.Solution.Status) != 'unknown':
            logger(f_log,'Phase 2 solved',time.clock())
            #phase2_inst.solutions.load_from(phase2_results)  # Put results in model instance
            logger(f_log,'Phase 2 results loaded',time.clock())
            
            
                
                
            # print phase2_results['Solution'][0]['Objective'][1]['Value']
            old_stdout = sys.stdout
            try:
                f2_sum = open(phase2_summary_file,'a')
                sys.stdout = f2_sum
            
                
                tot_cons = 0
                tot_vars = 0
                print("\n\nConstraint summary \n------------------")
                for c in phase2_inst.component_objects(Constraint, active=True):
                    conobj = getattr(phase2_inst, str(c))
                    print(c + " --> " + str(len(conobj)))
                    tot_cons += len(conobj)
                print("\n\nVariable summary \n------------------")
                for v in phase2_inst.component_objects(Var):
                    vobj = getattr(phase2_inst, str(v))
                    print(v + " --> " + str(len(vobj)))
                    tot_vars += len(vobj)
                    
                msg = "\ntotal cons = " + str(tot_cons)
                print(msg)
                msg = "total vars = " + str(tot_vars)
                print(msg)
                
                phase2_solution_status = phase2_results['Solution'][0]['Status']
                print(phase2_solution_status)
                print(str(value(phase2_inst.total_num_tours)))
                phase2_solution_value = value(phase2_inst.total_num_tours)
                
                f2_sum.close()
                f2_res = open(phase2_results_file,"w")
                sys.stdout = f2_res
                phase2_results.write()
                f2_res.close()
                logger(f_log,'Phase 2 summary and results written',time.clock())
            
            finally:
                sys.stdout = old_stdout
                
                
            
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
        #        for (i,t,p) in phase1_inst.ok_weekenddayswork,index[6],index[7]ed_index:
        #            phase2_inst.WeekendDaysWorked[i,t,p].value = phase1_inst.WeekendDaysWorked[i,t,p].value                                                                 
        #      
        
            
        
        
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
        
        #printf 'PP4\n';
        #printf {l in 1..n_tours, w in WEEKS, i in PERIODS, j in DAYS, d in DAYS, p in PERIODS, k in LENGTHS, t in TTYPES  : 
        #    k in tt_length_x[TT_x[l]] and t=TT_x[l] and 
        #     (i,j) in okWindowWepochs[WIN_x[l],d,k,TT_x[l]] and p=WIN_x[l] and tourshift[l,i,j,w,k,t,p,d] = 1}: 
        #    '%3i%3i%3i%3i%3i%3i%3i%3i\n',l,t,w,i,j,lengths[k],p,d;
            
         
        #for s in phase2_inst.TOURS:
        #    for w in phase2_inst.WEEKS:
        #        for i in phase2_inst.PERIODS: 
        #            for j in phase2_inst.DAYS:
        #                for d in phase2_inst.DAYS:
        #                    for p in [a for a in phase2_inst.PERIODS if a==phase2_inst.WIN_x[s].value]:
        #                        for k in [a for a in phase2_inst.LENGTHS if a == phase2_inst.length_x[phase2_inst.TT_x[s].value]]:
        #                            for t in [a for a in phase2_inst.TTYPES if a == phase2_inst.TT_x[s].value]:
        
        # Create the multi-week tour file
        indexcopy = []
        for index in phase2_inst.TourShifts:
            indexcopy.append(index)
        indexsorted = sorted(indexcopy)
        
        old_stdout = sys.stdout
        f2_tour = open(tour_file,"w")
        sys.stdout = f2_tour
        print('n_tours', phase2_inst.n_tours.value)
        print('n_weeks', phase2_inst.n_weeks.value)
        print('PP4')        
        #for index in phase2_inst.TourShifts:
        for index in indexsorted:
            if phase2_inst.TourShifts[index].value > 0:
        #        print index, index[0], phase2_inst.TourShifts[index], phase2_inst.TourShifts[index].value
                print(index[0],index[5],index[3],index[1],index[2],phase2_inst.lengths[index[4]],phase2_inst.WIN_x[index[0]])
        
        print('TTS')
        for s in phase2_inst.TOURS:
            print(phase2_inst.TT_x[s])
                
        f2_tour.close()
        sys.stdout = old_stdout 
    
        create_mwt(tour_file, scenario, path)
        logger(f_log,'Tour related output files created',time.clock())
        now = datetime.datetime.now()
        vtuple = (scenario,phase1_solution_value,tot_cap,
              str(phase1_solution_status),str(phase2_solution_status),us1_cost,us2_cost,phase2_solution_value,str(now))        
        logger(f_log,'Solution log record',str(vtuple))
    f_log.close()  

    # Connect to the problem solution log database.
    conn = sqlite3.connect(results_db)
    cur = conn.cursor()
    now = datetime.datetime.now()
    flist = '(problem,MIP_obj,tot_cap,phase1_sol_status,phase2_sol_status,us1_cost,us2_cost,phase2_obj,timestamp)'
    vtuple = (scenario,phase1_solution_value,tot_cap,
              str(phase1_solution_status),str(phase2_solution_status),us1_cost,us2_cost,phase2_solution_value,str(now))
    sql_insert = 'insert into solution_log ' + flist + ' values ' + str(vtuple)
    print(sql_insert)
    cur.execute(sql_insert)
    conn.commit()
    

    # Update the problem list table
    
    cur = conn.cursor()
    conn.row_factory = sqlite3.Row
    sql_update = 'update problem_list set sol_status="' + str(phase1_solution_status) + ':' + str(phase2_solution_status) + '" where problem="' + scenario + '"'
    cur.execute(sql_update)    
    print(sql_update)   
    conn.commit()
    conn.close()      
            
             
      
    

    
    
    
    
    
