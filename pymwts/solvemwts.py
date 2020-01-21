"""
Solve both phases of the mwts model. 

The solution from phase 1 is fed
to the phase 2 model to create tours. Phase 1 solution values for variables are converted
to GMPL DAT parameter statements and a phase 2 DAT file is created from a copy of the
phase 1 DAT file and these additional parameters.
"""

import sys
import sqlite3
import datetime
import logging

import pyomo.opt
from pyomo.environ import *

import mwts_phase1
import mwts_phase2
from mwts_utils import *
from pymwtsio.mwts_process_out_tour import create_mwt

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

def solvemwts(scenario, phase1_dat_file, path,
              which_solver, timelimit, mipgap,
              phase1_mod_file = 'mwts_phase1.py',
              phase2_mod_file = 'mwts_phase2.py',
              results_db = None,
              bWriteStartWinDebug = False,
              bWritePhase1Instance = False,
              bWritePhase2Instance = False,
              force_solve = False):
    """

    :param scenario:
    :param phase1_dat_file:
    :param path:
    :param which_solver:
    :param timelimit:
    :param mipgap:
    :param phase1_mod_file:
    :param phase2_mod_file:
    :param results_db:
    :param bWriteStartWinDebug:
    :param bWritePhase1Instance:
    :param bWritePhase2Instance:
    :return:
    """

    # Check to see if this problem has already been run
    # Connect to the problem solution log database.
    if not force_solve and (results_db is not None):
        conn = sqlite3.connect(results_db)
        cur = conn.cursor()
        row = cur.execute('SELECT Problem, sol_status FROM problem_list WHERE Problem=?', (scenario,))

        result = row.fetchone()
        if result[1] != 'not run':
            conn.close()
            print('Problem {} already run. sol_status = {}'.format(scenario, result[1]))
            exit(0)

        conn.close()

    # Create filenames
    phase2_dat_file = path + scenario + '_phase2.dat'
    
    phase1_inst_file = path + scenario + '_phase1_inst.txt'
    phase2_inst_file = path + scenario + '_phase2_inst.txt'

    phase1_summary_file = path + scenario + '_phase1_summary.txt'
    phase2_summary_file = path + scenario + '_phase2_summary.txt'

    phase1_results_file = path + scenario + '_phase1_results.yml'
    phase2_results_file = path + scenario + '_phase2_results.yml'

    phase1_shiftsum_file = path + scenario + '_phase1_shiftsum.csv'
    phase1_tourskeleton_file = path + scenario + '_phase1_tourskeleton.csv'

    tour_file = path + scenario + '.tur'

    # Setup logging
    log_file = path + scenario + '.log'
    logging.basicConfig(filename=log_file,
                        filemode='w',
                        format='%(levelname)s:%(message)s %(asctime)s',
                        level=logging.DEBUG)

    # Initialization
    phase1_solution_status = 'untried'
    phase2_solution_status = 'untried'
    phase1_solution_value = 0.0
    phase2_solution_value = 0.0

    logging.info('Scenario %s started', scenario)
    logging.info('DAT %s', phase1_dat_file)

    # Create Phase 1 model instance
    phase1_mdl = mwts_phase1.model
    phase1_inst = phase1_mdl.create_instance(filename=phase1_dat_file)
    phase1_inst.name = 'mwts_phase1_inst'
    logging.info('Phase 1 instance created')

    # Activate/deactivate constraints -----------------------------------------

    # Deactivate part-time fraction upper bound if all tour types are part-time
    tot_parttime_ttypes = sum(phase1_inst.tt_parttime[t] for t in phase1_inst.activeTT)
    if tot_parttime_ttypes == len(phase1_inst.activeTT) or phase1_inst.max_parttime_frac.value == 1.0:
        phase1_inst.max_ptfrac_con.deactivate()

    # Boolean indicators for debugging weekends related constraints
    b_weekend_subsets_5_4_con2_active = True
    b_weekend_subsets_4_3_con2_active = True
    b_weekend_subsets_3_2_con2_active = True
    b_weekend_subsets_2_1_con2_active = True

    # Boolean indicators for possible redundant constraints
    b_TTD_TT_weeklyconservation_active = False
    b_TTDS_TT_weeklyconservation_active = True

    # Conditional constraint deactivation
    if not b_weekend_subsets_5_4_con2_active:
        phase1_inst.weekend_subsets_5_4_con2.deactivate()

    if not b_weekend_subsets_4_3_con2_active:
        phase1_inst.weekend_subsets_4_3_con2.deactivate()

    if not b_weekend_subsets_3_2_con2_active:
        phase1_inst.weekend_subsets_3_2_con2.deactivate()

    if not b_weekend_subsets_2_1_con2_active:
        phase1_inst.weekend_subsets_2_1_con2.deactivate()

    # Possible redundant constraint deactivation
    if not b_TTD_TT_weeklyconservation_active:
        phase1_inst.TTD_TT_weeklyconservation_LB.deactivate()
        phase1_inst.TTD_TT_weeklyconservation_UB.deactivate()
        phase1_inst.TTD_TT_cumul_weeklyconservation_LB.deactivate()
        phase1_inst.TTD_TT_cumul_weeklyconservation_UB.deactivate()

    if not b_TTDS_TT_weeklyconservation_active:
        phase1_inst.TTDS_TT_weeklyconservation_LB.deactivate()
        phase1_inst.TTDS_TT_weeklyconservation_UB.deactivate()
        phase1_inst.TTDS_TT_cumul_weeklyconservation_LB.deactivate()
        phase1_inst.TTDS_TT_cumul_weeklyconservation_UB.deactivate()

    # Post Phase 1 construction taks ------------------------------------------

    # Optionally write out out phase 1 instance
    if bWritePhase1Instance:
        with open(phase1_inst_file, 'w') as f1_inst:
            phase1_inst.pprint(ostream=f1_inst)
            logging.info('Phase 1 instance written')

    # Write out phase 1 problem size info to phase 1 summary file
    with open(phase1_summary_file, 'w') as f1_sum:
        tot_cons = 0
        tot_vars = 0
        f1_sum.write("\n\nConstraint summary \n------------------\n")
        for c in phase1_inst.component_objects(Constraint, active=True):
            f1_sum.write(c.name + " --> " + str(len(c)) + "\n")
            tot_cons += len(c)
        f1_sum.write("\n\nVariable summary \n------------------\n")
        for v in phase1_inst.component_objects(Var):
            f1_sum.write(v.name + " --> " + str(len(v)) + "\n")
            tot_vars += len(v)

        msg = "\ntotal cons = " + str(tot_cons) + "\n"
        f1_sum.write(msg)
        msg = "total vars = " + str(tot_vars) + "\n"
        f1_sum.write(msg)

    # Optionally write out detailed debugging info for start windows (only if width > 0)
    if bWriteStartWinDebug:
        start_win_debug_file = path + scenario + '_debugwin.txt'
        with open(start_win_debug_file,"w") as f_windebug:
            for w in phase1_inst.WEEKS:
                for j in phase1_inst.DAYS:
                    for i in phase1_inst.PERIODS:
                        f_windebug.write(
                            'b_window_epoch[{0},{1},{2}] = {3}, '
                            'e_window_epoch[{0},{1},{2}] = {4}'.format(
                                i, j, w,phase1_inst.b_window_epoch[i, j, w].value,
                                phase1_inst.e_window_epoch[i, j, w].value))

            for (i, j, w) in phase1_inst.PotentialGlobalStartWindow_index:
                if phase1_inst.PotentialGlobalStartWindow[i, j, w]:
                    f_windebug.write(
                        'PotentialGlobalStartWindow[{},{},{}] = {}'.format(
                            i, j, w, phase1_inst.PotentialGlobalStartWindow[i, j, w].value))

            for (i, j, w, k, t) in phase1_inst.PotentialStartWindow_index:
                if phase1_inst.PotentialStartWindow[i, j, w, k, t]:
                    f_windebug.write(
                        'PotentialStartWindow[{},{},{},{},{}] = {}'.format(
                            i, j, w, k, t, phase1_inst.PotentialStartWindow[i, j, w, k, t].value))

            f_windebug.write('okStartWindowRoots_index = ')
            for (t, k) in phase1_inst.okStartWindowRoots_index:
                f_windebug.write(str((t, k)))

            for (t,k) in phase1_inst.okStartWindowRoots_index:
                f_windebug.write('okStartWindowRoots[{},{}] = '.format(t, k))
                for (i, j, w) in phase1_inst.okStartWindowRoots[t, k]:
                    f_windebug.write(str(i, j, w))

            f_windebug.write('okTourType = ')
            for (i,t) in phase1_inst.okTourType:
                f_windebug.write(str(i, t))

            f_windebug.write('okTourTypeDay = ')
            for (i,t,d) in phase1_inst.okTourTypeDay:
                f_windebug.write(str(i, t, d))

            f_windebug.write('bchain echain #links = ')
            for (t,k) in phase1_inst.okStartWindowRoots_index:
                f_windebug.write('(t,k) = [{},{}] = '.format(t, k))
                for (i, j, w) in phase1_inst.bchain[t, k]:
                    out = '({},{},{})'.format(i,j,w)
                    for (x, y, z) in phase1_inst.echain[t, k, i, j, w]:
                        out = out + '({},{},{}) {}\n'.format(x,y,z,phase1_inst.n_links[t, k, i, j, w].value)
                    f_windebug.write(out)

            for (t, k, i, j, w) in phase1_inst.chain_index:
                out = ''
                f_windebug.write('chain[{},{},{},{},{}]='.format(t, k, i, j, w))
                for (x, y, z) in phase1_inst.chain[t, k, i, j, w]:
                    out = out + '({},{},{})*'.format(x, y, z, phase1_inst.chain[t, k, i, j, w].value)
                f_windebug.write(out)

            for (t, k, i, j, w, m) in phase1_inst.link_index:
                f_windebug.write('link[{},{},{},{},{},{}]='.format(t, k, i, j, w, m))

            logging.info('Windows debug info written')

    # solver = pyomo.opt.SolverFactory('cplex')
    # results = solver.solve(self.m, tee=True, keepfiles=False,
    #                        options_string="mip_tolerances_integrality=1e-9 mip_tolerances_mipgap=0")
    #

    # TODO - add status messages during overall model generation and solution process

    # Solve Phase 1 -----------------------------------------------------------

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

    stream_solver = True
    phase1_results = solver.solve(phase1_inst, tee=stream_solver)

    # TODO - check if phase 1 solved
    if phase1_results.solver.status != pyomo.opt.SolverStatus.ok:
        logging.warning('Check solver not ok? Status = %s', phase1_results.solver.status)

    if phase1_results.solver.termination_condition != pyomo.opt.TerminationCondition.optimal:
        logging.warning('Check solver optimality? Term condition = %s',
                        phase1_results.solver.termination_condition)

    logging.info('Phase 1 solution status = %s', phase1_results.solver.status)

    # By default, results are automatically loaded into model instance
    # See https://groups.google.com/forum/#!topic/pyomo-forum/wjjY2XvmG2w

    try:
        phase1_solution_value = phase1_inst.total_cost()
        logging.info('Phase 1 solved successfully')
        logging.info('Phase 1 solution = %s', phase1_solution_value)
    except:
        logging.critical('Phase 1 problem not solved successfully.')
        phase1_solution_status = phase1_results.solver.status
        logging.critical('Status: %s', str(phase1_solution_status))
        sys.exit(1)

    # Write problem size summary file
    with open(phase1_summary_file, 'a') as f1_sum:
        tot_cons = 0
        tot_vars = 0
        f1_sum.write("\n\nConstraint summary \n------------------\n")
        for c in phase1_inst.component_objects(Constraint, active=True):
            # conobj = getattr(phase1_inst, str(c))
            f1_sum.write(c.name + " --> " + str(len(c)) + "\n")
            tot_cons += len(c)
        f1_sum.write("\n\nVariable summary \n------------------\n")
        for v in phase1_inst.component_objects(Var):
            # vobj = getattr(phase1_inst, str(v))
            f1_sum.write(v.name + " --> " + str(len(v)) + "\n")
            tot_vars += len(v)

        msg = "\ntotal cons = " + str(tot_cons) + "\n"
        f1_sum.write(msg)
        msg = "total vars = " + str(tot_vars) + "\n"
        f1_sum.write(msg)

        phase1_solution_status = phase1_results.solver.status
        f1_sum.write(str(phase1_solution_status))

    # Write results file
    with open(phase1_results_file, "w") as f1_res:
        phase1_inst.display(ostream=f1_res)
        logging.info('Phase 1 summary and results written')

    # Write shift summary
    phase1_shiftsummary = write_phase1_shiftsummary(phase1_inst)
    with open(phase1_shiftsum_file,'w') as f1_shiftsum:
        print(phase1_shiftsummary, file=f1_shiftsum)

    # Write tour skeleton
    phase1_tourskeleton = weekenddaysworked_to_tourskeleton(phase1_inst)
    with open(phase1_tourskeleton_file,'w') as f1_tourskeleton:
        print(phase1_tourskeleton, file=f1_tourskeleton)
        phase1_tourskeleton = tourtypeday_to_tourskeleton(phase1_inst)
        print(phase1_tourskeleton, file=f1_tourskeleton)

    # Phase 2 model construction ----------------------------------------------

    # If phase 1 solved, create phase 2 params from phase 1 vars and solve phase 2 problem.
    # TODO - how to check for status other than 'optimal'?

    # Need to convert the phase2 instance to concrete mode before we can add constraints, fix values, etc.
    # phase2_inst.concrete_mode()

    # In general, for concrete models, instance and model object share the same memory space (unlike
    # with an abstract model. Not sure if putting an abstract model into concrete_mode does the same thing.
    # We will reference the instance object to be safe.

    # phase_1_2_integrate(phase1_inst, phase2_inst)

    tot_cap = value(sum(phase1_inst.cov[i,j,w].value for (i,j,w) in phase1_inst.epoch_tuples))
    us1_cost = value(sum(phase1_inst.under1[i,j,w] * phase1_inst.cu1.value for (i,j,w) in phase1_inst.epoch_tuples))
    us2_cost = value(sum(phase1_inst.under2[i,j,w] * phase1_inst.cu2.value for (i,j,w) in phase1_inst.epoch_tuples))

    n_tours = int(round(sum((phase1_inst.TourType[i,t].value for (i,t) in phase1_inst.okTourType))))

    param_n_tours = scalar_to_param('n_tours',n_tours)

    param_Shift = shift_to_param('Shift',phase1_inst)
    param_TourType = tourtype_to_param('TourType',phase1_inst)
    param_TourTypeDay = tourtypeday_to_param('TourTypeDay',phase1_inst)
    param_TourTypeDayShift = tourtypedayshift_to_param('TourTypeDayShift',phase1_inst)
    param_WeekendDaysWorked = weekenddaysworked_to_param('WeekendDaysWorked',phase1_inst)
    param_MultiWeekDaysWorked = multiweekdaysworked_to_param('MultiWeekDaysWorked', phase1_inst)

    param_tour_WIN_TT = tour_WIN_TT_to_param(phase1_inst)

    with open(phase2_dat_file, 'w') as f2_out:
        with open(phase1_dat_file, 'r') as f1_in:
            f1_str = f1_in.read()
            f2_out.write(f1_str)

    # Put the pieces together
    dat = io.StringIO()

    print(param_n_tours, file=dat)
    print(param_Shift, file=dat)
    print(param_TourType, file=dat)
    print(param_TourTypeDay, file=dat)
    print(param_TourTypeDayShift, file=dat)
    print(param_WeekendDaysWorked, file=dat)
    print(param_MultiWeekDaysWorked, file=dat)
    print(param_tour_WIN_TT, file=dat)

    with open(phase2_dat_file, 'a') as f2_dat:
        print(dat.getvalue(), file=f2_dat)
        logging.info('Phase 2 dat file created')

    # Initialize the phase 2 instance

    # phase2_mdl = import_file(phase2_mod_file).model_phase2
    phase2_mdl = mwts_phase2.model
    phase2_inst = phase2_mdl.create_instance(filename = phase2_dat_file)
    phase2_inst.name = 'mwts_phase2_inst'
    logging.info('Phase 2 instance created')

    # Activate/deactivate constraints

    bTour_Weekend_conservation_active = True
    bTour_MWDW_conservation_active = True

    bOneWeekendPatternPerTour_active = True
    bOneMWDWPatternPerTour_active = True

    bTourShift_Weekend_integration1_active = True
    bTourShift_MWDW_integration1_active = True

    bTours_Daily_active = True
    bTours_Daily_conservation_active = True

    # The following four constraints shouldn't be needed now that we have
    # added the Tour_MWDW_con
    bTours_Weekly_LB_active = False
    bTours_Weekly_UB_active = False
    bTours_Total_LB_active = False
    bTours_Total_UB_active = False

    bTours_Shiftlen_Weekly_LB_active = True
    bTours_Shiftlen_Weekly_UB_active = True
    bTours_Shiftlen_Total_LB_active = True
    bTours_Shiftlen_Total_UB_active = True

    bTours_Weekly_Prds_LB_active = True
    bTours_Weekly_Prds_UB_active = True
    bTours_Total_Prds_LB_active = True
    bTours_Total_Prds_UB_active = True

    # Deactivate constraints per the above list of binaries

    if not bTour_Weekend_conservation_active:
        phase2_inst.Tour_Weekend_conservation.deactivate()

    if not bTour_MWDW_conservation_active:
        phase2_inst.Tour_MWDW_conservation.deactivate()

    if not bOneMWDWPatternPerTour_active:
        phase2_inst.OneMWDWPatternPerTour.deactivate()

    if not bOneWeekendPatternPerTour_active:
        phase2_inst.OneWeekendPatternPerTour.deactivate()

    if not bTourShift_Weekend_integration1_active:
        phase2_inst.Tour_MWDW_integration1.deactivate()

    if not bTourShift_MWDW_integration1_active:
        phase2_inst.Tour_MWDW_integration1.deactivate()

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

        with open(phase2_inst_file, 'w') as f2_inst:
            phase2_inst.pprint(ostream=f2_inst)
            logging.info('Phase 2 instance written')

    # Write out phase 2 problem size summary info

    with open(phase2_summary_file,'w') as f2_sum:
        tot_cons = 0
        tot_vars = 0
        f2_sum.write("\n\nConstraint summary \n------------------\n")
        for c in phase2_inst.component_objects(Constraint, active=True):
            # conobj = getattr(phase2_inst, str(c))
            f2_sum.write(c.name + " --> " + str(len(c)) + '\n')
            tot_cons += len(c)
        f2_sum.write("\n\nVariable summary \n------------------\n")
        for v in phase2_inst.component_objects(Var):
            # vobj = getattr(phase2_inst, str(v))
            f2_sum.write(v.name + " --> " + str(len(v)) + '\n')
            tot_vars += len(v)

        msg = "\ntotal cons = " + str(tot_cons)
        f2_sum.write(msg)
        msg = "\ntotal vars = " + str(tot_vars)
        f2_sum.write(msg)

    # Solve phase 2
    stream_solver = True
    phase2_results = solver.solve(phase2_inst, tee=stream_solver)
    phase2_solution_status = str(phase2_results.solver.status)
    logging.info('Phase 2 solution status = %s', str(phase2_solution_status))

    if str(phase2_results.Solution.Status) != 'unknown':

        phase2_solution_status = phase2_results.solver.status
        print(phase2_solution_status)
        print(str(value(phase2_inst.total_num_tours)))
        phase2_solution_value = value(phase2_inst.total_num_tours())

        # Write results file
        with open(phase2_results_file, "w") as f2_res:
            phase2_inst.display(ostream=f2_res)
            logging.info('Phase 2 summary and results written')

    # Create the multi-week tour file
    idx_copy = []
    for idx in phase2_inst.TourShift:
        idx_copy.append(idx)
    idx_sorted = sorted(idx_copy)

    with open(tour_file, "w") as f2_tour:
        f2_tour.write('n_tours {}\n'.format(phase2_inst.n_tours.value))
        f2_tour.write('n_weeks {}\n'.format(phase2_inst.n_weeks.value))
        f2_tour.write('PP4\n')

        for idx in idx_sorted:
            if phase2_inst.TourShift[idx].value > 0:
                f2_tour.write(
                    '{} {} {} {} {} {} {}\n'.format(idx[0], idx[5], idx[3], idx[1], idx[2],
                                                    phase2_inst.lengths[idx[4]],
                                                    phase2_inst.WIN_x[idx[0]]))

        f2_tour.write('TTS\n')
        for tour in phase2_inst.TOURS:
            f2_tour.write('{}\n'.format(phase2_inst.TT_x[tour]))

    create_mwt(tour_file, scenario, path)
    logging.info('Tour related output files created')
    ts_now = datetime.datetime.now()
    var_tuple = (scenario, phase1_solution_value, tot_cap,
                 str(phase1_solution_status), str(phase2_solution_status),
                 us1_cost, us2_cost, phase2_solution_value,
                 str(ts_now))
    logging.info('Solution log record %s', str(var_tuple))

    # Connect to the problem solution log database.
    if results_db is not None:
        conn = sqlite3.connect(results_db)
        cur = conn.cursor()
        now = datetime.datetime.now()
        field_name_tuple = '(problem,MIP_obj,tot_cap,phase1_sol_status,phase2_sol_status,' \
                           'us1_cost,us2_cost,phase2_obj,timestamp)'
        var_tuple = (scenario, phase1_solution_value, tot_cap,
                     str(phase1_solution_status), str(phase2_solution_status), us1_cost, us2_cost,
                     phase2_solution_value, str(now))
        sql_insert = 'insert into solution_log ' + field_name_tuple + ' values ' + str(var_tuple)
        print(sql_insert)
        cur.execute(sql_insert)
        conn.commit()

        # Update the problem list table
        cur = conn.cursor()
        conn.row_factory = sqlite3.Row
        sql_update = 'update problem_list set sol_status="' + str(phase1_solution_status) + ':' + str(
            phase2_solution_status) + '" where problem="' + scenario + '"'
        cur.execute(sql_update)
        print(sql_update)
        conn.commit()
        conn.close()


def probe_phase2(scenario, phase2_dat_file, path,
                 which_solver, timelimit, mipgap, wintt_filter=None):
    """
    Created this to make it easy to try out limited combinations of windows and tour
    types to try to debug model. Needed to manually add lines such as

    set activeTT := 8;
    set activeWIN := 44;

    to Phase 2 dat file and then resolve Phase 2. Couldn't find easy way to integrate with
    main model above, so created this probe_phase2() function.

    :param scenario:
    :param phase2_dat_file:
    :param path:
    :param which_solver:
    :param timelimit:
    :param mipgap:
    :param phase2_results_file:
    :return:
    """

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

    # Initialize the phase 2 instance

    # phase2_mdl = import_file(phase2_mod_file).model_phase2
    phase2_mdl = mwts_phase2.model_phase2
    phase2_inst = phase2_mdl.create_instance(filename=phase2_dat_file)
    phase2_inst.name = 'mwts_phase2_inst'
    logging.info('Phase 2 instance created')

    # Optionally limit which windows and tour types


    # Activate/deactivate constraints

    bTour_Weekend_conservation_active = True
    bTour_MWDW_conservation_active = True

    bOneWeekendPatternPerTour_active = True
    bOneMWDWPatternPerTour_active = True

    bTourShift_Weekend_integration1_active = True
    bTourShift_MWDW_integration1_active = True

    bTours_Daily_active = True
    bTours_Daily_conservation_active = True

    bTour_WkendDof_conservation_active = True


    bTours_Weekly_LB_active = False
    bTours_Weekly_UB_active = False
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

    if bTour_Weekend_conservation_active:
        phase2_inst.Tour_Weekend_conservation.deactivate()

    if bTour_MWDW_conservation_active:
        phase2_inst.Tour_MWDW_conservation.deactivate()

    if not bOneWeekendPatternPerTour_active:
        phase2_inst.OneWeekendPatternPerTour.deactivate()

    if not bOneMWDWPatternPerTour_active:
        phase2_inst.OneMWDWPatternPerTour.deactivate()

    if not bTourShift_Weekend_integration1_active:
        phase2_inst.TourShift_Weekend_integration1.deactivate()

    if not bTourShift_MWDW_integration1_active:
        phase2_inst.TourShift_MWDW_integration1.deactivate()

    if not bTours_Daily_active:
        phase2_inst.Tours_Daily.deactivate()



    if not bTour_WkendDof_conservation_active:
        phase2_inst.Tour_WkendDof_conservation.deactivate()





    if not bTours_Daily_conservation_active:
        phase2_inst.Tours_Daily_conservation.deactivate()

    if not bTours_Weekly_active:
        phase2_inst.Tours_Weekly.deactivate()

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

    # Solve phase 2
    stream_solver = True
    phase2_results = solver.solve(phase2_inst, tee=stream_solver)

    phase2_solution_status = phase2_results.solver.status
    print(phase2_solution_status)
    print(str(value(phase2_inst.total_num_tours)))
    phase2_solution_value = value(phase2_inst.total_num_tours())

    # Write results file
    phase2_results_file = path + scenario + '.yml'
    with open(phase2_results_file, "w") as f2_res:
        phase2_inst.display(ostream=f2_res)
        logging.info('Phase 2 summary and results written')

    # Create the multi-week tour file
    tour_file = path + scenario + '.tur'
    idxcopy = []
    for idx in phase2_inst.TourShift:
        idxcopy.append(idx)
    idxsorted = sorted(idxcopy)

    with open(tour_file, "w") as f2_tour:
        f2_tour.write('n_tours {}\n'.format(phase2_inst.n_tours.value))
        f2_tour.write('n_weeks {}\n'.format(phase2_inst.n_weeks.value))
        f2_tour.write('PP4\n')

        for idx in idxsorted:
            if phase2_inst.TourShift[idx].value > 0:
                f2_tour.write(
                    '{} {} {} {} {} {} {}\n'.format(idx[0], idx[5], idx[3], idx[1], idx[2],
                                                    phase2_inst.lengths[idx[4]],
                                                    phase2_inst.WIN_x[idx[0]]))

        f2_tour.write('TTS\n')
        for tour in phase2_inst.TOURS:
            f2_tour.write('{}\n'.format(phase2_inst.TT_x[tour]))

    create_mwt(tour_file, scenario, path)
