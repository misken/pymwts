"""
Solve both phases of the mwts model. 

The solution from phase 1 is fed
to the phase 2 model to create tours. Phase 1 solution values for variables are converted
to GMPL DAT parameter statements and a phase 2 DAT file is created from a copy of the
phase 1 DAT file and these additional parameters.
"""

import sys
import io
import sqlite3
import datetime
import logging

import pyomo.opt
import pyomo.environ as pyo

import pymwts.mwts_phase1 as phase1
import pymwts.mwts_phase2 as phase2
import pymwts.mwts_utils as mwts_utils
import pymwts.pymwtsio.mwts_output as mwts_output
from pymwts.pymwtsio.mwts_makedat import scalar_to_param

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
              prob_num=-1,
              results_db=None,
              debug_start_windows=False,
              write_phase1_instance=False,
              write_phase2_instance=False,
              force_solve=False):
    """

    :param scenario:
    :param phase1_dat_file:
    :param path:
    :param which_solver:
    :param timelimit:
    :param mipgap:
    :param results_db:
    :param debug_start_windows:
    :param write_phase1_instance:
    :param write_phase2_instance:
    :param force_solve:
    :return:
    """

    # Check to see if this problem has already been run
    # Connect to the problem solution log database.
    if not force_solve and (results_db is not None):
        conn = sqlite3.connect(results_db)
        cur = conn.cursor()
        row = cur.execute('SELECT Problem, sol_status FROM problem_list WHERE Problem=?',
                          (scenario,))

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
    phase1_capsum_file = path + scenario + '_phase1_capsum.csv'
    phase1_tourskeleton_file = path + scenario + '_phase1_tourskeleton.csv'

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

    print('\n*** Scenario {}\n'.format(scenario))
    logging.info('Scenario %s started', scenario)
    logging.info('DAT %s', phase1_dat_file)

    # Create Phase 1 model instance
    phase1_mdl = phase1.model
    phase1_inst = phase1_mdl.create_instance(filename=phase1_dat_file)
    phase1_inst.name = 'mwts_phase1_inst'
    logging.info('Phase 1 instance created')
    print('\n*** Phase 1 model instance created.\n')

    # Activate/deactivate constraints -----------------------------------------

    # Chain debugging
    # is_chains_sweep_l_con_active = False
    # is_chains_sweep_u_con_active = False
    # is_chains_tot_con_active = False
    #
    # if not is_chains_sweep_l_con_active:
    #     phase1_inst.chains_sweep_l_con.deactivate()
    #
    # if not is_chains_sweep_u_con_active:
    #     phase1_inst.chains_sweep_u_con.deactivate()
    #
    # if not is_chains_tot_con_active:
    #     phase1_inst.chains_tot_con.deactivate()

    # Deactivate part-time fraction upper bound if all tour types are part-time or all full-time
    tot_parttime_ttypes = sum(phase1_inst.tt_parttime[t] for t in phase1_inst.activeTT)
    if tot_parttime_ttypes == len(phase1_inst.activeTT) or \
            phase1_inst.max_parttime_frac.value == 1.0 or \
            tot_parttime_ttypes == 0:
        phase1_inst.max_ptfrac_con.deactivate()

    # Boolean indicators for debugging weekends related constraints
    b_weekend_subsets_5_4_con2_active = True
    b_weekend_subsets_4_3_con2_active = True
    b_weekend_subsets_3_2_con2_active = True
    b_weekend_subsets_2_1_con2_active = True

    # Boolean indicators for possible redundant constraints
    b_TTD_TT_weeklyconservation_active = False
    b_TTDS_TT_weeklyconservation_active = False
    b_prds_worked_shiflen_weekly_active = False

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

    if not b_prds_worked_shiflen_weekly_active:
        phase1_inst.prds_worked_shiflen_weekly_LB.deactivate()
        phase1_inst.prds_worked_shiflen_weekly_UB.deactivate()
        phase1_inst.prds_worked_cumul_shiflen_weekly_LB.deactivate()
        phase1_inst.prds_worked_cumul_shiflen_weekly_UB.deactivate()

    # Post Phase 1 construction tasks ------------------------------------------

    # Optionally write out out phase 1 instance
    if write_phase1_instance:
        with open(phase1_inst_file, 'w') as f1_inst:
            phase1_inst.pprint(ostream=f1_inst)
            logging.info('Phase 1 instance written')

    # Write out phase 1 problem size info to phase 1 summary file
    with open(phase1_summary_file, 'w') as f1_sum:
        tot_cons = 0
        tot_vars = 0
        f1_sum.write("\n\nConstraint summary \n------------------\n")
        for c in phase1_inst.component_objects(pyo.Constraint, active=True):
            f1_sum.write(c.name + " --> " + str(len(c)) + "\n")
            tot_cons += len(c)
        f1_sum.write("\n\nVariable summary \n------------------\n")
        for v in phase1_inst.component_objects(pyo.Var):
            f1_sum.write(v.name + " --> " + str(len(v)) + "\n")
            tot_vars += len(v)

        msg = "\ntotal cons = " + str(tot_cons) + "\n"
        f1_sum.write(msg)
        msg = "total vars = " + str(tot_vars) + "\n"
        f1_sum.write(msg)

    # Optionally write out detailed debugging info for start 
    # windows (only if width > 0)
    if debug_start_windows:
        start_win_debug_file = path + scenario + '_debugwin.txt'
        with open(start_win_debug_file, "w") as f_debug_win:
            f_debug_win.write('\nb_window_epochs and e_window_epochs\n')
            for w in phase1_inst.WEEKS:
                for j in phase1_inst.DAYS:
                    for i in phase1_inst.PERIODS:
                        f_debug_win.write(
                            'b_window_epoch[{0},{1},{2}] = {3}, '
                            'e_window_epoch[{0},{1},{2}] = {4}\n'.format(
                                i, j, w, phase1_inst.b_window_epoch[i, j, w],
                                phase1_inst.e_window_epoch[i, j, w]))
            f_debug_win.write('\nEND b_window_epochs and e_window_epochs\n')

            f_debug_win.write('\nPotentialGlobalStartWindow[i, j, w]\n')
            for (i, j, w) in phase1_inst.epoch_tuples:
                if phase1_inst.PotentialGlobalStartWindow[i, j, w]:
                    f_debug_win.write(
                        'PotentialGlobalStartWindow[{},{},{}] = {}\n'.format(
                            i, j, w, list(phase1_inst.PotentialGlobalStartWindow[i, j, w])))
            f_debug_win.write('\nEND PotentialGlobalStartWindow[i, j, w]\n')

            f_debug_win.write('\nPotentialStartWindow[i, j, w, k, t]\n')
            for (i, j, w, k, t) in phase1_inst.PotentialStartWindow_idx:
                if phase1_inst.PotentialStartWindow[i, j, w, k, t]:
                    f_debug_win.write(
                        '\nPotentialStartWindow[{},{},{},{},{}] = {}\n'.format(
                            i, j, w, k, t,
                            list(phase1_inst.PotentialStartWindow[i, j, w, k, t])))
            f_debug_win.write('\nEND PotentialStartWindow[i, j, w, k, t]\n')

            f_debug_win.write('okStartWindowRoots_idx = ')
            for (t, k) in phase1_inst.okStartWindowRoots_idx:
                f_debug_win.write('({},{})\n'.format(t, k))
            f_debug_win.write('\nEND okStartWindowRoots_idx = ')

            f_debug_win.write('\nokStartWindowRoots\n')
            for (t, k) in phase1_inst.okStartWindowRoots_idx:
                f_debug_win.write('okStartWindowRoots[{},{}] = \n'.format(t, k))
                for (i, j, w) in phase1_inst.okStartWindowRoots[t, k]:
                    f_debug_win.write('({},{},{})\n'.format(i, j, w))
            f_debug_win.write('\nEND okStartWindowRoots\n')

            f_debug_win.write('\nokTourType = ')
            for (i, t) in phase1_inst.okTourType:
                f_debug_win.write('({},{})\n'.format(i, t))
            f_debug_win.write('\nEND okTourType = ')

            f_debug_win.write('\nokTourTypeDay = ')
            for (i, t, d) in phase1_inst.okTourTypeDay:
                f_debug_win.write('({},{},{})\n'.format(i, t, d))
            f_debug_win.write('\nENDokTourTypeDay = ')

            f_debug_win.write('\nbchain echain links\n')
            for (t, k) in phase1_inst.okStartWindowRoots_idx:
                f_debug_win.write('\n(t,k) = [{},{}] = \n'.format(t, k))
                for (i, j, w) in phase1_inst.bchain[t, k]:
                    out = '({},{},{})'.format(i, j, w)
                    for (x, y, z) in phase1_inst.echain[t, k, i, j, w]:
                        out = out + '({},{},{}) {}\n'.format(
                            x, y, z, phase1_inst.n_links[t, k, i, j, w])
                    f_debug_win.write(out)
            f_debug_win.write('\nEND bchain echain links\n')

            f_debug_win.write('\nchains\n')
            for (t, k, i, j, w) in phase1_inst.chain_idx:
                out = ''
                f_debug_win.write('\nchain[{},{},{},{},{}]=\n'.format(t, k, i, j, w))
                for (x, y, z) in phase1_inst.chain[t, k, i, j, w]:
                    out = out + '({},{},{})*'.format(x, y, z, phase1_inst.chain[t, k, i, j, w])
                f_debug_win.write(out)
            f_debug_win.write('\nEND chains\n')

            f_debug_win.write('\nlinkspan\n')
            for (t, k, i, j, w, m) in phase1_inst.link_idx:
                f_debug_win.write('\nlinkspan[{},{},{},{},{},{}]={}\n'.format(t, k, i, j, w, m,
                                                                              list(
                                                                                  phase1_inst.linkspan[
                                                                                      t, k, i, j, w, m])))
            f_debug_win.write('\nEND linkspan\n')

            f_debug_win.write('\nchains_sweep_l_idx\n')
            for (t, k, b, j, w, p, v) in phase1_inst.chains_sweep_l_idx:
                f_debug_win.write('\n[{},{},{},{},{},{},{}]'.format(t, k, b, j, w, p, v))
            f_debug_win.write('\nEND chains_sweep_l_idx\n')
            # chains_sweep_l_idx

            f_debug_win.write('\nchains_sweep_u_idx\n')
            for (t, k, b, j, w, p, v) in phase1_inst.chains_sweep_l_idx:
                f_debug_win.write('\n[{},{},{},{},{},{},{}]'.format(t, k, b, j, w, p, v))
            f_debug_win.write('\nEND chains_sweep_u_idx\n')

            logging.info('Windows debug info written')

    # TODO - add status messages during overall model generation and solution process

    # Solve Phase 1 -----------------------------------------------------------

    # Setup the solver
    print('\n*** Setting up the solver.\n')
    solver = None
    if which_solver == 'cbc':
        solver = pyomo.opt.SolverFactory('cbc')
        solver.options.seconds = timelimit
        solver.options.ratioGap = mipgap
    # solver.options.solution = '../tests/solution.sol'
    # This is NOT the correct way to specify this option

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

    # By default, results are automatically loaded into model instance
    # See https://groups.google.com/forum/#!topic/pyomo-forum/wjjY2XvmG2w
    # In order to check if we found a solution before time limit reached,
    # need to override this default behavior with load_solutions=False.

    print('\n*** Starting to solve Phase 1 model\n')
    stream_solver = True
    phase1_results = solver.solve(phase1_inst, tee=stream_solver,
                                  load_solutions=False)

    if len(phase1_results.solution) == 0:
        print('\n*** No Phase 1 solution found. See log file for details.\n')
        logging.warning('No Phase 1 solution found. Status = %s Termination condition = %s',
                        phase1_results.solver.status, phase1_results.solver.termination_condition)

        ts_now = datetime.datetime.now()
        phase1_solution_status = 'no_int_soln_found'
        phase2_solution_status = 'na'
        var_tuple = (scenario, 0.0, 0.0,
                     phase1_solution_status, phase2_solution_status,
                     0.0, 0.0, 0.0,
                     str(ts_now))
        logging.info('Solution log record %s', str(var_tuple))

        # Connect to the problem solution log database.
        if results_db is not None:
            conn = sqlite3.connect(results_db)
            cur = conn.cursor()
            field_name_tuple = '(problem,MIP_obj,tot_cap,phase1_sol_status,phase2_sol_status,' \
                               'us1_cost,us2_cost,phase2_obj,timestamp)'
            sql_insert = 'insert into solution_log ' + field_name_tuple + ' values ' + str(var_tuple)
            print(sql_insert)
            cur.execute(sql_insert)
            conn.commit()

            # Update the problem list table
            cur = conn.cursor()
            conn.row_factory = sqlite3.Row
            sql_update = 'update problem_list set sol_status="' + phase1_solution_status + ':' + \
                phase2_solution_status + '" where problem="' + scenario + '"'
            cur.execute(sql_update)
            print(sql_update)
            conn.commit()
            conn.close()

        sys.exit(1)
    else:
        print('\n*** Phase 1 solution found\n')
        logging.info('Phase 1 solution found. Status = %s Termination condition = %s',
                     phase1_results.solver.status, phase1_results.solver.termination_condition)

    # Some sort of Phase 1 solution was found, let's load the results into model
    phase1_inst.solutions.load_from(phase1_results)

    if phase1_results.solver.status != pyomo.opt.SolverStatus.ok:
        logging.warning('Check solver not ok? Status = %s', phase1_results.solver.status)

    if phase1_results.solver.termination_condition != pyomo.opt.TerminationCondition.optimal:
        logging.warning('Check solver optimality? Term condition = %s',
                        phase1_results.solver.termination_condition)

    if pyo.value(phase1_inst.total_cost, exception=False) is not None:
        phase1_solution_value = pyo.value(phase1_inst.total_cost)
        logging.info('Phase 1 solved successfully')
        logging.info('Phase 1 solution = %s', phase1_solution_value)
    else:
        logging.critical('Phase 1 problem not solved successfully.')
        phase1_solution_status = phase1_results.solver.status
        logging.critical('Status: %s', str(phase1_solution_status))
        sys.exit(1)

    # Write problem size summary file
    with open(phase1_summary_file, 'a') as f1_sum:
        tot_cons = 0
        tot_vars = 0
        f1_sum.write("\n\nConstraint summary \n------------------\n")
        for c in phase1_inst.component_objects(pyo.Constraint, active=True):
            f1_sum.write(c.name + " --> " + str(len(c)) + "\n")
            tot_cons += len(c)
        f1_sum.write("\n\nVariable summary \n------------------\n")
        for v in phase1_inst.component_objects(pyo.Var):
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
    phase1_shiftsummary = mwts_output.write_phase1_shiftsummary(phase1_inst)
    with open(phase1_shiftsum_file, 'w') as f1_shiftsum:
        print(phase1_shiftsummary, file=f1_shiftsum)

    # Write capacity summary
    phase1_capsummary = mwts_output.write_phase1_capsummary(phase1_inst)
    with open(phase1_capsum_file, 'w') as f1_capsum:
        print(phase1_capsummary, file=f1_capsum)

    # Write tour skeleton
    phase1_tourskeleton = mwts_output.weekenddaysworked_to_tourskeleton(phase1_inst)
    with open(phase1_tourskeleton_file, 'w') as f1_tourskeleton:
        print(phase1_tourskeleton, file=f1_tourskeleton)
        phase1_tourskeleton = mwts_output.tourtypeday_to_tourskeleton(phase1_inst)
        print(phase1_tourskeleton, file=f1_tourskeleton)

    # Phase 2 model construction ----------------------------------------------

    # If phase 1 solved, create phase 2 params from phase 1 vars and solve
    # phase 2 problem.
    # TODO - how to check for status other than 'optimal'?

    # Need to convert the phase2 instance to concrete mode before we can add constraints, fix values, etc.
    # phase2_inst.concrete_mode()

    # In general, for concrete models, instance and model object share the same memory space (unlike
    # with an abstract model. Not sure if putting an abstract model into concrete_mode does the same thing.
    # We will reference the instance object to be safe.

    # phase_1_2_integrate(phase1_inst, phase2_inst)

    tot_cap = pyo.value(sum(phase1_inst.cov[i,j,w].value
                            for (i,j,w) in phase1_inst.epoch_tuples))

    tot_dmd = pyo.value(sum(phase1_inst.dmd_staff[i,j,w]
                            for (i,j,w) in phase1_inst.epoch_tuples))

    tot_us = pyo.value(sum(phase1_inst.under1[i,j,w] + phase1_inst.under2[i,j,w]
                            for (i,j,w) in phase1_inst.epoch_tuples))

    us1_cost = pyo.value(sum(phase1_inst.under1[i,j,w] * phase1_inst.cu1.value
                             for (i,j,w) in phase1_inst.epoch_tuples))
    us2_cost = pyo.value(sum(phase1_inst.under2[i,j,w] * phase1_inst.cu2.value
                             for (i,j,w) in phase1_inst.epoch_tuples))

    n_tours = int(round(sum((phase1_inst.TourType[i, t].value
                             for (i,t) in phase1_inst.okTourType))))

    param_n_tours = scalar_to_param('n_tours', n_tours)

    param_Shift = mwts_utils.shift_to_param('Shift', phase1_inst)
    param_TourType = mwts_utils.tourtype_to_param('TourType', phase1_inst)
    param_TourTypeDay = mwts_utils.tourtypeday_to_param(
        'TourTypeDay', phase1_inst)
    param_TourTypeDayShift = mwts_utils.tourtypedayshift_to_param(
        'TourTypeDayShift', phase1_inst)
    param_WeekendDaysWorked = mwts_utils.weekenddaysworked_to_param(
        'WeekendDaysWorked', phase1_inst)
    param_MultiWeekDaysWorked = mwts_utils.multiweekdaysworked_to_param(
        'MultiWeekDaysWorked', phase1_inst)

    param_tour_WIN_TT = mwts_utils.tour_WIN_TT_to_param(phase1_inst)

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
        logging.info('Phase 2 dat file created.')

    # Initialize the phase 2 instance
    phase2_mdl = phase2.model
    phase2_inst = phase2_mdl.create_instance(filename=phase2_dat_file)
    phase2_inst.name = 'mwts_phase2_inst'
    logging.info('Phase 2 instance created')
    print('\n*** Phase 2 model instance created.\n')

    # Activate/deactivate constraints
    bTour_Weekend_conservation_active = True
    bTour_MWDW_conservation_active = True

    bOneWeekendPatternPerTour_active = True
    bOneMWDWPatternPerTour_active = True

    bTourShift_Weekend_integration1_active = True
    bTourShift_MWDW_integration1_active = True

    bTours_Daily_active = True
    bTours_Daily_conservation_active = True

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
    if write_phase2_instance:

        with open(phase2_inst_file, 'w') as f2_inst:
            phase2_inst.pprint(ostream=f2_inst)
            logging.info('Phase 2 instance written')

    # Write out phase 2 problem size summary info

    with open(phase2_summary_file,'w') as f2_sum:
        tot_cons = 0
        tot_vars = 0
        f2_sum.write("\n\nConstraint summary \n------------------\n")
        for c in phase2_inst.component_objects(pyo.Constraint, active=True):
            f2_sum.write(c.name + " --> " + str(len(c)) + '\n')
            tot_cons += len(c)
        f2_sum.write("\n\nVariable summary \n------------------\n")
        for v in phase2_inst.component_objects(pyo.Var):
            f2_sum.write(v.name + " --> " + str(len(v)) + '\n')
            tot_vars += len(v)

        msg = "\ntotal cons = " + str(tot_cons)
        f2_sum.write(msg)
        msg = "\ntotal vars = " + str(tot_vars)
        f2_sum.write(msg)

    # Solve phase 2
    print('\n*** Starting to solve Phase 2 model.\n')
    stream_solver = True
    phase2_results = solver.solve(phase2_inst, tee=stream_solver,
                                  load_solutions=False)

    if len(phase2_results.solution) == 0:
        logging.warning('No Phase 2 solution found. Status = %s Termination condition = %s',
                        phase2_results.solver.status, phase2_results.solver.termination_condition)
    else:
        logging.info('Phase 2 solution found. Status = %s Termination condition = %s',
                     phase2_results.solver.status, phase2_results.solver.termination_condition)

    # Some sort of Phase 2 solution was found, let's load the results into model
    phase2_inst.solutions.load_from(phase2_results)

    if str(phase2_results.Solution.Status) != 'unknown':
        print('\n*** Phase 2 model solved.\n')
        phase2_solution_status = phase2_results.solver.status
        # print(phase2_solution_status)
        # print(str(pyo.value(phase2_inst.total_shifts)))
        phase2_solution_value = pyo.value(phase2_inst.total_shifts())
        # Write results file
        with open(phase2_results_file, "w") as f2_res:
            phase2_inst.display(ostream=f2_res)
            logging.info('Phase 2 summary and results written')

    # Create tour related output files
    prds_per_fte = 40.0 * phase2_inst.n_prds_per_day.value / 24.0
    mwts_output.write_phase2_tours(phase2_inst, prds_per_fte, tot_dmd, tot_us, scenario, path)

    logging.info('Tour related output files created')
    print('\n*** Output files created.\n')
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
        field_name_tuple = '(problem,MIP_obj,tot_cap,phase1_sol_status,phase2_sol_status,' \
                           'us1_cost,us2_cost,phase2_obj,timestamp)'
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


