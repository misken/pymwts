"""
Phase 1 for implicit multi-week tour scheduling model
"""

# Author: misken
# License: MIT License

import pyomo.environ as pyo

from pymwts.mwts_shared import epoch_to_tuple, epoch_increment
import pymwts.mwts_shared as mwts_shared


# TODO Would be nice if Phase 1 and Phase 2 could share base pyo.Parameters
# model = import_file('mwts_basepyo.Params.py').model

# Create Phase 1 abstract model
model = pyo.AbstractModel()
model.name = "mwts_phase1"

# Constants
infinity = float('inf')

# General temporal parameters -------------------------------------------------

model.n_prds_per_day = \
    pyo.Param(within=pyo.PositiveIntegers)  # n_P

model.n_days_per_week = \
    pyo.Param(within=pyo.PositiveIntegers)  # 7

model.n_weeks = \
    pyo.Param(within=pyo.PositiveIntegers)  # n_W

model.n_prds_per_week = pyo.Param(within=pyo.PositiveIntegers,
                                  initialize=mwts_shared.n_prds_per_week_init)

model.n_prds_per_cycle = \
    pyo.Param(within=pyo.PositiveIntegers,
              initialize=mwts_shared.n_prds_per_cycle_init)

# Temporal ranges used for various index sets
model.PERIODS = pyo.RangeSet(1, model.n_prds_per_day)
model.EPOCHS = pyo.RangeSet(1, model.n_prds_per_cycle)
model.WINDOWS = pyo.RangeSet(1, model.n_prds_per_day)
model.DAYS = pyo.RangeSet(1, model.n_days_per_week)
model.WEEKS = pyo.RangeSet(1, model.n_weeks)
model.WEEKENDS = pyo.RangeSet(1, 3) # For now, just implementing type 1 weekends (Sat and Sun)

model.epoch = pyo.Param(model.PERIODS, model.DAYS,
                        model.WEEKS, initialize=mwts_shared.epoch_init)

model.epoch_tuples = model.PERIODS * model.DAYS * model.WEEKS

# Coverage related parameters -------------------------------------------------

# Target and minimum staffing levels - this is week specific.
model.dmd_staff = pyo.Param(model.PERIODS, model.DAYS,
                            model.WEEKS)

model.min_staff = pyo.Param(model.PERIODS, model.DAYS,
                            model.WEEKS)

# ### Tour type related parameters

# Shift Lengths ---------------------------------------------------------------
# Number of different shift lengths
model.n_lengths = pyo.Param(within=pyo.PositiveIntegers)
# Range of length indexes
model.LENGTHS = pyo.RangeSet(1, model.n_lengths)
# Vector of shift lengths
model.lengths = pyo.Param(model.LENGTHS)

# Tour Types ------------------------------------------------------------------
# Number of different tour types
model.n_tts = pyo.Param(within=pyo.PositiveIntegers)
# Range of tour type indexes
model.TTYPES = pyo.RangeSet(1, model.n_tts)
# Set of allowable length indices by tour type
model.tt_length_x = pyo.Set(model.TTYPES, ordered=True, )
# Bounds on tour type pyo.Variables
model.tt_lb = pyo.Param(model.TTYPES)
model.tt_ub = pyo.Param(model.TTYPES, default=infinity)

model.activeTT = pyo.Set(dimen=1,
                         ordered=True, initialize=mwts_shared.activeTT_init)

# Bounds on days and shifts worked over the week ------------------------------

# 2a. Min and max number of days worked by week by shiftlen by tour type
model.tt_shiftlen_min_dys_weeks = pyo.Param(model.TTYPES,
                                            model.LENGTHS,
                                            model.WEEKS, default=0)

model.tt_shiftlen_max_dys_weeks = pyo.Param(model.TTYPES,
                                            model.LENGTHS,
                                            model.WEEKS, default=1e+6)

# 2b. Min and max number of days worked by cumulative weeks by shiftlen by tour type
model.tt_shiftlen_min_cumul_dys_weeks = pyo.Param(model.TTYPES,
                                                  model.LENGTHS,
                                                  model.WEEKS, default=0)

model.tt_shiftlen_max_cumul_dys_weeks = pyo.Param(model.TTYPES,
                                                  model.LENGTHS,
                                                  model.WEEKS, default=1e+6)

# 3a. Min and max number of periods worked by week by tour type
model.tt_min_prds_weeks = pyo.Param(model.TTYPES,
                                    model.WEEKS, default=0)

model.tt_max_prds_weeks = pyo.Param(model.TTYPES,
                                    model.WEEKS, default=1e+6)

# 3b. Min and max number of periods worked by cumulative weeks by tour type
model.tt_min_cumul_prds_weeks = pyo.Param(model.TTYPES,
                                          model.WEEKS, default=0)

model.tt_max_cumul_prds_weeks = pyo.Param(model.TTYPES,
                                          model.WEEKS, default=1e+6)


# 1a. Min and max number of days worked by week by tour type. These are computed
# convenience parameters and are used in the weekend subset constraints.

def tt_min_dys_weeks_init(M, t, w):
    return sum(M.tt_shiftlen_min_dys_weeks[t, k, w] for k in M.tt_length_x[t])

def tt_max_dys_weeks_init(M, t, w):
    return sum(M.tt_shiftlen_max_dys_weeks[t, k, w] for k in M.tt_length_x[t])

model.tt_min_dys_weeks = pyo.Param(model.TTYPES,
                                   model.WEEKS, initialize=tt_min_dys_weeks_init)

model.tt_max_dys_weeks = pyo.Param(model.TTYPES,
                                   model.WEEKS, initialize=tt_max_dys_weeks_init)

# Legacy parameters - the following are holdovers from previous model versions.
#    Leaving these in so that model still works with older data files
#    containing these parameters.



# 1b. Min and max number of days worked by cumulative weeks by tour type
model.tt_min_cumul_dys_weeks = pyo.Param(model.TTYPES,
                                         model.WEEKS, default=0.0)

model.tt_max_cumul_dys_weeks = pyo.Param(model.TTYPES,
                                         model.WEEKS, default=1e+6)

# 4a. Min and max number of periods worked by week by tour type by shift length
model.tt_shiftlen_min_prds_weeks = pyo.Param(model.TTYPES,
                                             model.LENGTHS,
                                             model.WEEKS, default=0)

model.tt_shiftlen_max_prds_weeks = pyo.Param(model.TTYPES,
                                             model.LENGTHS,
                                             model.WEEKS, default=1e+6)

# 4b. Min and max number of periods worked by cumulative weeks by tour type by shift length
model.tt_shiftlen_min_cumul_prds_weeks = pyo.Param(model.TTYPES,
                                                   model.LENGTHS,
                                                   model.WEEKS, default=0)

model.tt_shiftlen_max_cumul_prds_weeks = pyo.Param(model.TTYPES,
                                                   model.LENGTHS,
                                                   model.WEEKS, default=1e+6)

# ------ End of legacy parameters

# Indicator for part-time tour types: 1 for part-time, 0 for full-time
model.tt_parttime = pyo.Param(model.TTYPES)

# Allowable shift start times -------------------------------------------------

# Note that these are tour type and shift length specific
model.allow_start = pyo.Param(model.PERIODS,
                              model.DAYS,
                              model.LENGTHS,
                              model.TTYPES,
                              default=0.0)

model.okShifts = pyo.Set(dimen=5, ordered=True,
                         initialize=mwts_shared.okShifts_rule)

model.okShiftTypes = pyo.Set(
    dimen=4, ordered=True, initialize=mwts_shared.okShiftTypes_rule)

# Limits on part time labor and limits on total labor -------------------------

# Maximum fraction of labor hours covered by part-time employees
model.max_parttime_frac = pyo.Param()
# Maximum labor expenditure
model.labor_budget = pyo.Param()

# Cost related parameters -----------------------------------------------------

# Tour type differential
model.tt_cost_multiplier = pyo.Param(model.TTYPES)

# Understaffing cost
model.cu1 = pyo.Param()  # Tier 1 understaffing cost
model.cu2 = pyo.Param()  # Tier 2 understaffing cost
model.usb = pyo.Param()  # Tier 1 understaffing upper bound

# Weekend related parameters --------------------------------------------------

# Need to generalize for any number of periods per day
model.midnight_thresh = pyo.Param(model.TTYPES, default=1e+6)

model.weekend = pyo.Set(model.WINDOWS, model.TTYPES,
                        ordered=True, initialize=mwts_shared.weekend_init)

model.weekend_type = pyo.Param(model.WINDOWS, model.TTYPES,
                               initialize=mwts_shared.weekend_type_init)

model.max_weekend_patterns = pyo.Param(
    initialize=mwts_shared.max_wkend_patterns_init)

model.num_weekend_patterns = pyo.Param(model.WEEKENDS,
                                       model.TTYPES, default=0)

model.A_wkend_days_idx = pyo.Set(
    dimen=5, ordered=True, initialize=mwts_shared.A_wkend_days_idx_rule)

model.A_wkend_days = pyo.Param(model.A_wkend_days_idx,
                               within=pyo.Boolean, default=0)

# Multiweek days worked patterns ----------------------------------------------

model.max_mwdw_patterns = pyo.Param(initialize=mwts_shared.max_mwdw_init)
model.num_mwdw_patterns = pyo.Param(model.TTYPES)

model.A_mwdw_idx = pyo.Set(dimen=3, ordered=True,
                           initialize=mwts_shared.A_mwdw_idx_rule)

model.A_mwdw = pyo.Param(model.A_mwdw_idx,
                         within=pyo.NonNegativeIntegers, default=0)


# Phase 1 weekend related elements -------------------------------------------

def num_wkend_days_idx_rule(M):
    """
    Construct index for num weekend days pattern

    :param M: Model
    :return: list of tuples of indexes
    """

    return [(i, w, t, e) for i in pyo.sequence(M.max_weekend_patterns)
            for w in M.WEEKS
            for t in M.TTYPES
            for e in pyo.sequence(2) if i <= M.num_weekend_patterns[e, t]]


model.num_wkend_days_idx = pyo.Set(dimen=4, ordered=True,
                                   initialize=num_wkend_days_idx_rule)


def A_num_wkend_days_init(M, i, w, t, e):
    """
    Initialize number of weekend days worked in the i'th pattern, for week k,
    tour type t and weekend type e.

    :param M: Model
    :param i: weekend worked pattern
    :param w: week
    :param t: tour type
    :param e: weekend type (1=Sun and Sat, 2=Fri and Sat)
    :return: number of weekend days worked (int)
    """

    if e == 1:
        return M.A_wkend_days[i, 1, w, t, e] + M.A_wkend_days[i, 7, w, t, e]

    else:
        return M.A_wkend_days[i, 6, w, t, e] + M.A_wkend_days[i, 7, w, t, e]


model.A_num_wkend_days = pyo.Param(model.num_wkend_days_idx,
                                   initialize=A_num_wkend_days_init)


def A_tot_wkend_days_idx_rule(M):
    """
    Index for number of days worked over scheduling cycle for each weekend pattern
    is (pattern, tour type, weekend type).

    :param M: Model
    :return: Parameter index (list of tuples)
    """
    return [(i, t, e) for i in pyo.sequence(M.max_weekend_patterns)
            for t in M.TTYPES
            for e in pyo.sequence(2) if i <= M.num_weekend_patterns[e, t]]


model.A_tot_wkend_days_idx = pyo.Set(dimen=3, ordered=True,
                                     initialize=A_tot_wkend_days_idx_rule)


def A_tot_wkend_days_init(M, i, t, e):
    """
    Initialize number of weekend days worked in the i'th pattern, tour type t,
    and weekend type e.

    :param M: Model
    :param i: weekend worked pattern
    :param t: tour type
    :param e: weekend type (1=Sun and Sat, 2=Fri and Sat)
    :return: number of weekend days worked over entire schedule cycle (int)
    """

    if e == 1:
        return sum(M.A_wkend_days[i, 1, w, t, e] + M.A_wkend_days[i, 7, w, t, e] for w in M.WEEKS)
    else:
        return sum(M.A_wkend_days[i, 6, w, t, e] + M.A_wkend_days[i, 7, w, t, e] for w in M.WEEKS)


model.A_tot_wkend_days = pyo.Param(model.A_tot_wkend_days_idx,
                                   initialize=A_tot_wkend_days_init)


def A_is_two_wkend_days_init(M, i, w, t, e):
    """
    Initialize indicator for each weekend worked pattern as to whether each week is
    sandwiched by consec weekend days worked.

    :param M: Model
    :param i: weekend worked pattern
    :param w: week
    :param t: tour type
    :param e: weekend type (1=Sun and Sat, 2=Fri and Sat)
    :return: 1 if True, 0 if False
    """

    if e == 1:
        if M.A_wkend_days[i, 1, w, t, e] == 1 and M.A_wkend_days[i, 7, w, t, e] == 1:
            return 1
        else:
            return 0
    else:
        if M.A_wkend_days[i, 6, w, t, e] == 1 and M.A_wkend_days[i, 7, w, t, e] == 1:
            return 1
        else:
            return 0


model.A_is_two_wkend_days = pyo.Param(model.num_wkend_days_idx,
                                      initialize=A_is_two_wkend_days_init)


def A_is_one_wkend_days_init(M, i, w, t, e):
    """
    Initialize indicator for each weekend worked pattern as to whether each week is
    sandwiched by consec weekend days worked.

    :param M: Model
    :param i: weekend worked pattern
    :param w: week
    :param t: tour type
    :param e: weekend type (1=Sun and Sat, 2=Fri and Sat)
    :return: 1 if True, 0 if False
    """

    if e == 1:
        if (M.A_wkend_days[i, 1, w, t, e] + M.A_wkend_days[i, 7, w, t, e]) == 1:
            return 1
        else:
            return 0

    else:
        if (M.A_wkend_days[i, 6, w, t, e] + M.A_wkend_days[i, 7, w, t, e]) == 1:
            return 1
        else:
            return 0


model.A_is_one_wkend_days = pyo.Param(model.num_wkend_days_idx,
                                      initialize=A_is_one_wkend_days_init)


def A_is_Sunday_init(M, i, w, t, e):
    """
    Initialize indicator for each weekend worked pattern as to whether it's
    a Sunday only pattern

    :param M: Model
    :param i: weekend worked pattern
    :param w: week
    :param t: tour type
    :param e: weekend type (1=Sun and Sat, 2=Fri and Sat)
    :return: 1 if True, 0 if False
    """

    if e == 1:
        if M.A_wkend_days[i, 1, w, t, e] == 1 and M.A_wkend_days[i, 7, w, t, e] == 0:
            return 1
        else:
            return 0

    else:
        # TODO - FS weekend not implemented
        if M.A_wkend_days[i, 1, w, t, e] == 1 and M.A_wkend_days[i, 7, w, t, e] == 0:
            return 1
        else:
            return 0


model.A_is_Sunday = pyo.Param(model.num_wkend_days_idx,
                              initialize=A_is_Sunday_init)


def A_is_Saturday_init(M, i, w, t, e):
    """
    Initialize indicator for each weekend worked pattern as to whether it's
    a Saturday only pattern

    :param M: Model
    :param i: weekend worked pattern
    :param w: week
    :param t: tour type
    :param e: weekend type (1=Sun and Sat, 2=Fri and Sat)
    :return: 1 if True, 0 if False
    """

    if e == 1:
        if M.A_wkend_days[i, 1, w, t, e] == 0 and M.A_wkend_days[i, 7, w, t, e] == 1:
            return 1
        else:
            return 0

    else:
        # TODO - FS weekend not implemented
        if M.A_wkend_days[i, 1, w, t, e] == 0 and M.A_wkend_days[i, 7, w, t, e] == 1:
            return 1
        else:
            return 0


model.A_is_Saturday = pyo.Param(model.num_wkend_days_idx,
                                initialize=A_is_Saturday_init)

# TODO: start windows ---------------------------------------------------------
# START WINDOWS - should these be tour type specific since allow start is tour type specific?
# - my inclination is yes since segregation by tour types makes constraint creation easier and
#   better isolates windows. Also, the more things that are tour type specific, the more
#   modeling flexibility we have.
#
# To start with, I made them week specific but not sure if they should be. Seems they should be
# period, day, tour type.
#
# Maybe best strategy is to do the variables and constraints first and then work backwards to
# define windows appropriately.
# -----------------------------------------------------------------------

# Start window width should be tour type specific
model.g_start_window_width = pyo.Param(default=0)  # Width of start-time windows

model.tt_start_window_width = pyo.Param(model.TTYPES, default=0)


# #/**** Beginning of each start window (in total periods from Sunday @ midnight)****/
# #pyo.Param b_window_wepoch{i in PERIODS,j in DAYS} := n_prds_per_day*(j-1)+i;
# #
# #/**** End of each start window (in total periods from Sunday @ midnight) ****/
# #pyo.Param e_window_wepoch{i in PERIODS,j in DAYS} :=git push origin --delete
# # ( if n_prds_per_day*(j-1)+i+width <= n_prds_per_day*n_days then
# #    n_prds_per_day*(j-1)+i+width
# #  else
# #    (n_prds_per_day*(j-1)+i+width )-n_prds_per_day*n_days);


# if model.g_start_window_width > 0:
# This version spans multiple weeks. I think I need to leave these week specific. The fact that allow_start is
# not week specific will force consistency across the weeks. I can always remove week if it becomes apparent
# that it's the right approach. For now, I don't recall exactly how this parameter gets used.
def b_window_epoch_init(M, i, j, w):
    return M.epoch[i, j, w]


model.b_window_epoch = pyo.Param(model.epoch_tuples,
                                 initialize=b_window_epoch_init)


def e_window_epoch_init(M, i, j, w):
    epoch = M.epoch[i, j, w] + M.g_start_window_width.value

    if epoch <= M.n_prds_per_cycle.value:
        e_window = epoch
    else:
        e_window = epoch - M.n_prds_per_cycle.value

    return e_window


model.e_window_epoch = pyo.Param(model.epoch_tuples,
                                 initialize=e_window_epoch_init)


# ##### PotentialGlobalStartWindow ######

# Index: PotentialGlobalStartWindow is defined for all (period,day,week) epoch_tuples.
# Defn:  PotentialGlobalStartWindow[i,j,w] contains all the bin triplets within g_start_window_width periods of
#    bin (i,j,w). The "Potential" indicates that many of these will be eliminated because [i,j,w] may
#    not be an allowable shift start time. Notice that these windows are NOT yet shift length and tour
#    type specific - hence the "Global" in the name.

# Single week model equivalent: WindowWepochs

# EXAMPLE: Let g_start_window_width=2. Then PotentialGlobalStartWindow[5,2,1] = [(5,2,1),(6,2,1),(7,2,1)]

def PotentialGlobalStartWindow_init(M, i, j, w):
    window_list = []
    for (l, m, n) in M.epoch_tuples:

        test_epoch = M.epoch[l, m, n]

        test1 = (test_epoch >= M.b_window_epoch[i, j, w])

        if M.b_window_epoch[i, j, w] <= M.e_window_epoch[i, j, w]:
            test2rhs = M.e_window_epoch[i, j, w]
        else:
            test2rhs = M.n_prds_per_cycle

        test2 = (test_epoch <= test2rhs)

        if M.b_window_epoch[i, j, w] <= M.e_window_epoch[i, j, w]:
            test3rhs = M.b_window_epoch[i, j, w]
        else:
            test3rhs = 1

        test3 = (test_epoch >= test3rhs)

        test4 = (test_epoch <= M.e_window_epoch[i, j, w])

        if (test1 and test2) or (test3 and test4):
            window_list.append((l, m, n))

    return window_list


model.PotentialGlobalStartWindow = pyo.Set(model.PERIODS, model.DAYS,
                                           model.WEEKS, dimen=3,
                                           ordered=True,
                                           initialize=PotentialGlobalStartWindow_init)


# ###/**** The set okWindowWepochs{i in PERIODS,j in DAYS,k in LENGTHS,t in TTYPES } creates shift
# ###length and tour type specific sets of windows that start in (i,j) pairs which have
# ###(i,j,k) as an allowable shift start time. ****/
# ###
# ###set okWindowWepochs{i in PERIODS,j in DAYS,k in LENGTHS,t in TTYPES}
# ### :={(l,m) in WindowWepochs[i,j]: allow_start[i,j,k,t]>0 and allow_start[l,m,k,t]>0};
#

# ##### PotentialStartWindow ######

# Index: PotentialStartWindow is defined for all (period,day,week,shift length,tour type).
# Defn:  PotentialStartWindow[i,j,k,t] contains all the bin triplets within g_start_window_width periods of
#    bin (i,j,w) and having (i,j,k,t) be an allowable start time and all of its elements being allowable start times.
#    The "Potential" indicates that some of these will be eliminated because they may
#    be subsets of other, nearby, start windows.

# Single week model equivalent: okWindowWepochs

# EXAMPLE: Let g_start_window_width=2 and only odd periods be allowable start time periods.
#    Then PotentialStartWindow[5,2,1,1,1] = [(5,2,1),(7,2,1)]

def PotentialStartWindow_idx_rule(M):
    return [(i, j, w, k, t) for i in M.PERIODS
            for j in M.DAYS
            for w in M.WEEKS
            for k in M.LENGTHS
            for t in M.activeTT]


model.PotentialStartWindow_idx = pyo.Set(dimen=5, ordered=True,
                                         initialize=PotentialStartWindow_idx_rule)


def PotentialStartWindow_init(M, i, j, w, k, t):
    return [(l, m, n) for (l, m, n) in M.PotentialGlobalStartWindow[i, j, w]
            if M.allow_start[i, j, k, t] > 0 and M.allow_start[l, m, k, t] > 0]


model.PotentialStartWindow = pyo.Set(model.PotentialStartWindow_idx,
                                     ordered=True,
                                     initialize=PotentialStartWindow_init,
                                     dimen=3)


# ##### okStartWindowRoots ######

# Index: okStartWindowRoots is defined for all (tour type, shift length) pairs such that the shift length
#        is allowed for the tour type.

# Defn:  okStartWindowRoots[t,k] contains all the bin triplets in PotentialStartWindow

# Single week model equivalent: ok_window_beginnings

# EXAMPLE: Let g_start_window_width=2 and only odd periods be allowable start time periods.
#    Then PotentialStartWindow[5,2,1,1,1] = [(5,2,1),(7,2,1)]


def okStartWindowRoots_idx_rule(M):
    return [(t, k) for k in M.LENGTHS
            for t in M.activeTT
            if k in M.tt_length_x[t]]


model.okStartWindowRoots_idx = pyo.Set(dimen=2, ordered=True,
                                       initialize=okStartWindowRoots_idx_rule)


#
# ##/**** ok_window_beginnings is the set of start windows in which there is
# ## at least one period in which a shift of length k can start
# ## and which are not subsets of some other window.
# ##
# ##*/
# #
# ##set okWindowBeginnings{t in okTTYPES, k in tt_length_x[t]} :=
# ##  setof{(p,q) in {PERIODS,DAYS}: (p,q,k,t) in ok_shifts and
# ##   forall{(i,j) in {PERIODS,DAYS} diff {(p,q)}:allow_start[i,j,k,t]>0 }
# ##    (not
# ##     ({(l,m) in okWindowWepochs[p,q,k,t]}
# ##      within {(n,o) in okWindowWepochs[i,j,k,t]}))  } (p,q);
# #
def okStartWindowRoots_init(M, t, k):
    window_list = []
    for (p, q, w) in M.epoch_tuples:
        # test1 = False

        # Create a list of all possible window beginnings and then eliminate those that
        # have window epochs that are subsets of other window epochs
        test1 = ((p, q, k, t) in M.okShiftTypes)
        if test1:
            window = []
            for (i, j, m) in M.PotentialStartWindow[p, q, w, k, t]:
                window.append((i, j, m))
            window_list.append(window)

    # Get rid of subsets
    window_list_copy = window_list[:]

    for s1 in window_list:
        for s2 in window_list:
            if set(s1) < set(s2) and s1 != s2:
                window_list_copy.remove(s1)
                break

    #   Now find the earliest bin in each element (list) of the list
    window_root_list = []
    for w in window_list_copy:
        onSaturday = False
        onSunday = False
        early = M.n_prds_per_cycle
        sat_early = M.n_prds_per_cycle
        late = 0
        for btup in w:
            if btup[1] == 7:
                onSaturday = True
            if btup[1] == 1:
                onSunday = True

            epoch_of_week = M.epoch[btup[0], btup[1], btup[2]]
            if epoch_of_week < early:
                early = epoch_of_week
            if epoch_of_week > late:
                late = epoch_of_week
            if btup[1] == 7 and epoch_of_week < sat_early:
                sat_early = epoch_of_week

        if onSaturday and onSunday:
            early = sat_early

        early_bin = epoch_to_tuple(M, early)
        window_root_list.append(early_bin)

    return window_root_list


model.okStartWindowRoots = pyo.Set(model.okStartWindowRoots_idx,
                                   dimen=3,
                                   ordered=True,
                                   initialize=okStartWindowRoots_init)


# #
# #
# ## CHAINS - NON-OVERLAPPING SEQUENCES OF (i,j) period,day PAIRS THAT
# ##          CAN BE ISOLATED FOR COORDINATING Shift AND TourTypeDayShift variables
# ##
# ## Let's wait on generalizing these for multiple weeks. Get model working
# ## for start window width of 0.
# ## -----------------------------------------------------------------------
# #
# #
# ##set bchain {t in okTTYPES, k in tt_length_x[t]} := setof{(w,j) in (okWindowBeginnings[t,k]):
# ##    forall{(p,q) in (okWindowBeginnings[t,k] diff {(w,j)})} (w,j) not in (okWindowWepochs[p,q,k,t])} (w,j) ;
# #
# #


def bchain_init(M, t, k):
    window_list = []
    if M.g_start_window_width > 0:
        for (i, j, w) in M.okStartWindowRoots[t, k]:
            is_contained = False
            for (p, q, r) in M.okStartWindowRoots[t, k]:
                #            for (p, q, r) in (M.okStartWindowRoots[t, k] - pyo.Set(initialize=[(i, j, w)])):
                if (i, j, w) in M.PotentialStartWindow[p, q, r, k, t]:
                    if (i, j, w) != (p, q, r):
                        is_contained = True
                        break
            if not is_contained:
                if (i, j, w) not in window_list:
                    window_list.append((i, j, w))

    return window_list


model.bchain = pyo.Set(model.okStartWindowRoots_idx, dimen=3,
                       ordered=True, initialize=bchain_init)


# #set echain {t in okTTYPES,k in LENGTHS,i in PERIODS,j in DAYS:
# # (i,j) in bchain[t,k]}
# # := setof{(w,x) in {PERIODS,DAYS}: (w,x) not in (bchain[t,k] diff {(i,j)}) and
# # (w,x) in okWindowBeginnings[t,k] and
# # forall{(p,q) in (okWindowBeginnings[t,k] diff {(w,x)})} (p,q) not in (okWindowWepochs[w,x,k,t]) and
# #     (period[w,x]>=period[i,j] and
# #     forall{(n,o) in bchain[t,k]: period[n,o]>period[i,j]} period[w,x]<
# #  period[n,o] or
# #  (period[w,x]<period[i,j] and
# #  forall{(n,o) in bchain[t,k]} (period[w,x]<
# #   period[n,o] and period[n,o]<=period[i,j]) ) )
# #  } (w,x) ;
#
#


def chain_idx_rule(M):
    return [(t, k, i, j, w) for t in M.activeTT
            for k in M.LENGTHS
            for i in M.PERIODS
            for j in M.DAYS
            for w in M.WEEKS
            if (t, k) in M.okStartWindowRoots_idx and (i, j, w) in M.bchain[t, k]]


model.chain_idx = pyo.Set(dimen=5, ordered=True, initialize=chain_idx_rule)


def echain_init(M, t, k, i, j, w):
    window_list = []
    # Compute epoch of (i,j,w)
    g_prd = M.epoch[i, j, w]
    prd = epoch_to_tuple(M, g_prd)

    steps = 1
    while steps <= M.g_start_window_width:
        g_prd_next = epoch_increment(M, g_prd, steps)
        prd_next = epoch_to_tuple(M, g_prd_next)
        if prd_next in M.PotentialStartWindow[prd[0], prd[1], prd[2], k, t]:
            g_prd = g_prd_next
            prd = epoch_to_tuple(M, g_prd)
            steps = 1
        else:
            steps = steps + 1

    # Once we leave the while loop, g_prd should correspond to echain. Need
    # to convert it back to a tuple from a epoch number.

    if prd not in window_list:
        window_list.append(prd)

    return window_list


model.echain = pyo.Set(model.chain_idx, ordered=True,
                       dimen=3, initialize=echain_init)


def n_links_init(M, t, k, i, j, w):
    # Compute epoch of (i,j,w)
    b_g_prd = M.epoch[i, j, w]
    e_prd = M.echain[t, k, i, j, w][1]
    e_g_prd = M.epoch[e_prd[0], e_prd[1], e_prd[2]]

    return mwts_shared.epoch_difference(M, b_g_prd, e_g_prd)


model.n_links = pyo.Param(model.chain_idx, initialize=n_links_init)


def chain_init(M, t, k, i, j, w):
    window_list = [(i, j, w)]
    # Compute epoch of (i,j,w)
    g_prd = M.epoch[i, j, w]
    # done = False

    steps = 1
    while steps <= M.n_links[t, k, i, j, w]:
        g_prd_next = epoch_increment(M, g_prd, steps)
        prd_next = epoch_to_tuple(M, g_prd_next)
        if prd_next in M.okStartWindowRoots[t, k]:
            window_list.append(prd_next)
        steps = steps + 1

    return window_list


model.chain = pyo.Set(model.chain_idx, ordered=True,
                      dimen=3, initialize=chain_init)


def link_idx_rule(M):
    # TODO Check the m index to see what the upper index limit should be
    return [(t, k, i, j, w, m) for t in M.activeTT
            for k in M.LENGTHS
            for i in M.PERIODS
            for j in M.DAYS
            for w in M.WEEKS
            for m in M.EPOCHS
            if (t, k) in M.okStartWindowRoots_idx and (i, j, w) in M.bchain[t, k]
            and m <= M.n_links[t, k, i, j, w]]


model.link_idx = pyo.Set(dimen=6, ordered=True, initialize=link_idx_rule)


def link_init(M, t, k, i, j, w, m):
    """
    Returns the m'th link (a (period,day,week) tuple) of the chain
    starting in period (i,j,w).
    """
    window_list = []
    g_prd = M.epoch[i, j, w]
    g_prd_next = epoch_increment(M, g_prd, m - 1)
    prd_next = epoch_to_tuple(M, g_prd_next)
    window_list.append(prd_next)

    return window_list


model.link = pyo.Set(model.link_idx, ordered=True,
                     dimen=3, initialize=link_init)


def linkspan_init(M, t, k, i, j, w, m):
    """
    Returns the start windows spanned by the m'th link (a (period,day,week) tuple)
    of the chain starting in period (i,j,w)
    """
    window_list = []
    for (p, d, q) in M.link[t, k, i, j, w, m]:
        for (a, b, c) in M.PotentialGlobalStartWindow[p, d, q]:
            if (a, b, c) not in window_list:
                window_list.append((a, b, c))
        if m > 1:
            for (a, b, c) in M.linkspan[t, k, i, j, w, m - 1]:
                if (a, b, c) not in window_list:
                    window_list.append((a, b, c))

    return window_list


model.linkspan = pyo.Set(model.link_idx, dimen=3,
                         ordered=True, initialize=linkspan_init)


# Variables ===================================================================

# Tour type
# Tour type daily
# Tour type daily shift
# Shift
# Weekend pattern
# Multiweek number of days worked pattern
# Multiweek patterns (weekend x number of days worked)
# Coverage components
# Understaffing

# Tour type Variables ---------------------------------------------------------

# TourType[i,t] -  Number of people assigned to tour type t in start window i


def okTourType_rule(M):
    """
    List of (window,tour type) tuples that are allowable.

    To be allowable,
    for every week there must be at least the minumum required number of
    days worked having an allowable shift (of any length allowed for that 
    tour type).

    :param M: Model
    """
    index_list = []
    for (i, t) in M.WINDOWS * M.activeTT:
        n_ok_weeks = 0
        for w in M.WEEKS:
            n_ok_days = 0
            for j in M.DAYS:
                for k in M.tt_length_x[t]:
                    if (i, j, w) in M.okStartWindowRoots[t, k]:
                        n_ok_days += 1
                        # Break out of the inner for since this day is ok
                        # for at least one shift length
                        break
            if n_ok_days >= sum(M.tt_shiftlen_min_dys_weeks[t, k, w] for k in M.tt_length_x[t]):
                n_ok_weeks += 1
        if n_ok_weeks == M.n_weeks:
            index_list.append((i, t))

    return index_list


model.okTourType = pyo.Set(dimen=2, initialize=okTourType_rule)


def TourType_idx_rule(M):
    """
    Index is (start window, tour type)

    :param M: Model
    :return: Variable index (list of tuples)
    """
    return [(i, t) for i in M.WINDOWS
            for t in M.activeTT
            if (i, t) in M.okTourType]


model.TourType_idx = pyo.Set(dimen=2, initialize=TourType_idx_rule)

model.TourType = pyo.Var(model.TourType_idx,
                         within=pyo.NonNegativeIntegers)

# Tour type daily variables ---------------------------------------------------

def okTourTypeDay_rule(M):
    """
    List of (window, tour type, day) tuples that are allowable.

    To be allowable,
    for every week and day there must be an allowable shift (of any length allowed for that
    tour type).

    :param M: Model
    """
    index_list = []
    for (i, t, j) in M.WINDOWS * M.activeTT * M.DAYS:
        n_ok_weeks = 0
        for w in M.WEEKS:
            n_ok_days = 0
            for k in M.tt_length_x[t]:
                    if (i, j, w) in M.okStartWindowRoots[t, k]:
                        n_ok_days += 1
                        # Break out of the inner for since this day is ok
                        # for at least one shift length
                        break
            if n_ok_days == 1:
                n_ok_weeks += 1
        if n_ok_weeks == M.n_weeks:
            index_list.append((i, t, j))

    return index_list


model.okTourTypeDay = pyo.Set(dimen=3, initialize=okTourTypeDay_rule)

# model.okTourTypeDay = model.okTourType * model.DAYS


# TourTypeDay[i,t,d] - Number of people assigned to tour type t in start window i
#                      and working day d in week w.

def TourTypeDay_idx_rule(M):
    """
    Index is (start window, tour type, day, week)

    :param M:
    :return: Constraint index (list of tuples)
    """

    return [(i, t, j, w) for i in M.WINDOWS
            for t in M.activeTT
            for j in M.DAYS
            for w in M.WEEKS
            if (i, t, j) in M.okTourTypeDay]


model.TourTypeDay_idx = pyo.Set(dimen=4, initialize=TourTypeDay_idx_rule)

model.TourTypeDay = pyo.Var(model.TourTypeDay_idx,
                            within=pyo.NonNegativeIntegers)


# Tour type daily shift variables ---------------------------------------------

def okTourTypeDayShift_rule(M):
    """
    List of (window, tour type, shift length, day) tuples that are allowable.

    To be allowable,
    for every week and day there must be an allowable shift a given length allowed for that
    tour type).

    :param M: Model
    """
    index_list = []
    for (i, t, j) in M.WINDOWS * M.activeTT * M.DAYS:
        # n_ok_weeks = 0
        for w in M.WEEKS:
            n_ok_days = 0
            for k in M.tt_length_x[t]:
                    if (i, j, w) in M.okStartWindowRoots[t, k] and (i, t, k, j) not in index_list:
                        index_list.append((i, t, k, j))

    return index_list


model.okTourTypeDayShift = pyo.Set(dimen=4, initialize=okTourTypeDayShift_rule)

def TourTypeDayShift_idx_rule(M):
    """
    Index is (start window, tour type, shift length, day, week)

    :param M: Model
    :return: Constraint index (list of tuples)
    """

    return [(i, t, k, j, w) for i in M.WINDOWS
            for t in M.activeTT
            for k in M.tt_length_x[t]
            for j in M.DAYS
            for w in M.WEEKS
            if (i, t, k, j) in M.okTourTypeDayShift]


model.TourTypeDayShift_idx = pyo.Set(dimen=5, initialize=TourTypeDayShift_idx_rule)

model.TourTypeDayShift = pyo.Var(model.TourTypeDayShift_idx,
                                 within=pyo.NonNegativeIntegers)

# Shift variables -------------------------------------------------------------

# Shift[i,j,w,k,t] = Number of shifts of length k starting in period i
# of day j in week w for a tour of type t

model.Shift = pyo.Var(model.okShifts,
                      within=pyo.NonNegativeIntegers)


# Weekend Days off variables --------------------------------------------------

# WeekendDaysWorked[d,i,t] = Number of people assigned days-off patterns d
# in start window i and of tour type t. These are convenience variables as
# the weekend patterns are part of the MultiWeekDaysWorked variables. Makes
# it easy to coordinate with Shift variables assigned to weekend days.


def weekend_days_worked_idx_rule(M):
    """
    Index is (start window, tour type, weekend pattern)

    :param M: Model
    :return: Constraint index (list of tuples)
    """

    index_list = []
    for (i, t) in M.okTourType:
        for p in pyo.sequence(M.max_weekend_patterns):
            weekend_type = M.weekend_type[i, t]
            if p <= M.num_weekend_patterns[weekend_type, t]:
                index_list.append((i, t, p))

    return index_list


model.weekend_days_worked_idx = pyo.Set(dimen=3,
                                        initialize=weekend_days_worked_idx_rule)

model.WeekendDaysWorked = pyo.Var(model.weekend_days_worked_idx,
                                  within=pyo.NonNegativeIntegers)


# Multiweek days worked variables --------------------------------------------------

# MultiWeekDaysWorked[i, t, p1, p2] = Number of people assigned to tour type t
# with mwdw pattern p1 and with weekend pattern p2 in start window i


def multiweekdaysworked_idx_rule(M):
    """
    Index is (start window, tour type, multiweek num days worked pattern,
    weekend pattern)

    :param M: Model
    :return: Constraint index (list of tuples)
    """

    index_list = []
    for (i, t) in M.okTourType:
        for p1 in pyo.sequence(M.max_mwdw_patterns):
            if p1 <= M.num_mwdw_patterns[t]:
                for p2 in pyo.sequence(M.max_weekend_patterns):
                    weekend_type = M.weekend_type[i, t]
                    if p2 <= M.num_weekend_patterns[weekend_type, t]:
                        index_list.append((i, t, p1, p2))

    return index_list


model.multiweekdaysworked_idx = \
    pyo.Set(dimen=4, initialize=multiweekdaysworked_idx_rule)

model.MultiWeekDaysWorked = \
    pyo.Var(model.multiweekdaysworked_idx,
            within=pyo.NonNegativeIntegers)

# Coverage and understaffing variables ----------------------------------------

# For computational convenience, we broke up the calculation of coverage in
# each period into four separate cases related
# different types of overlap (or lack of). See coverage constraints.

model.cov1 = pyo.Var(model.PERIODS, model.DAYS,
                     model.WEEKS, within=pyo.NonNegativeReals)

model.cov2 = pyo.Var(model.PERIODS, model.DAYS,
                     model.WEEKS, within=pyo.NonNegativeReals)

model.cov3 = pyo.Var(model.PERIODS, model.DAYS,
                     model.WEEKS, within=pyo.NonNegativeReals)

model.cov4 = pyo.Var(model.PERIODS, model.DAYS,
                     model.WEEKS, within=pyo.NonNegativeReals)

model.cov = pyo.Var(model.PERIODS, model.DAYS,
                    model.WEEKS, within=pyo.NonNegativeReals)


# Upper bound on tier 1 understaffing level
def under1_bounds(M, i, j, w):
    lb = 0.0
    ub = M.usb.value
    return lb, ub


model.under1 = pyo.Var(model.PERIODS, model.DAYS,
                       model.WEEKS, bounds=under1_bounds)


# Upper bound on tier 2 understaffing level
def under2_bounds(M, i, j, w):
    lb = 0.0
    ub = infinity
    return lb, ub


model.under2 = pyo.Var(model.PERIODS, model.DAYS,
                       model.WEEKS, bounds=under2_bounds)


# Objective function ==========================================================


def objective_rule(M):
    """
    Sum of labor costs and understaffing costs.

    :param M: Model
    :return: Objective function rule
    """
    obj1 = sum(M.Shift[i, j, w, k, t] * M.lengths[k] * M.tt_cost_multiplier[t]
               for (i, j, w, k, t) in M.okShifts)

    obj2 = sum(M.under1[i, j, w] * M.cu1.value + M.under2[i, j, w] * M.cu2.value
               for (i, j, w) in M.epoch_tuples)

    return obj1 + obj2


model.total_cost = pyo.Objective(rule=objective_rule, sense=pyo.minimize)


# Constraints =================================================================


# Budget constraints ----------------------------------------------------------

def max_labor_budget_rule(M):
    """
    Put upper bound on total labor cost.
    
    Using individual shift variables multiplied by their length and a tour
    type specific cost multiplier. Could easily generalize this to make costs
    be complex function of time of day, day of week, shift length, tour type,
    or some combination of.
    
    :param M: Model
    :return: Constraint rule
    """
    return sum(M.Shift[i, j, w, k, t] * M.lengths[k] * M.tt_cost_multiplier[t]
               for (i, j, w, k, t) in M.okShifts) <= M.labor_budget


model.max_labor_budget = pyo.Constraint(rule=max_labor_budget_rule)


# Coverage constraints --------------------------------------------------------

# Breaking them up into four different constraints, one for each case
# in terms of handling end of day horizon wrapping


def coverage1_rule(M, i, j, w):
    """
    Shift does not span beginning of any day.

    :param M: Model
    :param i: period
    :param j: day
    :param w: week
    :return: Constraint rule
    """

    return sum(M.Shift[(i - p), j, w, l, t]
               for l in M.LENGTHS
               for t in M.activeTT
               for p in range(0, M.lengths[l])
               if (i - p) > 0 and
               M.allow_start[(i - p), j, l, t] > 0) == M.cov1[i, j, w]


model.coverage1 = pyo.Constraint(model.PERIODS, model.DAYS,
                                 model.WEEKS,
                                 rule=coverage1_rule)


def coverage2_rule(M, i, j, w):
    """
    Shift spans beginning of some day other than the first day.

    :param M: Model
    :param i: period
    :param j: day
    :param w: week
    :return: Constraint rule
    """

    return sum(M.Shift[M.n_prds_per_day.value + (i - p), j - 1, w, l, t]
               for l in M.LENGTHS
               for t in M.activeTT
               for p in range(0, M.lengths[l])
               if (i - p) <= 0 and j > 1 and
               M.allow_start[M.n_prds_per_day.value + (i - p), j - 1, l, t] > 0) == M.cov2[i, j, w]


model.coverage2 = pyo.Constraint(model.PERIODS, model.DAYS,
                                 model.WEEKS,
                                 rule=coverage2_rule)


def coverage3_rule(M, i, j, w):
    """
    Shift does span beginning of the first day of week but not the
    first week (i.e. doesn't span beginning of cycle).

    :param M: Model
    :param i: period
    :param j: day
    :param w: week
    :return: Constraint rule
    """

    return sum(M.Shift[M.n_prds_per_day.value + (i - p), M.n_days_per_week.value, w - 1, l, t]
               for l in M.LENGTHS
               for t in M.activeTT
               for p in range(0, M.lengths[l])
               if (i - p) <= 0 and j == 1 and w > 1 and
               M.allow_start[M.n_prds_per_day.value + (i - p), M.n_days_per_week.value, l, t] > 0
               ) == M.cov3[i, j, w]


model.coverage3 = pyo.Constraint(model.PERIODS, model.DAYS, model.WEEKS,
                                 rule=coverage3_rule)


def coverage4_rule(M, i, j, w):
    """
    Shift spans beginning of the scheduling cycle (day 1, week 1).

    :param M: Model
    :param i: period
    :param j: day
    :param w: week
    :return: Constraint rule
    """

    return sum(
        M.Shift[M.n_prds_per_day.value + (i - p), M.n_days_per_week.value, M.n_weeks.value, l, t]
        for l in M.LENGTHS
        for t in M.activeTT
        for p in range(0, M.lengths[l])
        if (i - p) <= 0 and j == 1 and w == 1 and
        M.allow_start[M.n_prds_per_day.value + (i - p), M.n_days_per_week.value, l, t] > 0
    ) == M.cov4[i, j, w]


model.coverage4 = pyo.Constraint(model.PERIODS, model.DAYS, model.WEEKS,
                                 rule=coverage4_rule)


def tot_coverage_rule(M, i, j, w):
    """
    Total staffing coverage

    :param M: Model
    :param i: period
    :param j: day
    :param w: week
    :return: Constraint rule
    """
    return M.cov1[i, j, w] + M.cov2[i, j, w] + M.cov3[i, j, w] + M.cov4[i, j, w] == M.cov[i, j, w]


model.tot_coverage = pyo.Constraint(model.PERIODS,
                                    model.DAYS,
                                    model.WEEKS, rule=tot_coverage_rule)


def coverage_rule(M, i, j, w):
    """
    Total staffing coverage constraint

    :param M: Model
    :param i: period
    :param j: day
    :param w: week
    :return: Constraint rule
    """
    return M.cov[i, j, w] + M.under1[i, j, w] + M.under2[i, j, w] >= M.dmd_staff[i, j, w]


model.coverage = pyo.Constraint(model.PERIODS, model.DAYS,
                                model.WEEKS, rule=coverage_rule)


def minstaff_rule(M, i, j, w):
    """
    Minimum staffing levels required

    :param M: Model
    :param i: period
    :param j: day
    :param w: week
    :return: Constraint rule
    """

    # Note no understaffing on lhs
    return M.cov[i, j, w] >= M.min_staff[i, j, w]


model.minstaff = pyo.Constraint(model.PERIODS, model.DAYS,
                                model.WEEKS, rule=minstaff_rule)


# Bounds on tour type variables -----------------------------------------------
def TourType_LB_rule(M, t):
    """
    Lower bound on tour type

    :param M: Model
    :param t: tour type
    :return: Constraint rule
    """

    return sum(M.TourType[i, t] for (i, s) in M.okTourType if s == t) >= M.tt_lb[t]


model.TourType_LB_con = pyo.Constraint(model.activeTT, rule=TourType_LB_rule)


def TourType_UB_rule(M, t):
    """
    Upper bound on tour type

    :param M: Model
    :param t: tour type
    :return: Constraint rule
    """

    return sum(M.TourType[i, t] for (i, s) in M.okTourType if s == t) <= M.tt_ub[t]


model.TourType_UB_con = pyo.Constraint(model.activeTT, rule=TourType_UB_rule)


# Each tour type variable must get a mwdw pattern assigned to it
def MWDW_total_idx_rule(M):
    """
    The index for the MWDW constraints is (window, tour type)
    in the set okTourType.

    :param M: Model
    :return: Constraint index rule
    """

    return [(i, t) for (i, t) in M.okTourType]


def MWDW_total_rule(M, i, t):
    """
    Each tour type variable must get a multi-week pattern assigned to it.

    :param M: Model
    :param i: start window
    :param t: tour type
    :return: Constraint rule
    """

    weekend_type = M.weekend_type[i, t]
    return sum(M.MultiWeekDaysWorked[i, t, p1, p2]
               for p1 in pyo.sequence(M.num_mwdw_patterns[t])
               for p2 in pyo.sequence(M.num_weekend_patterns[weekend_type, t])
               ) == M.TourType[i, t]


model.MWDW_total_idx = pyo.Set(dimen=2, initialize=MWDW_total_idx_rule)

model.MWDW_total_con = pyo.Constraint(model.MWDW_total_idx,
                                      rule=MWDW_total_rule)


# Coordinate WeekendDaysWorked variables and multi-week pattern variables.
def MWDW_weekend_integration_rule(M, i, t, p):
    """
    Each weekend variable must be consistent with multi-week pattern variables.

    Could actually eliminate the weekend variables but convenient for use in constraints.

    :param M: Model
    :param i: start window
    :param t: tour type
    :param p: weekend pattern
    :return: Constraint rule
    """

    return sum(M.MultiWeekDaysWorked[i, t, p1, p]
               for p1 in pyo.sequence(
        M.num_mwdw_patterns[t])) == M.WeekendDaysWorked[i, t, p]


model.MWDW_weekend_integration_con = pyo.Constraint(
    model.weekend_days_worked_idx, rule=MWDW_weekend_integration_rule)


# Integrate WeekendDaysWorked variables with TourTypeDay variables
def weekend_integration_1_SS_idx_rule(M):
    """
    The index for type 1 weekend-tour type integration constraints is
    (day, week, window, tour type).

    :param M: Model
    :return: Constraint index rule
    """

    return [(j, w, i, t) for j in M.DAYS
            for w in M.WEEKS
            for i in M.WINDOWS
            for t in M.activeTT
            if j in M.weekend[i, t]
            and 1 in M.weekend[i, t]
            and 7 in M.weekend[i, t]
            and (i, t, j) in M.okTourTypeDay]


model.weekend_integration_1_SS_idx = pyo.Set(
    dimen=4, initialize=weekend_integration_1_SS_idx_rule)


def weekend_integration_1_SS_rule(M, j, w, i, t):
    """
    The number of people working on Sun or Sat as part of a
    type 1 weekend worked pattern must be equal to the
    number of people working on Sun or Sat as specified
    by TourTypeDay variables.

    :param M: Model
    :param j: day
    :param w: week
    :param i: window
    :param t: tour type
    :return: Constraint rule
    """

    return sum(M.A_wkend_days[p, j, w, t, 1] * M.WeekendDaysWorked[i, t, p]
               for p in pyo.sequence(
        M.num_weekend_patterns[1, t])) == M.TourTypeDay[i, t, j, w]


model.weekend_integration_1_SS_con = pyo.Constraint(
    model.weekend_integration_1_SS_idx,
    rule=weekend_integration_1_SS_rule)


def weekend_integration_2_FS_idx_rule(M):
    """
    The index for type 2 weekend-tour type integration constraints is
    (day, week, window, tour type).

    :param M: Model
    :return: Constraint index rule
    """
    return [(j, w, i, t) for j in M.DAYS
            for w in M.WEEKS
            for i in M.WINDOWS
            for t in M.activeTT
            if j in M.weekend[i, t] and 6 in M.weekend[i, t] and 7 in M.weekend[i, t]
            and (i, t, j) in M.okTourTypeDay]


model.weekend_integration_2_FS_idx = pyo.Set(
    dimen=4, initialize=weekend_integration_2_FS_idx_rule)


def weekend_integration_2_FS_rule(M, j, w, i, t):
    """
    The number of people working on Fri or Sat as part of a
    type 2 weekend worked pattern must be equal to the
    number of people working on Sun or Sat as specified
    by TourTypeDay variables.

    :param M: Model
    :param j: day
    :param w: week
    :param i: window
    :param t: tour type
    :return: Constraint rule
    """
    return sum(M.A_wkend_days[p, j, w, t, 2] * M.WeekendDaysWorked[i, t, p]
               for p in pyo.sequence(
        M.num_weekend_patterns[2, t])) == M.TourTypeDay[i, t, j, w]


model.weekend_integration_2_FS_con = pyo.Constraint(
    model.weekend_integration_2_FS_idx,
    rule=weekend_integration_2_FS_rule)


# Integrate shift, tour type days, and tour type variables
def shift_TTD_dailyconservation_idx_rule(M):
    """
    Index is (window, day, week, tour type).

    The index set gets used in the TTD_TTDS_con and TTD_TT_UB constraints.

    :param M:
    :return: Constraint index rule
    """
    return [(i, j, w, t) for i in M.WINDOWS
            for j in M.DAYS
            for w in M.WEEKS
            for t in M.activeTT
            if (i, t, j) in M.okTourTypeDay]


model.shift_TTD_dailyconservation_idx = pyo.Set(
    dimen=4, initialize=shift_TTD_dailyconservation_idx_rule)


def TTD_TT_UB_rule(M, i, j, w, t):
    """
    Every day of every week for each (window, ttype), there can
    be no more people scheduled (TourTypeDay) than
    number of people assigned to TourType[window, ttype].

    :param M: Model
    :param i: window
    :param j: day
    :param w: week
    :param t: tour type
    :return: Constraint rule
    """
    return M.TourTypeDay[i, t, j, w] <= M.TourType[i, t]


model.TTD_TT_UB = \
    pyo.Constraint(model.shift_TTD_dailyconservation_idx,
                   rule=TTD_TT_UB_rule)


def TTD_TTDS_rule(M, i, j, w, t):
    """
    Coordinate TourTypeDayShift and TourTypeDay variables for each day of each week

    :param M: Model
    :param i: window
    :param j: day
    :param w: week
    :param t: tour type
    :return: Constraint rule
    """

    return sum(M.TourTypeDayShift[i, t, k, j, w]
               for k in M.tt_length_x[t] if (i, t, k, j) in M.okTourTypeDayShift) == M.TourTypeDay[i, t, j, w]


model.TTD_TTDS_con = pyo.Constraint(
    model.shift_TTD_dailyconservation_idx, rule=TTD_TTDS_rule)


# Now that mwdw vars added, we can directly determine TTD sum over each
# week instead of using the LB, UBs.
def TTD_MWDW_idx_rule(M):
    """
    Index is (window, tour type, week).

    :param M: Model
    :return: Constraint index rule
    """

    index_list = []
    for t in M.activeTT:
        n_patterns = M.num_mwdw_patterns[t]
        if n_patterns > 0:
            for i in M.WINDOWS:
                if (i, t) in M.okTourType:
                    for w in M.WEEKS:
                        index_list.append((i, t, w))
    return index_list


def TTD_MWDW_rule(M, i, t, w):
    """
    Coordinate TTD and MWDW vars for each week by summing over the days.

    :param i: start window
    :param t: tour type
    :param w: week
    :return: Constraint rule
    """

    weekend_type = M.weekend_type[i, t]
    return sum(M.TourTypeDay[i, t, j, w] for j in M.DAYS) == sum(
        M.MultiWeekDaysWorked[i, t, p1, p2] * M.A_mwdw[t, p1, w]
        for p1 in pyo.sequence(M.num_mwdw_patterns[t])
        for p2 in pyo.sequence(M.num_weekend_patterns[weekend_type, t]))


model.TTD_MWDW_idx = pyo.Set(dimen=3, initialize=TTD_MWDW_idx_rule)

model.TTD_MWDW_con = pyo.Constraint(model.TTD_MWDW_idx,
                                    rule=TTD_MWDW_rule)


# Shiftlen versions of the TTD_TT constraints using the TTDS vars.
# These are needed as they are shift length specific.

def TTDS_TT_shiftlen_weeklyconservation_idx_rule(M):
    """
    Index is (window, tour type, shift length, week) if more than one shift length
    in tour type.

    :param M: Model
    :return: Constraint index rule
    """
    index = []
    for (i, t) in M.okTourType:
        for k in M.tt_length_x[t]:
            if len(M.tt_length_x[t]) > 1:
                for w in M.WEEKS:
                    for d in M.DAYS:
                        if (i, t, k, d) in M.okTourTypeDayShift and (i, t, k, w) not in index:
                            index.append((i, t, k, w))

    return index


model.TTDS_TT_shiftlen_weeklyconservation_idx = pyo.Set(
    dimen=4, initialize=TTDS_TT_shiftlen_weeklyconservation_idx_rule)


def TTDS_TT_shiftlen_weeklyconservation_LB_rule(M, i, t, k, w):
    """
    Coordinate TTDS and TT vars for num days worked
    each week - lower bounds.

    :param M: Model
    :param i: start window
    :param t: tour type
    :param k: shift length
    :param w: week
    :return: Constraint rule
    """

    return sum(M.TourTypeDayShift[i, t, k, d, w]
               for d in M.DAYS if (i, t, k, d) in M.okTourTypeDayShift) >= \
           M.TourType[i, t] * M.tt_shiftlen_min_dys_weeks[t, k, w]


model.TTDS_TT_shiftlen_weeklyconservation_LB = \
    pyo.Constraint(
        model.TTDS_TT_shiftlen_weeklyconservation_idx,
        rule=TTDS_TT_shiftlen_weeklyconservation_LB_rule)


def TTDS_TT_shiftlen_weeklyconservation_UB_rule(M, i, t, k, w):
    """
    Coordinate TTDS and TT vars for num days worked
    each week - upper bounds.

    :param M: Model
    :param i: start window
    :param t: tour type
    :param k: shift length
    :param w: week
    :return: Constraint rule
    """

    return sum(M.TourTypeDayShift[i, t, k, d, w]
               for d in M.DAYS if (i, t, k, d) in M.okTourTypeDayShift) <= \
           M.TourType[i, t] * M.tt_shiftlen_max_dys_weeks[t, k, w]


model.TTDS_TT_shiftlen_weeklyconservation_UB_con = \
    pyo.Constraint(
        model.TTDS_TT_shiftlen_weeklyconservation_idx,
        rule=TTDS_TT_shiftlen_weeklyconservation_UB_rule)


def TTDS_TT_shiftlen_cumul_weeklyconservation_LB_rule(M, i, t, k, w):
    """
    Coordinate TTDS and TT vars for num days worked
    each week - cumulative lower bounds.

    :param M: Model
    :param i: start window
    :param t: tour type
    :param k: shift length
    :param w: week
    :return: Constraint rule
    """

    return sum(M.TourTypeDayShift[i, t, k, d, z]
               for d in M.DAYS for z in pyo.sequence(w) if (i, t, k, d) in M.okTourTypeDayShift) >= \
           M.TourType[i, t] * M.tt_shiftlen_min_cumul_dys_weeks[t, k, w]


model.TTDS_TT_shiftlen_cumul_weeklyconservation_LB = \
    pyo.Constraint(model.TTDS_TT_shiftlen_weeklyconservation_idx,
                   rule=TTDS_TT_shiftlen_cumul_weeklyconservation_LB_rule)


def TTDS_TT_shiftlen_cumul_weeklyconservation_UB_rule(M, i, t, k, w):
    """
    Coordinate TTDS and TT vars for num days worked
    each week - cumulative upper bounds.

    :param M: Model
    :param i: start window
    :param t: tour type
    :param k: shift length
    :param w: week
    :return: Constraint rule
    """
    return sum(M.TourTypeDayShift[i, t, k, d, z]
               for d in M.DAYS for z in pyo.sequence(w) if (i, t, k, d) in M.okTourTypeDayShift) <= \
           M.TourType[i, t] * M.tt_shiftlen_max_cumul_dys_weeks[t, k, w]


model.TTDS_TT_shiftlen_cumul_weeklyconservation_UB = \
    pyo.Constraint(model.TTDS_TT_shiftlen_weeklyconservation_idx,
                   rule=TTDS_TT_shiftlen_cumul_weeklyconservation_UB_rule)


def prds_worked_weekly_idx_rule(M):
    """
    Index is (window, tour type, week) if more than one shift length
    in tour type.

    :param M: Model
    :return: Constraint index rule
    """
    return [(i, t, w) for (i, t) in M.okTourType
            for w in M.WEEKS
            if len(M.tt_length_x[t]) > 1]


model.prds_worked_weekly_idx = pyo.Set(
    dimen=3, initialize=prds_worked_weekly_idx_rule)


def prds_worked_weekly_LB_rule(M, i, t, w):
    """
    Coordinate TTDS and TT vars for num periods worked
    each week - lower bounds.

    :param M: Model
    :param i: start window
    :param t: tour type
    :param w: week
    :return: Constraint rule
    """

    return sum(M.TourTypeDayShift[i, t, k, d, w] * M.lengths[k]
               for d in M.DAYS for k in M.tt_length_x[t] if
               (i, t, k, d, w) in M.TourTypeDayShift_idx) >= \
           M.TourType[i, t] * M.tt_min_prds_weeks[t, w]


model.prds_worked_weekly_LB = \
    pyo.Constraint(model.prds_worked_weekly_idx, rule=prds_worked_weekly_LB_rule)


def prds_worked_weekly_UB_rule(M, i, t, w):
    """
    Coordinate TTDS and TT vars for num periods worked
    each week - upper bounds.

    :param M: Model
    :param i: start window
    :param t: tour type
    :param w: week
    :return: Constraint rule
    """

    return sum(M.TourTypeDayShift[i, t, k, d, w] * M.lengths[k]
               for d in M.DAYS for k in M.tt_length_x[t]
               if (i, t, k, d, w) in M.TourTypeDayShift_idx) <= \
           M.TourType[i, t] * M.tt_max_prds_weeks[t, w]


model.prds_worked_weekly_UB = pyo.Constraint(
    model.prds_worked_weekly_idx, rule=prds_worked_weekly_UB_rule)


# # Cumulative versions of the above 2 constraints
def prds_worked_cumul_weekly_LB_rule(M, i, t, w):
    """
    Coordinate TTDS and TT vars for num periods worked
    each week - cumulative lower bounds.

    :param M: Model
    :param i: start window
    :param t: tour type
    :param w: week
    :return: Constraint rule
    """

    return sum(M.TourTypeDayShift[i, t, k, d, z] * M.lengths[k]
               for d in M.DAYS for k in M.tt_length_x[t]
               for z in pyo.sequence(w)
               if (i, t, k, d, z) in M.TourTypeDayShift_idx) >= \
           M.TourType[i, t] * M.tt_min_cumul_prds_weeks[t, w]


model.prds_worked_cumul_weekly_LB = \
    pyo.Constraint(model.prds_worked_weekly_idx,
                   rule=prds_worked_cumul_weekly_LB_rule)


def prds_worked_cumul_weekly_UB_rule(M, i, t, w):
    """
    Coordinate TTDS and TT vars for num periods worked each
    week - cumulative upper bounds.

    :param M: Model
    :param i: start window
    :param t: tour type
    :param w: week
    :return: Constraint rule
    """
    return sum(M.TourTypeDayShift[i, t, k, d, z] * M.lengths[k]
               for d in M.DAYS for k in M.tt_length_x[t]
               for z in pyo.sequence(w)
               if (i, t, k, d, z) in M.TourTypeDayShift_idx) <= \
           M.TourType[i, t] * M.tt_max_cumul_prds_weeks[t, w]


model.prds_worked_cumul_weekly_UB = \
    pyo.Constraint(model.prds_worked_weekly_idx,
                   rule=prds_worked_cumul_weekly_UB_rule)


# Weekend subset constraints --------------------------------------------------
# Need constraints to prevent cases such as the following. Consider two people
# with same tour type working 5 days per week. Assume that consecutive weekends
# are allowed and that the two weekend patterns for week 1
# are:  1 0 0 0 0 0 1 and 0 0 0 0 0 0 0.
# Now consider the following TourTypeDay solution:

#  1 2 0 2 2 2 1

# The employee with the 1 0 0 0 0 0 1 pattern would be forced to work 6 days.


def weekend_subsets_5_4_idx_rule(M):
    """
    TODO: Write me
    Only implemented for weekend type 1 (Sat, Sun)
    :param M: Model
    :return:
    """

    index_list = []

    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_max_dys_weeks[t, w] == 5 and e == 1:
                    index_list.append((i, t, w, e, 2, 3, 4, 5))
                    index_list.append((i, t, w, e, 2, 3, 4, 6))
                    index_list.append((i, t, w, e, 2, 3, 5, 6))
                    index_list.append((i, t, w, e, 2, 4, 5, 6))
                    index_list.append((i, t, w, e, 3, 4, 5, 6))

    return index_list


model.weekend_subsets_5_4_idx = pyo.Set(dimen=8, initialize=weekend_subsets_5_4_idx_rule)


def weekend_subsets_5_4_rule(M, i, t, w, e, d1, d2, d3, d4):
    """
    TODO: Write me

    :param M: Model
    :param i:
    :param t:
    :param w:
    :param e:
    :param d1:
    :param d2:
    :param d3:
    :param d4:
    :return:
    """

    days_subset = [d1, d2, d3, d4]
    return sum(M.TourTypeDay[i, t, d, w] for d in days_subset) <= \
           sum(M.MultiWeekDaysWorked[i, t, p1, p2] *
               min(len(days_subset),
                   (M.A_mwdw[t, p1, w] - M.A_num_wkend_days[p2, w, t, e]))
               for p1 in pyo.sequence(M.num_mwdw_patterns[t])
               for p2 in pyo.sequence(M.num_weekend_patterns[e, t]))


model.weekend_subsets_5_4_con2 = pyo.Constraint(
    model.weekend_subsets_5_4_idx,
    rule=weekend_subsets_5_4_rule)


def weekend_subsets_4_3_idx_rule(M):
    """
    TODO: Write me

    :param M: Model
    :return:
    """

    index_list = []
    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] <= 4 <= M.tt_max_dys_weeks[t, w] and e == 1:
                    index_list.append((i, t, w, e, 2, 3, 4))
                    index_list.append((i, t, w, e, 2, 3, 5))
                    index_list.append((i, t, w, e, 2, 3, 6))
                    index_list.append((i, t, w, e, 2, 4, 5))
                    index_list.append((i, t, w, e, 2, 4, 6))
                    index_list.append((i, t, w, e, 2, 5, 6))

                    index_list.append((i, t, w, e, 3, 4, 5))
                    index_list.append((i, t, w, e, 3, 4, 6))
                    index_list.append((i, t, w, e, 3, 5, 6))

                    index_list.append((i, t, w, e, 4, 5, 6))

    return index_list


model.weekend_subsets_4_3_idx = pyo.Set(
    dimen=7,
    initialize=weekend_subsets_4_3_idx_rule)


def weekend_subsets_4_3_rule(M, i, t, w, e, d1, d2, d3):
    """
    TODO: Write me

    :param M: Model
    :param i:
    :param t:
    :param w:
    :param e:
    :param d1:
    :param d2:
    :param d3:
    :return:
    """

    days_subset = [d1, d2, d3]
    return sum(M.TourTypeDay[i, t, d, w] for d in days_subset) <= \
           sum(M.MultiWeekDaysWorked[i, t, p1, p2] *
               min(len(days_subset),
                   (M.A_mwdw[t, p1, w] - M.A_num_wkend_days[p2, w, t, e]))
               for p1 in pyo.sequence(M.num_mwdw_patterns[t])
               for p2 in pyo.sequence(M.num_weekend_patterns[e, t]))


model.weekend_subsets_4_3_con2 = pyo.Constraint(
    model.weekend_subsets_4_3_idx,
    rule=weekend_subsets_4_3_rule)


def weekend_subsets_3_2_idx_rule(M):
    """
    TODO: Write me

    :param M: Model
    :return:
    """

    index_list = []

    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] <= 3 <= M.tt_max_dys_weeks[t, w] and e == 1:
                    index_list.append((i, t, w, e, 2, 3))
                    index_list.append((i, t, w, e, 2, 4))
                    index_list.append((i, t, w, e, 2, 5))
                    index_list.append((i, t, w, e, 2, 6))

                    index_list.append((i, t, w, e, 3, 4))
                    index_list.append((i, t, w, e, 3, 5))
                    index_list.append((i, t, w, e, 3, 6))

                    index_list.append((i, t, w, e, 4, 5))
                    index_list.append((i, t, w, e, 4, 6))
                    index_list.append((i, t, w, e, 5, 6))

    return index_list


model.weekend_subsets_3_2_idx = pyo.Set(dimen=6, initialize=weekend_subsets_3_2_idx_rule)


def weekend_subsets_3_2_rule(M, i, t, w, e, d1, d2):
    """
    TODO: Write me

    :param M: Model
    :param i:
    :param t:
    :param w:
    :param e:
    :param d1:
    :param d2:
    :return:
    """

    days_subset = [d1, d2]
    return sum(M.TourTypeDay[i, t, d, w] for d in days_subset) <= \
           sum(M.MultiWeekDaysWorked[i, t, p1, p2] *
               min(len(days_subset),
                   (M.A_mwdw[t, p1, w] - M.A_num_wkend_days[p2, w, t, e]))
               for p1 in pyo.sequence(M.num_mwdw_patterns[t])
               for p2 in pyo.sequence(M.num_weekend_patterns[e, t]))


model.weekend_subsets_3_2_con2 = pyo.Constraint(
    model.weekend_subsets_3_2_idx,
    rule=weekend_subsets_3_2_rule)


def weekend_subsets_2_1_idx_rule(M):
    """
    TODO: Write me

    :param M: Model
    :return:
    """

    index_list = []

    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] <= 2 <= M.tt_max_dys_weeks[t, w] and e == 1:
                    index_list.append((i, t, w, e, 2))
                    index_list.append((i, t, w, e, 3))
                    index_list.append((i, t, w, e, 4))
                    index_list.append((i, t, w, e, 5))
                    index_list.append((i, t, w, e, 6))

    return index_list


model.weekend_subsets_2_1_idx = pyo.Set(
    dimen=5,
    initialize=weekend_subsets_2_1_idx_rule)


def weekend_subsets_2_1_rule(M, i, t, w, e, d1):
    """
    TODO: Write me

    :param M: Model
    :param i:
    :param t:
    :param w:
    :param e:
    :param d1:
    :return:
    """

    days_subset = [d1]
    return sum(M.TourTypeDay[i, t, d, w] for d in days_subset) <= \
           sum(M.MultiWeekDaysWorked[i, t, p1, p2] *
               min(len(days_subset),
                   (M.A_mwdw[t, p1, w] - M.A_num_wkend_days[p2, w, t, e]))
               for p1 in pyo.sequence(M.num_mwdw_patterns[t])
               for p2 in pyo.sequence(M.num_weekend_patterns[e, t]))


model.weekend_subsets_2_1_con2 = pyo.Constraint(
    model.weekend_subsets_2_1_idx,
    rule=weekend_subsets_2_1_rule)


def max_ptfrac_rule(M):
    """
    Max fraction of part-time periods covered.

    :param M: Model
    :return: Constraint rule
    """

    return sum(M.Shift[i, j, w, k, t] * M.lengths[k]
               for (i, j, w, k, t) in M.okShifts if M.tt_parttime[t] > 0) <= \
           M.max_parttime_frac.value * sum(M.Shift[i, j, w, k, t] * M.lengths[k]
                                           for (i, j, w, k, t) in M.okShifts)


model.max_ptfrac_con = pyo.Constraint(rule=max_ptfrac_rule)


# Chains for intra-tour start time flexibility --------------------------------

# The follow proxy constraints are only case of no
# intra-tour start time flexibility
# Coordinate Shift variables with TourTypeDayShift for start window width=0.
# In this case, shift start time periods are same as start time windows

def chains_tot_proxy1_rule(M, w, t, k, i, j):
    """
    Coordinate shift and tour type day shift variables.

    :param M: Model
    :param w: week
    :param t: tour type
    :param k: shift length
    :param i: window
    :param j: day
    :return: Constraint rule
    """

    return M.Shift[i, j, w, k, t] == M.TourTypeDayShift[i, t, k, j, w]


def chains_tot_proxy1_idx_rule(M):
    """
    Index is (week, tour type, shift length, window, day)

    :param M: Model
    :return: Constraint index rule
    """

    if M.g_start_window_width == 0:
        return [(w, t, k, i, j) for w in M.WEEKS
                for t in M.activeTT
                for k in M.tt_length_x[t]
                for i in M.PERIODS
                for j in M.DAYS
                if (i, j, w, k, t) in M.okShifts
                if (i, t, k, j, w) in M.TourTypeDayShift_idx]
    else:
        return []


model.chains_tot_proxy1_idx = pyo.Set(
    dimen=5,
    initialize=chains_tot_proxy1_idx_rule)

model.chains_tot_proxy1_con = pyo.Constraint(
    model.chains_tot_proxy1_idx,
    rule=chains_tot_proxy1_rule)


def chains_tot_proxy2_rule(M, w, t, k, i, j):
    """
    Tour type day shift variables must respect allowable start times.

    :param M: Model
    :param w: week
    :param t: tour type
    :param k: shift length
    :param i: period/window
    :param j: day
    :return: Constraint rule
    """

    return M.TourTypeDayShift[i, t, k, j, w] == 0


def chains_tot_proxy2_idx_rule(M):
    """
    Index (tour type, shift length, period/window, day)
    :param M:
    :return: Constraint index rule
    """

    if M.g_start_window_width == 0:
        return [(w, t, k, i, j) for w in M.WEEKS
                for t in M.activeTT
                for k in M.tt_length_x[t]
                for i in M.PERIODS
                for j in M.DAYS
                if (i, t, k, j, w) in M.TourTypeDayShift_idx
                if (i, t, j) in M.okTourTypeDay and M.allow_start[i, j, k, t] == 0]
    else:
        return []


model.chains_tot_proxy2_idx = pyo.Set(
    dimen=5,
    initialize=chains_tot_proxy2_idx_rule)

model.chains_tot_proxy2_con = pyo.Constraint(
    model.chains_tot_proxy2_idx, rule=chains_tot_proxy2_rule)



# Chain debugging
is_chains_sweep_l_con_active = True
is_chains_sweep_u_con_active = True
is_chains_tot_con_active = True


def chains_sweep_l_rule(M, t, k, b, j, w, p, v):

    return sum(
        M.Shift[prd, day, wk, k, t]
        for (prd, day, wk) in M.linkspan[t, k, b, j, w, v + 1]
        if (prd, day, wk, k, t) in M.okShifts) >= \
           sum(M.TourTypeDayShift[epoch_to_tuple(M, u)[0], t, k, epoch_to_tuple(M, u)[1], epoch_to_tuple(M, u)[2]]
               for u in [vv for vv in range(p, p + v + 1)
                         if (epoch_to_tuple(M, vv)[0], epoch_to_tuple(M, vv)[1], epoch_to_tuple(M, vv)[2])
                         in M.okStartWindowRoots[t, k] and sum(M.allow_start[x, y, k, t] for (x, y, z)
                                                               in M.PotentialGlobalStartWindow[
                                                                   epoch_to_tuple(M, vv)[0], epoch_to_tuple(M, vv)[1],
                                                                   epoch_to_tuple(M, vv)[2]]) > 0])


def chains_sweep_l_idx_rule(M):
    index_list = []
    if M.g_start_window_width > 0:
        for t in M.activeTT:
            for k in [len for len in M.LENGTHS if len in M.tt_length_x[t]]:
                for b in M.PERIODS:
                    for j in M.DAYS:
                        for w in M.WEEKS:
                            for p in range(M.epoch[b, j, w], M.epoch[b, j, w] + 1):
                                if (b, j, w) in M.bchain[t, k]:
                                    for m in range(0,
                                                   M.n_links[t, k, b, j, w] - 1 - (p - M.epoch[b, j, w]) + 1):
                                        index_list.append((t, k, b, j, w, p, m))

    return index_list


model.chains_sweep_l_idx = pyo.Set(dimen=7, ordered=True,
                                   initialize=chains_sweep_l_idx_rule)

if is_chains_sweep_l_con_active:
    model.chains_sweep_l_con = pyo.Constraint(
        model.chains_sweep_l_idx, rule=chains_sweep_l_rule)


def chains_sweep_u_rule(M, t, k, b, j, w, p, v):
    return sum(M.Shift[epoch_to_tuple(M, i)[0], epoch_to_tuple(M, i)[1], epoch_to_tuple(M, i)[2], k, t]
               for i in range(p, p + v + 1) if
               (epoch_to_tuple(M, i)[0], epoch_to_tuple(M, i)[1], epoch_to_tuple(M, i)[2], k, t) in M.okShifts) <= \
           sum(M.TourTypeDayShift[epoch_to_tuple(M, u)[0], t, k, epoch_to_tuple(M, u)[1], epoch_to_tuple(M, u)[2]]
               for u in [vv for vv in range(p, p + v + 1)
                         if (epoch_to_tuple(M, vv)[0], epoch_to_tuple(M, vv)[1], epoch_to_tuple(M, vv)[2])
                         in M.okStartWindowRoots[t, k] and sum(M.allow_start[x, y, k, t] for (x, y, z)
                                                               in M.PotentialGlobalStartWindow[
                                                                   epoch_to_tuple(M, vv)[0], epoch_to_tuple(M, vv)[1],
                                                                   epoch_to_tuple(M, vv)[2]]) > 0])


def chains_sweep_u_idx_rule(M):
    index_list = []
    if M.g_start_window_width > 0:
        for t in M.activeTT:
            for k in [len for len in M.LENGTHS if len in M.tt_length_x[t]]:
                for b in M.PERIODS:
                    for j in M.DAYS:
                        for w in M.WEEKS:
                            for p in range(M.epoch[b, j, w], M.epoch[b, j, w] + 1):
                                if (b, j, w) in M.bchain[t, k]:
                                    for m in range(0, M.n_links[t, k, b, j, w] - 1):
                                        index_list.append((t, k, b, j, w, p, m))

    return index_list


model.chains_sweep_u_idx = pyo.Set(dimen=7, ordered=True,
                                   initialize=chains_sweep_u_idx_rule)

if is_chains_sweep_u_con_active:
    model.chains_sweep_u_con = pyo.Constraint(
        model.chains_sweep_l_idx, rule=chains_sweep_u_rule)


def chains_tot_rule(M, t, k, b, j, w):

    return sum(M.Shift[l, m, n, k, t] for (l, m, n) in
               [(p, q, r) for (p, q, r) in M.linkspan[t, k, b, j, w, M.n_links[t, k, b, j, w]] if
                (p, q, r, k, t) in M.okShifts]) - \
           sum(M.TourTypeDayShift[x, t, k, y, z]
               for (x, y, z) in [(u, v, xx) for (u, v, xx) in M.chain[t, k, b, j, w]
                                 if sum(
                   M.allow_start[p, q, k, t] for (p, q, r) in M.PotentialGlobalStartWindow[u, v, xx]) > 0]) == 0


def chains_tot_idx_rule(M):
    return [(t, k, b, j, w) for t in M.activeTT
            for k in [length for length in M.LENGTHS if length in M.tt_length_x[t]]
            for b in M.PERIODS
            for j in M.DAYS
            for w in M.WEEKS
            if M.g_start_window_width.value > 0 and (b, j, w) in M.bchain[t, k]]


model.chains_tot_idx = pyo.Set(dimen=5, ordered=True,
                               initialize=chains_tot_idx_rule)

if is_chains_tot_con_active:
    model.chains_tot_con = pyo.Constraint(
        model.chains_tot_idx, rule=chains_tot_rule)





def main():
    pass


if __name__ == '__main__':
    main()
