"""
Phase 1 for implicit multi-week tour scheduling model
"""

# Author: misken
# License: TBD

import pyomo.environ as pyo

# TODO Would be nice if Phase 1 and Phase 2 could share base pyo.Parameters
# model_phase1 = import_file('mwts_basepyo.Params.py').model

# Create Phase 1 abstract model
model_phase1 = pyo.AbstractModel()
model_phase1.name = "mwts_phase1"

# Constants
infinity = float('inf')

# General temporal parameters
model_phase1.n_prds_per_day = \
    pyo.Param(within=pyo.PositiveIntegers)  # n_P

model_phase1.n_days_per_week = \
    pyo.Param(within=pyo.PositiveIntegers)  # 7

model_phase1.n_weeks = \
    pyo.Param(within=pyo.PositiveIntegers)  # n_W


def n_prds_per_week_init(M):
    """
    Initialize convenience Parameter n_prds_per_week
    """
    return M.n_days_per_week() * M.n_prds_per_day()


model_phase1.n_prds_per_week = pyo.Param(within=pyo.PositiveIntegers,
                                         initialize=n_prds_per_week_init)


def n_prds_per_cycle_init(M):
    """
    Initialize convenience parameter n_prds_per_cycle where cycle may include
    one or more weeks.
    """
    return M.n_weeks() * M.n_days_per_week() * M.n_prds_per_day()


model_phase1.n_prds_per_cycle = \
    pyo.Param(within=pyo.PositiveIntegers, initialize=n_prds_per_cycle_init)


def g_period_init(M, i, j, w):
    """
    Initialize epoch index from daily period, day of week and week of cycle.

    :param M: Model
    :param i: period of day
    :param j: day of week
    :param w: week of cycle
    :return: period of cycle in 1..n_prds_per_cycle
    """
    return (w - 1) * M.n_days_per_week() * M.n_prds_per_day() + (j - 1) * M.n_prds_per_day() + i


# def g_tuple_to_period(i, j, w, n_prds_per_day, n_days_per_week):
#     return ((w - 1) * n_days_per_week * n_prds_per_day +
#             (j - 1) * n_prds_per_day + i)


# def period_increment(M, i, j, w, incr):
#     p = M.g_period[i, j, w]
#     if p + incr <= M.n_prds_per_cycle:
#         return p + incr
#     else:
#         return p + incr - M.n_prds_per_cycle


def g_period_increment(M, p, incr):
    """
    Increment period of cycle, circling as necessary

    :param M: Model
    :param p: period
    :param incr: Number of periods to increment by
    :return: incremented period
    """
    if p + incr <= M.n_prds_per_cycle:
        return p + incr
    else:
        return p + incr - M.n_prds_per_cycle


def g_period_difference(M, b_prd, e_prd):
    """
    Compute difference between two epochs accounting for end of horizon circling

    :param M: Model
    :param b_prd: global start period
    :param e_prd: global end period
    :return: end - start + 1 (+ n_prds_per_cycle if spans end of horizon)
    """

    if e_prd >= b_prd:
        return e_prd - b_prd + 1
    else:
        return M.n_prds_per_cycle + e_prd - b_prd + 1


def g_prd_to_tuple(M, p):
    """
    Convert epoch to (daily period, day, week) tuple

    :param M: Model
    :param p: epoch
    :return: (daily period, day, week) tuple
    """

    week = ((p - 1) // M.n_prds_per_week.value) + 1
    prds_remainder = p - (week - 1) * M.n_prds_per_week
    if prds_remainder == 0:
        day = 1
    else:
        day = ((prds_remainder - 1) // M.n_prds_per_day.value) + 1

    prds_remainder = prds_remainder - (day - 1) * M.n_prds_per_day
    if prds_remainder == 0:
        period = 1
    else:
        period = prds_remainder

    return period, day, week


# Temporal ranges used for various index sets
model_phase1.PERIODS = pyo.RangeSet(1, model_phase1.n_prds_per_day)
model_phase1.EPOCHS = pyo.RangeSet(1, model_phase1.n_prds_per_cycle)
model_phase1.WINDOWS = pyo.RangeSet(1, model_phase1.n_prds_per_day)
model_phase1.DAYS = pyo.RangeSet(1, model_phase1.n_days_per_week)
model_phase1.WEEKS = pyo.RangeSet(1, model_phase1.n_weeks)
model_phase1.WEEKENDS = pyo.RangeSet(1, 2)

model_phase1.g_period = pyo.Param(model_phase1.PERIODS, model_phase1.DAYS,
                                  model_phase1.WEEKS, initialize=g_period_init)

model_phase1.epoch_tuples = model_phase1.PERIODS * model_phase1.DAYS * model_phase1.WEEKS  # B


# ### Tour type related parameters

# -- Shift Lengths
# Number of different shift lengths
model_phase1.n_lengths = pyo.Param(within=pyo.PositiveIntegers)
# Range of length indexes
model_phase1.LENGTHS = pyo.RangeSet(1, model_phase1.n_lengths)
# Vector of shift lengths
model_phase1.lengths = pyo.Param(model_phase1.LENGTHS)

# -- Tour Types
# Number of different tour types
model_phase1.n_tts = pyo.Param(within=pyo.PositiveIntegers)
# Range of tour type indexes
model_phase1.TTYPES = pyo.RangeSet(1, model_phase1.n_tts)
# Set of allowable length indices by tour type
model_phase1.tt_length_x = pyo.Set(model_phase1.TTYPES, ordered=True,)
# Bounds on tour type pyo.Variables
model_phase1.tt_lb = pyo.Param(model_phase1.TTYPES)
model_phase1.tt_ub = pyo.Param(model_phase1.TTYPES,  default=infinity)


def activeTT_init(M):
    """
    Initialize list of tour types that can be used in this problem instance

    :param M: Model
    :return: list of tour type indexes
    """
    return [t for t in M.TTYPES if M.tt_ub[t] > 0]


model_phase1.activeTT = pyo.Set(dimen=1, ordered=True, initialize=activeTT_init)


# -- Weekend patterns
def max_wkend_patterns(M):
    """
    Compute maximum number of weekend worked patterns

    :param M:
    :return:
    """
    max_patterns = 2 ** (2 * M.n_weeks.value)
    return max_patterns


model_phase1.max_weekend_patterns = pyo.Param(initialize=max_wkend_patterns)

model_phase1.num_weekend_patterns = pyo.Param(model_phase1.WEEKENDS,
                                              model_phase1.TTYPES)


def A_idx_rule(M):
    """
    Construct index for weekend pattern variables

    A[i,j,w,t,e] = 1 if weekend pattern i calls for work on day j of week k
    for tour type t having weekend type e, and 0 otherwise

    :param M:
    :return: list of tuples of indexes
    """
    return [(i, j, w, t, e) for i in pyo.sequence(M.max_weekend_patterns)
            for j in M.DAYS
            for w in M.WEEKS
            for t in M.TTYPES
            for e in pyo.sequence(2) if i <= M.num_weekend_patterns[e, t]]


model_phase1.A_idx = pyo.Set(dimen=5, ordered=True, initialize=A_idx_rule)
model_phase1.A = pyo.Param(model_phase1.A_idx, within=pyo.Boolean, default=0)


# def A_wkend_days_idx_rule(M):
#     return [(i, w, t, e) for i in pyo.sequence(M.max_weekend_patterns)
#             for w in M.WEEKS
#             for t in M.TTYPES
#             for e in pyo.sequence(2) if i <= M.num_weekend_patterns[e, t]]
#
#
# model_phase1.A_wkend_days_idx = pyo.Set(dimen=4,ordered=True,
#                                         initialize=A_wkend_days_idx_rule)


# -- Multiweek days worked patterns

# TODO: Review max_mwdw_patterns
def max_mwdw_init(M):
    max_mwdw = 4 ** M.n_weeks.value
    return max_mwdw


model_phase1.max_mwdw_patterns = pyo.Param(initialize=max_mwdw_init)
model_phase1.num_mwdw_patterns = pyo.Param(model_phase1.TTYPES)


def A_mwdw_idx_rule(M):
    return [(t, p, w) for t in M.TTYPES
            for p in pyo.sequence(M.max_mwdw_patterns)
            for w in M.WEEKS
            if p <= M.num_mwdw_patterns[t]]


model_phase1.A_mwdw_idx = pyo.Set(dimen=3, ordered=True, initialize=A_mwdw_idx_rule)
model_phase1.A_mwdw = pyo.Param(model_phase1.A_mwdw_idx,
                                within=pyo.NonNegativeIntegers, default=0)


# Number of weekend days worked within a week
def A_num_wkend_days_init(M, i, w, t, e):
    """
    Initialize number of weekend days worked in the i'th pattern, for week k, ttype t,
    and weekend type e.

    :param M: Model
    :param i: weekend worked pattern
    :param w: week
    :param t: tour type
    :param e: weekend type (1=Sun and Sat, 2=Fri and Sat)
    :return:
    """

    if e == 1:
        return M.A[i, 1, w, t, e] + M.A[i, 7, w, t, e]

    else:
        return M.A[i, 6, w, t, e] + M.A[i, 7, w, t, e]
        
    
model_phase1.A_num_wkend_days = pyo.Param(model_phase1.A_wkend_days_idx,
                                          initialize=A_num_wkend_days_init)


# Number of weekend days worked over entire scheduling cycle
def A_tot_wkend_days_idx_rule(M):
    return [(i, t, e) for i in pyo.sequence(M.max_weekend_patterns)
            for t in M.TTYPES
            for e in pyo.sequence(2) if i <= M.num_weekend_patterns[e, t]]


model_phase1.A_tot_wkend_days_idx = pyo.Set(dimen=3, ordered=True,
                                            initialize=A_tot_wkend_days_idx_rule)


def A_tot_wkend_days_init(M, i, t, e):
    """
    Initialize number of weekend days worked in the i'th pattern, ttype t,
    and weekend type e.

    :param M: Model
    :param i: weekend worked pattern
    :param t: tour type
    :param e: weekend type (1=Sun and Sat, 2=Fri and Sat)
    :return:
    """

    if e == 1:
        return sum(M.A[i, 1, w, t, e] + M.A[i, 7, w, t, e] for w in M.WEEKS)
    else:
        return sum(M.A[i, 6, w, t, e] + M.A[i, 7, w, t, e] for w in M.WEEKS)


model_phase1.A_tot_wkend_days = pyo.Param(model_phase1.A_tot_wkend_days_idx,
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
    """

    """
    if e == 1:
        if M.A[i, 1, w, t, e] == 1 and M.A[i, 7, w, t, e] == 1:
            return 1
        else:
            return 0

    else:
        if M.A[i, 6, w, t, e] == 1 and M.A[i, 7, w, t, e] == 1:
            return 1
        else:
            return 0
        
    
model_phase1.A_is_two_wkend_days = pyo.Param(model_phase1.A_wkend_days_idx,
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
        if (M.A[i, 1, w, t, e] + M.A[i, 7, w, t, e]) == 1:
            return 1
        else:
            return 0

    else:
        if (M.A[i, 6, w, t, e] + M.A[i, 7, w, t, e]) == 1:
            return 1
        else:
            return 0
        
    
model_phase1.A_is_one_wkend_days = pyo.Param(model_phase1.A_wkend_days_idx,
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
        if M.A[i, 1, w, t, e] == 1 and M.A[i, 7, w, t, e] == 0:
            return 1
        else:
            return 0

    else:
        # TODO - FS weekend not implemented
        if M.A[i, 1, w, t, e] == 1 and M.A[i, 7, w, t, e] == 0:
            return 1
        else:
            return 0


model_phase1.A_is_Sunday = pyo.Param(model_phase1.A_wkend_days_idx,
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
        if M.A[i, 1, w, t, e] == 0 and M.A[i, 7, w, t, e] == 1:
            return 1
        else:
            return 0

    else:
        # TODO - FS weekend not implemented
        if M.A[i, 1, w, t, e] == 0 and M.A[i, 7, w, t, e] == 1:
            return 1
        else:
            return 0


model_phase1.A_is_Saturday = pyo.Param(model_phase1.A_wkend_days_idx,
                                       initialize=A_is_Saturday_init)


# Weekends consisting of a Fri and Sat imply updated lower bounds on
# some of the daily tour type pyo.Variables.
# These were the key to modeling weekends worked patterns.

# def FriSat_idx_rule(M):
#     # return [(t,e,w) for t in M.activeTT for e in range(1,M.num_weekend_patterns[2,t] + 1) for w in model_phase1.WEEKS]
#     index_list = []
#     for t in M.activeTT:
#         numpats = M.num_weekend_patterns[2, t]
#         for p in pyo.sequence(numpats):
#             for w in M.WEEKS:
#                 index_list.append((t, p, w))
#
#     return index_list
#
#
# model_phase1.FriSat_idx = pyo.Set(dimen=3, ordered=True, initialize=FriSat_idx_rule)
#
# model_phase1.FriSat_min_dys_weeks = pyo.Param(model_phase1.FriSat_idx, default=0)
# model_phase1.FriSat_min_cumul_dys_weeks = pyo.Param(model_phase1.FriSat_idx, default=0)

# Bounds on days and shifts worked over the week
# 1a. Min and max number of days worked by week by tour type
model_phase1.tt_min_dys_weeks = pyo.Param(model_phase1.TTYPES,
                                          model_phase1.WEEKS, default=0.0)

model_phase1.tt_max_dys_weeks = pyo.Param(model_phase1.TTYPES,
                                          model_phase1.WEEKS, default=1e+6)

# 1b. Min and max number of days worked by cumulative weeks by tour type
model_phase1.tt_min_cumul_dys_weeks = pyo.Param(model_phase1.TTYPES,
                                                model_phase1.WEEKS, default=0.0)

model_phase1.tt_max_cumul_dys_weeks = pyo.Param(model_phase1.TTYPES,
                                                model_phase1.WEEKS, default=1e+6)

# 2a. Min and max number of days worked by week by shiftlen by tour type
model_phase1.tt_shiftlen_min_dys_weeks = pyo.Param(model_phase1.TTYPES,
                                                   model_phase1.LENGTHS,
                                                   model_phase1.WEEKS, default=0)

model_phase1.tt_shiftlen_max_dys_weeks = pyo.Param(model_phase1.TTYPES,
                                                   model_phase1.LENGTHS,
                                                   model_phase1.WEEKS, default=1e+6)

# 2b. Min and max number of days worked by cumulative weeks by shiftlen by tour type
model_phase1.tt_shiftlen_min_cumul_dys_weeks = pyo.Param(model_phase1.TTYPES,
                                                         model_phase1.LENGTHS,
                                                         model_phase1.WEEKS, default=0)

model_phase1.tt_shiftlen_max_cumul_dys_weeks = pyo.Param(model_phase1.TTYPES,
                                                         model_phase1.LENGTHS,
                                                         model_phase1.WEEKS, default=1e+6)

# 3a. Min and max number of periods worked by week by tour type
model_phase1.tt_min_prds_weeks = pyo.Param(model_phase1.TTYPES,
                                           model_phase1.WEEKS, default=0)

model_phase1.tt_max_prds_weeks = pyo.Param(model_phase1.TTYPES,
                                           model_phase1.WEEKS, default=1e+6)

# 3b. Min and max number of periods worked by cumulative weeks by tour type
model_phase1.tt_min_cumul_prds_weeks = pyo.Param(model_phase1.TTYPES,
                                                 model_phase1.WEEKS, default=0)

model_phase1.tt_max_cumul_prds_weeks = pyo.Param(model_phase1.TTYPES,
                                                 model_phase1.WEEKS, default=1e+6)

# 4a. Min and max number of periods worked by week by tour type by shift length
model_phase1.tt_shiftlen_min_prds_weeks = pyo.Param(model_phase1.TTYPES,
                                                    model_phase1.LENGTHS,
                                                    model_phase1.WEEKS, default=0)

model_phase1.tt_shiftlen_max_prds_weeks = pyo.Param(model_phase1.TTYPES,
                                                    model_phase1.LENGTHS,
                                                    model_phase1.WEEKS, default=1e+6)

# 4b. Min and max number of periods worked by cumulative weeks by tour type by shift length
model_phase1.tt_shiftlen_min_cumul_prds_weeks = pyo.Param(model_phase1.TTYPES,
                                                          model_phase1.LENGTHS,
                                                          model_phase1.WEEKS, default=0)

model_phase1.tt_shiftlen_max_cumul_prds_weeks = pyo.Param(model_phase1.TTYPES,
                                                          model_phase1.LENGTHS,
                                                          model_phase1.WEEKS, default=1e+6)


# To the above, we'll add parameters and sets to allow direct modeling of side constraints
# of the form sum{subset of tour types} =, >=, <= some bound

# Indicator for part-time tour types: 1 for part-time, 0 for full-time
model_phase1.tt_parttime = pyo.Param(model_phase1.TTYPES)

# Allowable shift start times  - note that these are tour type and
# shift length specific
model_phase1.allow_start = pyo.Param(model_phase1.PERIODS,
                                     model_phase1.DAYS,
                                     model_phase1.LENGTHS,
                                     model_phase1.TTYPES,
                                     default=0.0)


def okShifts_rule(M):
    return [(i, j, w, k, t) for i in M.PERIODS
            for j in M.DAYS
            for w in M.WEEKS
            for t in M.activeTT
            for k in M.tt_length_x[t]
            if M.allow_start[i, j, k, t] > 0]


model_phase1.okShifts = pyo.Set(dimen=5,ordered=True,initialize=okShifts_rule)


def okShiftTypes_rule(M):
    return [(i, j, k, t) for i in M.PERIODS
            for j in M.DAYS
            for t in M.activeTT
            for k in M.tt_length_x[t]
            if M.allow_start[i, j, k, t] > 0]


model_phase1.okShiftTypes = pyo.Set(dimen=4, ordered=True, initialize=okShiftTypes_rule)


# Limits on part time labor and limits on total labor

model_phase1.max_parttime_frac = pyo.Param()    # Maximum fraction of labor hours covered by part-time employees
model_phase1.labor_budget = pyo.Param()         # Maximum labor expenditure

# -----------------------------------------------------------------------
# Cost related parameters
# -----------------------------------------------------------------------

model_phase1.tt_cost_multiplier = pyo.Param(model_phase1.TTYPES)           # Tour type differential

model_phase1.cu1 = pyo.Param()
model_phase1.cu2 = pyo.Param()        
model_phase1.usb = pyo.Param()

# -----------------------------------------------------------------------
# Weekend related Parameters
# -----------------------------------------------------------------------

# Need to generalize for any number of periods per day
model_phase1.midnight_thresh = pyo.Param(model_phase1.TTYPES, default=1e+6)


def weekend_init(M, i, t):
    result = []
    lens = [M.lengths[k] for k in M.tt_length_x[t]]
    maxlen = max(lens)
    if i + maxlen - 1 >= M.midnight_thresh[t]:
        result.append(6)
    else:
        result.append(1)
    result.append(7)
    return result


model_phase1.weekend = pyo.Set(model_phase1.WINDOWS, model_phase1.TTYPES, ordered=True, initialize=weekend_init)


def weekend_type_init(M, i, t):
    result = 1
    lens = [M.lengths[k] for k in M.tt_length_x[t]]
    maxlen = max(lens)
    if i + maxlen - 1 >= M.midnight_thresh[t]:
        result = 2

    return result
        
          
model_phase1.weekend_type = pyo.Param(model_phase1.WINDOWS, model_phase1.TTYPES, initialize=weekend_type_init)



def zero_wkend_day_init(M,w,t,e):
    set_list = []
    pattern_list = [(p,W,T,E) for (p,W,T,E) in M.A_wkend_days_idx if W==w and T==t and E==e]
    
    for (p,W,T,E) in pattern_list:
        if M.A_is_one_wkend_days[p,W,T,E] == 0 and M.A_is_two_wkend_days[p,W,T,E] == 0:
            set_list.append(p)
        
    return set_list


model_phase1.zero_wkend_day = pyo.Set(model_phase1.WEEKS, model_phase1.TTYPES, model_phase1.WEEKENDS,
                                  initialize=zero_wkend_day_init)


def one_wkend_day_init(M, w, t, e):
    set_list = []
    for (p, W, T, E) in M.A_wkend_days_idx:
        if W == w and T == t and E == e and M.A_is_one_wkend_days[p, w, t, e] == 1:
            set_list.append(p)
        
    return set_list


model_phase1.one_wkend_day = pyo.Set(model_phase1.WEEKS, model_phase1.TTYPES, model_phase1.WEEKENDS,
                                 initialize=one_wkend_day_init)


def Sun_wkend_day_init(M, w, t, e):
    set_list = []
    for (p, W, T, E) in M.A_wkend_days_idx:
        if W == w and T == t and E == e and M.A_is_Sunday[p, w, t, e]:
            set_list.append(p)
        
    return set_list


model_phase1.Sun_wkend_day = pyo.Set(model_phase1.WEEKS, model_phase1.TTYPES, model_phase1.WEEKENDS,
                                 initialize=Sun_wkend_day_init)


def Sat_wkend_day_init(M, w, t, e):
    set_list = []
    for (p, W, T, E) in M.A_wkend_days_idx:
        if W == w and T == t and E == e and M.A_is_Saturday[p, w, t, e]:
            set_list.append(p)
        
    return set_list


model_phase1.Sat_wkend_day = pyo.Set(model_phase1.WEEKS, model_phase1.TTYPES, model_phase1.WEEKENDS,
                                 initialize=Sat_wkend_day_init)


def oneorzero_wkend_day_init(M, w, t, e):
    set_list = []
    for (p, W, T, E) in M.A_wkend_days_idx:
        if W == w and T == t and E == e and M.A_is_two_wkend_days[p, w, t, e] == 0:
            set_list.append(p)
        
    return set_list


model_phase1.oneorzero_wkend_day = pyo.Set(model_phase1.WEEKS, model_phase1.TTYPES, model_phase1.WEEKENDS,
                                       initialize=oneorzero_wkend_day_init)


def two_wkend_days_init(M, w, t, e):
    set_list = []
    for (p, W, T, E) in M.A_wkend_days_idx:
        if W == w and T == t and E == e and M.A_is_two_wkend_days[p, w, t, e] == 1:
            set_list.append(p)

    return set_list


model_phase1.two_wkend_days = pyo.Set(model_phase1.WEEKS, model_phase1.TTYPES, model_phase1.WEEKENDS,
                                  initialize=two_wkend_days_init)



  
# -----------------------------------------------------------------------
# Coverage related pyo.ParamETERS
# -----------------------------------------------------------------------

# Target and minimum staffing levels - this is week specific. We can always allow user to input
# a single week and then repeat it for the other weeks.

model_phase1.dmd_staff = pyo.Param(model_phase1.PERIODS,model_phase1.DAYS,model_phase1.WEEKS)
model_phase1.min_staff = pyo.Param(model_phase1.PERIODS,model_phase1.DAYS,model_phase1.WEEKS)

# -----------------------------------------------------------------------
# START WINDOWS - should these be tour type specific since allow start is tour type specific?
#
# To start with,I made them week specific but not sure if they should be. Seems they should be 
# period, day, tour type.
#
# Maybe best strategy is to do the pyo.Variables and constraints first and then work backwards to
# define windows appropriately.
# -----------------------------------------------------------------------

model_phase1.g_start_window_width = pyo.Param()                            # Width of start-time windows



##/**** Beginning of each start window (in total periods from Sunday @ midnight)****/
##pyo.Param b_window_wepoch{i in PERIODS,j in DAYS} := n_prds_per_day*(j-1)+i;
##
##/**** End of each start window (in total periods from Sunday @ midnight) ****/
##pyo.Param e_window_wepoch{i in PERIODS,j in DAYS} :=
## ( if n_prds_per_day*(j-1)+i+width <= n_prds_per_day*n_days then
##    n_prds_per_day*(j-1)+i+width
##  else
##    (n_prds_per_day*(j-1)+i+width )-n_prds_per_day*n_days);


# if model_phase1.g_start_window_width > 0:
# This version spans multiple weeks
def b_window_epoch_init(M, i, j, w):
    return M.n_days_per_week.value*M.n_prds_per_day.value*(w-1) + M.n_prds_per_day.value*(j-1)+i

model_phase1.b_window_epoch = pyo.Param(model_phase1.epoch_tuples, initialize=b_window_epoch_init)


def e_window_epoch_init(M,i,j,w):
    epoch = M.n_days_per_week.value * M.n_prds_per_day.value*(w-1) + \
            M.n_prds_per_day.value * (j-1) + i + M.g_start_window_width.value

    if epoch <= M.n_prds_per_cycle.value:
        e_window = epoch
    else:
        e_window = epoch-M.n_prds_per_cycle.value

    return e_window

model_phase1.e_window_epoch = pyo.Param(model_phase1.epoch_tuples, initialize=e_window_epoch_init)



###### PotentialGlobalStartWindow ######

# Index: PotentialGlobalStartWindow is defined for all (period,day,week) epoch_tuples.
# Defn:  PotentialGlobalStartWindow[i,j,w] contains all the bin triplets within g_start_window_width periods of
#    bin (i,j,w). The "Potential" indicates that many of these will be eliminated because [i,j,w] may
#    not be an allowable shift start time. Notice that these windows are NOT yet shift length and tour
#    type specific - hence the "Global" in the name.

# Single week model equivalent: WindowWepochs

# EXAMPLE: Let g_start_window_width=2. Then PotentialGlobalStartWindow[5,2,1] = [(5,2,1),(6,2,1),(7,2,1)]

def PotentialGlobalStartWindow_init(M,i,j,w):
    window_list =[]
    for (l,m,n) in M.bins:

        test_epoch = M.g_period[l,m,n]

        test1 = (test_epoch >= M.b_window_epoch[i,j,w])

        if M.b_window_epoch[i,j,w] <= M.e_window_epoch[i,j,w]:
            test2rhs = M.e_window_epoch[i,j,w]
        else:
            test2rhs = M.n_prds_per_cycle

        test2 = (test_epoch <= test2rhs)

        if M.b_window_epoch[i,j,w] <= M.e_window_epoch[i,j,w]:
            test3rhs = M.b_window_epoch[i,j,w]
        else:
            test3rhs = 1

        test3 = (test_epoch >= test3rhs)

        test4 = (test_epoch <= M.e_window_epoch[i,j,w])

        if (test1 and test2) or (test3 and test4):
            window_list.append((l,m,n))

    return window_list

model_phase1.PotentialGlobalStartWindow = pyo.Set(model_phase1.PERIODS,model_phase1.DAYS,model_phase1.WEEKS,dimen=3,ordered=True,initialize=PotentialGlobalStartWindow_init)




#
#
####/**** The set okWindowWepochs{i in PERIODS,j in DAYS,k in LENGTHS,t in TTYPES } creates shift
####length and tour type specific sets of windows that start in (i,j) pairs which have
####(i,j,k) as an allowable shift start time. ****/
####
####set okWindowWepochs{i in PERIODS,j in DAYS,k in LENGTHS,t in TTYPES}
#### :={(l,m) in WindowWepochs[i,j]: allow_start[i,j,k,t]>0 and allow_start[l,m,k,t]>0};
#

###### PotentialStartWindow ######

# Index: PotentialStartWindow is defined for all (period,day,week,shift length,tour type).
# Defn:  PotentialStartWindow[i,j,k,t] contains all the bin triplets within g_start_window_width periods of
#    bin (i,j,w) and having (i,j,k,t) be an allowable start time and all of its elements being allowable start times.
#    The "Potential" indicates that some of these will be eliminated because they may
#    be subsets of other, nearby, start windows.

# Single week model equivalent: okWindowWepochs

# EXAMPLE: Let g_start_window_width=2 and only odd periods be allowable start time periods.
#    Then PotentialStartWindow[5,2,1,1,1] = [(5,2,1),(7,2,1)]

def PotentialStartWindow_idx_rule(M):
    return [(i,j,w,k,t) for i in M.PERIODS
                        for j in M.DAYS
                        for w in M.WEEKS
                        for k in M.LENGTHS
                        for t in M.activeTT]


model_phase1.PotentialStartWindow_idx = pyo.Set(dimen=5,ordered=True,initialize=PotentialStartWindow_idx_rule)

def PotentialStartWindow_init(M,i,j,w,k,t):
    return [(l,m,n) for (l,m,n) in M.PotentialGlobalStartWindow[i,j,w]
                  if M.allow_start[i,j,k,t] > 0 and M.allow_start[l,m,k,t] > 0]

model_phase1.PotentialStartWindow = pyo.Set(model_phase1.PotentialStartWindow_idx,ordered=True,initialize=PotentialStartWindow_init,dimen=3)
#
#
###### okStartWindowRoots ######

# Index: okStartWindowRoots is defined for all (tour type,shift length) pairs such that the shift length
#        is allowed for the tour type.

# Defn:  okStartWindowRoots[t,k] contains all the bin triplets in PotentialStartWindow

# Single week model equivalent: ok_window_beginnings

# EXAMPLE: Let g_start_window_width=2 and only odd periods be allowable start time periods.
#    Then PotentialStartWindow[5,2,1,1,1] = [(5,2,1),(7,2,1)]


def okStartWindowRoots_idx_rule(M):
    return [(t,k) for k in M.LENGTHS
                         for t in M.activeTT
                         if k in M.tt_length_x[t]]


model_phase1.okStartWindowRoots_idx = pyo.Set(dimen=2,ordered=True,initialize=okStartWindowRoots_idx_rule)
#
###/**** ok_window_beginnings is the set of start windows in which there is
###at least one period in which a shift of length k can start
###and which are not subsets of some other window.
###
###*/
##
###set okWindowBeginnings{t in okTTYPES, k in tt_length_x[t]} :=
###  setof{(p,q) in {PERIODS,DAYS}: (p,q,k,t) in ok_shifts and
###   forall{(i,j) in {PERIODS,DAYS} diff {(p,q)}:allow_start[i,j,k,t]>0 }
###    (not
###     ({(l,m) in okWindowWepochs[p,q,k,t]}
###      within {(n,o) in okWindowWepochs[i,j,k,t]}))  } (p,q);
##
def okStartWindowRoots_init(M,t,k):
    window_list =[]
    for (p,q,w) in M.bins:
        test1 = False


        # Create a list of all possible window beginnings and then eliminate those that
        # have window wepochs that are subsets of other window wepochs
        test1 = ((p,q,k,t) in M.okShiftTypes)
        if test1:
            window = []
            for (i,j,m) in M.PotentialStartWindow[p,q,w,k,t]:
                window.append((i,j,m))
            window_list.append(window)
#
#    # Get rid of subsets
    window_list_copy = window_list[:]
#
    for s1 in window_list:
        for s2 in window_list:
            if set(s1) < set(s2) and s1!=s2:
                window_list_copy.remove(s1)
                break

#   Now find the earliest bin in each element (list) of the list
    windowroot_list = []
    for w in window_list_copy:
        onSaturday = False
        onSunday = False
        early = M.n_prds_per_cycle
        satearly = M.n_prds_per_cycle
        late = 0
        for btup in w:
            if btup[1] == 7:
                onSaturday = True
            if btup[1] == 1:
                onSunday = True

            binofweek = M.g_period[btup[0],btup[1],btup[2]]
            if binofweek < early:
                early = binofweek
            if binofweek > late:
                late = binofweek
            if btup[1] == 7 and binofweek < satearly:
                satearly = binofweek

        if (onSaturday and onSunday):
            early = satearly

        earlybin = g_prd_to_tuple(M,early)
        windowroot_list.append(earlybin)
        #print earlybin

    return windowroot_list


model_phase1.okStartWindowRoots = pyo.Set(model_phase1.okStartWindowRoots_idx,dimen=3,ordered=True,initialize=okStartWindowRoots_init)

##
##
### CHAINS - NON-OVERLAPPING SEQUENCES OF (i,j) period,day PAIRS THAT
###          CAN BE ISOLATED FOR COORDINATING ix AND DWT pyo.VarIABLES.
###
### Let's wait on generalizing these for multiple weeks. Get model working
### for start window width of 0.
### -----------------------------------------------------------------------
##
##
###set bchain {t in okTTYPES, k in tt_length_x[t]} := setof{(w,j) in (okWindowBeginnings[t,k]):
###    forall{(p,q) in (okWindowBeginnings[t,k] diff {(w,j)})} (w,j) not in (okWindowWepochs[p,q,k,t])} (w,j) ;
##
##
def bchain_init(M,t,k):
    window_list =[]
    if M.g_start_window_width>0:
        for (i,j,w) in M.okStartWindowRoots[t,k]:
            for (p,q,r) in (M.okStartWindowRoots[t,k] - Set(initialize=[(i,j,w)])):
                if (i,j,w) not in M.PotentialStartWindow[p,q,r,k,t]:
                    window_list.append((i,j,w))

    return window_list
#
#


model_phase1.bchain = pyo.Set(model_phase1.okStartWindowRoots_idx, dimen=3, ordered=True, initialize=bchain_init)

#
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
    return [(t,k,i,j,w) for t in M.activeTT
                      for k in M.LENGTHS
                      for i in M.PERIODS
                      for j in M.DAYS
                      for w in M.WEEKS
                      if (t,k) in M.okStartWindowRoots_idx and (i,j,w) in M.bchain[t,k]]


model_phase1.chain_idx = pyo.Set(dimen=5, ordered=True, initialize=chain_idx_rule)
#
def echain_init(M,t,k,i,j,w):
    window_list =[]
    # Compute epoch of (i,j,w)
    g_prd = M.g_period[i, j, w]
    prd = g_prd_to_tuple(M, g_prd)
    #done = False

    steps = 1
    while steps <= M.g_start_window_width:
        g_prd_next = g_period_increment(M,g_prd,steps)
        prd_next = g_prd_to_tuple(M, g_prd_next)
        if prd_next in M.PotentialStartWindow[prd[0],prd[1],prd[2],k,t]:
            g_prd = g_prd_next
            prd = g_prd_to_tuple(M, g_prd)
            steps = 1
        else:
            steps = steps + 1


    # Once we leave the while loop, g_prd should correspond to echain. Need
    # to convert it back to a tuple from a epoch number.


    window_list.append(prd)
    return window_list


model_phase1.echain = pyo.Set(model_phase1.chain_idx,ordered=True,dimen=3,initialize=echain_init)


def n_links_init(M,t,k,i,j,w):

    # Compute epoch of (i,j,w)
    b_g_prd = M.g_period[i, j, w]
    e_prd = M.echain[t, k, i, j, w][1]
    e_g_prd = M.g_period[e_prd[0], e_prd[1], e_prd[2]]

    return g_period_difference(M,b_g_prd,e_g_prd)


model_phase1.n_links = pyo.Param(model_phase1.chain_idx,initialize=n_links_init)


def chain_init(M,t,k,i,j,w):
    window_list =[(i,j,w)]
    # Compute epoch of (i,j,w)
    g_prd = M.g_period[i, j, w]
    #done = False

    steps = 1
    while steps <= M.n_links[t,k,i,j,w]:
        g_prd_next = g_period_increment(M,g_prd,steps)
        prd_next = g_prd_to_tuple(M, g_prd_next)
        if prd_next in M.okStartWindowRoots[t,k]:
            window_list.append(prd_next)
        steps = steps + 1


    return window_list


model_phase1.chain = pyo.Set(model_phase1.chain_idx,ordered=True,dimen=3,initialize=chain_init)


def link_idx_rule(M):

# TODO Check the m index to see what the upper index limit should be
    return [(t,k,i,j,w,m) for t in M.activeTT
                          for k in M.LENGTHS
                          for i in M.PERIODS
                          for j in M.DAYS
                          for w in M.WEEKS
                          for m in M.CYCLEPERIODS
                          if (t,k) in M.okStartWindowRoots_idx and (i,j,w) in M.bchain[t,k]
                              and m <= M.n_links[t,k,i,j,w]]



model_phase1.link_idx = pyo.Set(dimen=6, ordered=True, initialize=link_idx_rule)


def link_init(M,t,k,i,j,w,m):
    """
    Returns the m'th link (a (period,day,week) tuple) of the chain starting in period (i,j,w)
    """
    window_list =[]
    g_prd = M.g_period[i, j, w]
    g_prd_next = g_period_increment(M,g_prd,m-1)
    prd_next = g_prd_to_tuple(M, g_prd_next)
    window_list.append(prd_next)

    return window_list


model_phase1.link = pyo.Set(model_phase1.link_idx,ordered=True,dimen=3,initialize=link_init)


def linkspan_init(M,t,k,i,j,w,m):
    """
    Returns the start windows spanned by the m'th link (a (period,day,week) tuple) of the chain
    starting in period (i,j,w)
    """
    window_list =[]
    for (p,d,q) in M.link[t,k,i,j,w,m]:
        for (a,b,c) in M.PotentialGlobalStartWindow[p,d,q]:
            window_list.append((a,b,c))
        if m > 1:
            for (a,b,c) in M.linkspan[t,k,i,j,w,m-1]:
                window_list.append((a,b,c))

    return window_list


model_phase1.linkspan = pyo.Set(model_phase1.link_idx, dimen=3, ordered=True, initialize=linkspan_init)




# #### Shift Variables

model_phase1.Shift = pyo.Var(model_phase1.okShifts, within=pyo.NonNegativeIntegers)

# Shift[i,j,w,k,t] = Number of shifts of length k starting in period i
# of day j in week w for a tour of type t

# #### Tour type Variables
# TourType[i,j] Number of employees working tour type j starting in window i  

def okTourType_rule(M):
    """
    List of (window,tour type) tuples that are allowable. To be allowable,
    for every week there must be at least the minumum required number of
    days worked having an allowable shift (of any length allowed for that 
    tour type).
    """
    index_list =[]
    for (i,t) in M.WINDOWS * M.activeTT:
        n_ok_weeks = 0
        for w in M.WEEKS:
            n_ok_days = 0
            for j in M.DAYS:
                for k in M.tt_length_x[t]:
                    if (i,j,w) in M.okStartWindowRoots[t,k]:
                        n_ok_days += 1
                        # Break out of the inner for since this day is ok
                        # for at least one shift length
                        break
            if n_ok_days >= M.tt_min_dys_weeks[t,w]:
                n_ok_weeks += 1
        if n_ok_weeks == M.n_weeks:
            index_list.append((i,t))
    
    return index_list


model_phase1.okTourType = pyo.Set(dimen=2, initialize=okTourType_rule)


def TourType_idx_rule(M):
    return [(i, t) for i in M.WINDOWS
            for t in M.activeTT
            if (i, t) in M.okTourType]
                        

model_phase1.TourType_idx = pyo.Set(dimen=2, initialize=TourType_idx_rule)

model_phase1.TourType = pyo.Var(model_phase1.TourType_idx, within=pyo.NonNegativeIntegers)


# /*set okDailyTourType := setof{i in WINDOWS,t in okTTYPES, d in DAYS :
#    (i,d) in ok_window_beginnings[t]} (i,t,d);*/
#
#
# set okDailyTourType := setof{i in WINDOWS,t in okTTYPES, d in DAYS} (i,t,d);
    
# #### Daily tour type pyo.Variables


model_phase1.okDailyTourType = model_phase1.okTourType * model_phase1.DAYS


#    /* DailyTourType[i,t,d] Number of employees working tour type t
#       starting in window i and working day d in week w*/

def DailyTourType_idx_rule(M):
    return [(i, t, j, w) for i in M.WINDOWS
            for t in M.activeTT
            for j in M.DAYS
            for w in M.WEEKS
            if (i, t, j) in M.okDailyTourType]
                        

model_phase1.DailyTourType_idx = pyo.Set(dimen=4, initialize=DailyTourType_idx_rule)

model_phase1.DailyTourType = pyo.Var(model_phase1.DailyTourType_idx, within=pyo.NonNegativeIntegers)



##### Daily shift worked pyo.Variables

#pyo.Var DailyShiftWorked{i in 1..n_windows,t in okTTYPES,k in tt_length_x[t],d in DAYS,w in WEEKS : (i,t,d) in okDailyTourType}
#    >= 0, integer;

# TODO - modify this index when get width>0 working. These are windows, not period start times.
# Just using okShifts for now for w=0 case.  
def DailyShiftWorked_idx_rule(M):
    return [(i, t, k, j, w) for i in M.WINDOWS
            for t in M.activeTT
            for k in M.tt_length_x[t]
            for j in M.DAYS
            for w in M.WEEKS
            if (i, t, j) in M.okDailyTourType]
                        

model_phase1.DailyShiftWorked_idx = pyo.Set(dimen=5, initialize=DailyShiftWorked_idx_rule)

model_phase1.DailyShiftWorked = pyo.Var(model_phase1.DailyShiftWorked_idx, within=pyo.NonNegativeIntegers)

# #### Weekend Days off pyo.Variables

# pyo.Var WeekendDaysWorked{p in 1..max_weekend_patterns,i in 1..n_windows,t in okTTYPES :
#                      (i,t) in okTourType and p <= num_weekend_patterns[weekend_type[i,t],t]} 
#       >= 0, integer ;


def weekenddaysworked_idx_rule(M):
    index_list =[]
    for (i, t) in M.okTourType:
        for p in pyo.sequence(M.max_weekend_patterns):
            weekendtype = M.weekend_type[i,t]
            if p <= M.num_weekend_patterns[weekendtype,t]:
                index_list.append((i, t, p))
    
    return index_list


model_phase1.weekenddaysworked_idx = pyo.Set(dimen=3, initialize=weekenddaysworked_idx_rule)

model_phase1.WeekendDaysWorked = pyo.Var(model_phase1.weekenddaysworked_idx, within=pyo.NonNegativeIntegers)
    
    # WeekendDaysWorked[d,i,t] = Number of employees working days-off patterns d
    # in start window i and of tour type t


def multiweekpattern_idx_rule(M):
    index_list = []
    for (i,t) in M.okTourType:
        for p1 in pyo.sequence(M.max_mwdw_patterns):
            if p1 <= M.num_mwdw_patterns[t]:
                for p2 in pyo.sequence(M.max_weekend_patterns):
                    weekendtype = M.weekend_type[i, t]
                    if p2 <= M.num_weekend_patterns[weekendtype, t]:
                        index_list.append((i, t, p1, p2))

    return index_list


model_phase1.multiweekpattern_idx = pyo.Set(dimen=4, initialize=multiweekpattern_idx_rule)

model_phase1.MultiWeekPattern = pyo.Var(model_phase1.multiweekpattern_idx, within=pyo.NonNegativeIntegers)

# MultiWeekPattern[i, t, p1, p2] = Number of employees of tour type t working mwdw pattern p1 with
# weekend pattern p2 in start window i

# #### Coverage related pyo.Variables

# For computational convenience, we broke up the calculation of coverage in each period into four separate cases related
# to pyo.Various types of overlap (or lack of). See coverage constraints.

model_phase1.cov1 = pyo.Var(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, within=pyo.NonNegativeReals)
model_phase1.cov2 = pyo.Var(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, within=pyo.NonNegativeReals)
model_phase1.cov3 = pyo.Var(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, within=pyo.NonNegativeReals)
model_phase1.cov4 = pyo.Var(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, within=pyo.NonNegativeReals)
model_phase1.cov = pyo.Var(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, within=pyo.NonNegativeReals)


def under1_bounds(M, i, j, w):
    lb = 0.0
    ub = M.usb.value
    return (lb, ub)


model_phase1.under1 = pyo.Var(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, bounds=under1_bounds)


def under2_bounds(M, i, j, w):
    lb = 0.0
    ub = infinity
    return (lb, ub)

model_phase1.under2 = pyo.Var(model_phase1.PERIODS,model_phase1.DAYS,model_phase1.WEEKS, bounds=under2_bounds)

# Objective function

def objective_rule(M):
    obj1 = sum(M.Shift[i,j,w,k,t] * M.lengths[k] * M.tt_cost_multiplier[t] for (i, j, w, k, t) in M.okShifts)
    obj2 = sum(M.under1[i,j,w] * M.cu1.value + M.under2[i,j,w] * M.cu2.value for (i,j,w) in M.bins)
    # obj3 = sum(M.WeekendDaysWorked[i, t, p] * M.A_tot_wkend_days[p, t, e] ** 2 for (i, t, p) in M.weekenddaysworked_idx for e in M.WEEKENDS if M.weekend_type[i, t] == e)
    return obj1 + obj2


model_phase1.total_cost = pyo.Objective(rule=objective_rule, sense=pyo.minimize)


##### Budget constraints

def max_labor_budget_rule(M): 
    return sum(M.Shift[i,j,w,k,t] * M.lengths[k] * M.tt_cost_multiplier[t] for (i, j, w, k, t) in M.okShifts) <= M.labor_budget

model_phase1.max_labor_budget = pyo.Constraint(rule=max_labor_budget_rule)

# #### Coverage constraints

# Breaking them up into four different constraints, one for each case
# in terms of handling end of day horizon wrapping


def coverage1_rule(M, i, j, w):
    return sum(M.Shift[(i - p), j, w, l, t] for l in M.LENGTHS
               for t in M.activeTT
               for p in range(0, M.lengths[l])
               if (i - p) > 0 and M.allow_start[(i - p), j, l, t] > 0) == M.cov1[i, j, w]


model_phase1.coverage1 = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS,
                                        rule=coverage1_rule)


def coverage2_rule(M, i, j, w):
    return sum(M.Shift[M.n_prds_per_day.value + (i - p), j - 1, w, l, t] for l in M.LENGTHS
               for t in M.activeTT
               for p in range(0, M.lengths[l])
               if (i - p) <= 0 and j > 1 and
               M.allow_start[M.n_prds_per_day.value + (i - p), j - 1, l, t] > 0) == M.cov2[i, j, w]


model_phase1.coverage2 = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS,
                                        rule=coverage2_rule)


def coverage3_rule(M, i, j, w):
    return sum(M.Shift[M.n_prds_per_day.value + (i - p), M.n_days_per_week.value, w - 1, l, t] for l in M.LENGTHS
               for t in M.activeTT
               for p in range(0, M.lengths[l])
               if (i - p) <= 0 and j == 1 and w > 1 and
               M.allow_start[M.n_prds_per_day.value + (i - p), M.n_days_per_week.value, l, t] > 0) == M.cov3[i, j, w]


model_phase1.coverage3 = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS,
                                        rule=coverage3_rule)


def coverage4_rule(M, i, j, w):
    return sum(
        M.Shift[M.n_prds_per_day.value + (i - p), M.n_days_per_week.value, M.n_weeks.value, l, t] for l in M.LENGTHS
        for t in M.activeTT
        for p in range(0, M.lengths[l])
        if (i - p) <= 0 and j == 1 and w == 1 and
        M.allow_start[M.n_prds_per_day.value + (i - p), M.n_days_per_week.value, l, t] > 0) == M.cov4[i, j, w]


model_phase1.coverage4 = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS,
                                        rule=coverage4_rule)


def tot_coverage_rule(M, i, j, w):
    return M.cov1[i, j, w] + M.cov2[i, j, w] + M.cov3[i, j, w] + M.cov4[i, j, w] == M.cov[i, j, w]


model_phase1.tot_coverage = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, rule=tot_coverage_rule)


def coverage_rule(M, i, j, w):
    return M.cov[i, j, w] + M.under1[i, j, w] + M.under2[i, j, w] >= M.dmd_staff[i, j, w]


model_phase1.coverage = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, rule=coverage_rule)


def minstaff_rule(M, i, j, w):
    return M.cov[i, j, w] >= M.min_staff[i, j, w]


model_phase1.minstaff = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, rule=minstaff_rule)

# Lower and upper bounds on tour type variables


def TourType_LB_rule(M, t):
    return sum(M.TourType[i, t] for (i, s) in M.okTourType if s == t) >= M.tt_lb[t]


model_phase1.TourType_LB_con = pyo.Constraint(model_phase1.activeTT, rule=TourType_LB_rule)


def TourType_UB_rule(M, t):
    return sum(M.TourType[i, t] for (i, s) in M.okTourType if s == t) <= M.tt_ub[t]


model_phase1.TourType_UB_con = pyo.Constraint(model_phase1.activeTT, rule=TourType_UB_rule)


# Each tour Variable must get a mwdw pattern and weekday pattern assigned to it


def MWP_total_idx_rule(M):
    """
    The index for the MWP constraints is (window, ttype) tuples
    in the set okTourType.
    """
    return [(i, t) for (i, t) in M.okTourType]


def MWP_total_rule(M, i, t):
    """
    Each tour Variable must get a multi-week pattern assigned to it
    """
    weekendtype = M.weekend_type[i, t]
    return sum(M.MultiWeekPattern[i, t, p1, p2] for p1 in pyo.sequence(M.num_mwdw_patterns[t]) for p2
        in pyo.sequence(M.num_weekend_patterns[weekendtype, t])) == M.TourType[i,t]


model_phase1.MWP_total_idx = pyo.Set(dimen=2, initialize=MWP_total_idx_rule)

model_phase1.MWP_total_con = pyo.Constraint(model_phase1.MWP_total_idx, rule=MWP_total_rule)


# Coordinate weekend days worked variables and multi-week pattern variables.
# Could actually eliminate the weekend variables but convenient for use in constraints.

def MWP_weekend_integration_rule(M, i, t, p):
    """
    Each weekend variable must be consistent with multi-week pattern variables
    """
    return sum(M.MultiWeekPattern[i, t, p1, p] for p1 in pyo.sequence(M.num_mwdw_patterns[t])) \
           == M.WeekendDaysWorked[i, t, p]


model_phase1.MWP_weekend_integration_con = pyo.Constraint(model_phase1.weekenddaysworked_idx,
                                                          rule=MWP_weekend_integration_rule)

# Integrate WeekendDaysWorked variables with DailyTourType variables

def weekend_integration_1_SS_idx_rule(M):
    """
    The number of people working on Sun or Sat as part of a type 1 weekend worked pattern must be equal to the
    number of people working on Sun or Sat as specified by DailyTourType variables.
    """
    return [(j, w, i, t) for j in M.DAYS
            for w in M.WEEKS
            for i in M.WINDOWS
            for t in M.activeTT
            if j in M.weekend[i, t] and 1 in M.weekend[i, t] and 7 in M.weekend[i, t]
            and (i, t, j) in M.okDailyTourType]


model_phase1.weekend_integration_1_SS_idx = pyo.Set(dimen=4, initialize=weekend_integration_1_SS_idx_rule)


def weekend_integration_1_SS_rule(M, j, w, i, t):
    """
    The number of people working on Sun or Sat as part of a type 1 weekend worked pattern must be equal to the
    number of people working on Sun or Sat as specified by DailyTourType pyo.Variables.
    """
    return sum(M.A[p, j, w, t, 1] * M.WeekendDaysWorked[i, t, p] for p in pyo.sequence(M.num_weekend_patterns[1, t])) \
           == M.DailyTourType[i, t, j, w]


model_phase1.weekend_integration_1_SS_con = pyo.Constraint(model_phase1.weekend_integration_1_SS_idx,
                                                           rule=weekend_integration_1_SS_rule)


def weekend_integration_2_FS_rule(M, j, w, i, t):
    """
    The number of people working on Fri or Sat as part of a type 2 weekend worked pattern must be equal to the
    number of people working on Sun or Sat as specified by DailyTourType pyo.Variables.
    """
    return sum(M.A[p, j, w, t, 2] * M.WeekendDaysWorked[i, t, p] for p in pyo.sequence(M.num_weekend_patterns[2, t])) \
           == M.DailyTourType[i, t, j, w]


def weekend_integration_2_FS_idx_rule(M):
    return [(j, w, i, t) for j in M.DAYS
            for w in M.WEEKS
            for i in M.WINDOWS
            for t in M.activeTT
            if j in M.weekend[i, t] and 6 in M.weekend[i, t] and 7 in M.weekend[i, t]
            and (i, t, j) in M.okDailyTourType]


model_phase1.weekend_integration_2_FS_idx = pyo.Set(dimen=4, initialize=weekend_integration_2_FS_idx_rule)

model_phase1.weekend_integration_2_FS_con = pyo.Constraint(model_phase1.weekend_integration_2_FS_idx,
                                                           rule=weekend_integration_2_FS_rule)

# Integrate shift, days worked, and tour type variables

# The index set gets used in the DTT_DSW_con and DTT_TT_UB

def shift_DTT_dailyconservation_idx_rule(M):
    return [(i, j, w, t) for i in M.WINDOWS
            for j in M.DAYS
            for w in M.WEEKS
            for t in M.activeTT
            if (i, t, j) in M.okDailyTourType]


model_phase1.shift_DTT_dailyconservation_idx = pyo.Set(dimen=4, initialize=shift_DTT_dailyconservation_idx_rule)

def DTT_TT_UB_rule(M, i, j, w, t):
    """
    Every day of every week for each (window,ttype), there can be no more people scheduled (DailyTourType) than
    number of people assigned to TourType[window,ttype]
    """
    return M.DailyTourType[i, t, j, w] <= M.TourType[i, t]
    
 
model_phase1.DTT_TT_UB = \
    pyo.Constraint(model_phase1.shift_DTT_dailyconservation_idx, rule=DTT_TT_UB_rule)


# Coordinate DailyShiftWorked and DailyTourType variables for each day of each week

def DTT_DSW_rule(M, i, j, w, t):
    return sum(M.DailyShiftWorked[i, t, k, j, w] for k in M.tt_length_x[t]) == M.DailyTourType[i, t, j, w]


model_phase1.DTT_DSW_con = pyo.Constraint(model_phase1.shift_DTT_dailyconservation_idx, rule=DTT_DSW_rule)

# ---- Min and max bounds on days worked each week

# TODO - actually need to think through possible redundancy of the weekly and cumulative weekly
# constraints on DTT and TT now that the mwdw vars explicitly determine the number of people
# working each day - i.e. the DTT values. Still need these constraints on DSW_TT since those
# are shift length specific variables.

# Index for both lower and upper bound versions of these constraints are the same


def DTT_TT_weeklyconservation_idx_rule(M):
    return [(i,t,w) for (i,t) in M.okTourType
                        for w in M.WEEKS]


model_phase1.DTT_TT_weeklyconservation_idx = pyo.Set(dimen=3, initialize=DTT_TT_weeklyconservation_idx_rule)

# Now that mwdw vars added, we can directly determine DTT sum over each week instead of using the LB, UBs.


def DTT_MWP_idx_rule(M):
    index_list = []
    for t in M.activeTT:
        numpats = M.num_mwdw_patterns[t]
        if numpats > 0:
            for i in M.WINDOWS:
                if (i, t) in M.okTourType:
                    for w in M.WEEKS:
                        index_list.append((i, t, w))
    return index_list


def DTT_MWP_rule(M, i, t, w):
    weekendtype = M.weekend_type[i, t]
    return sum(M.DailyTourType[i, t, j, w] for j in M.DAYS) == sum(
        M.MultiWeekPattern[i, t, p1, p2] * M.A_mwdw[t, p1, w] for p1 in pyo.sequence(M.num_mwdw_patterns[t]) for p2
        in pyo.sequence(M.num_weekend_patterns[weekendtype, t]))


model_phase1.DTT_MWP_idx = pyo.Set(dimen=3, initialize=DTT_MWP_idx_rule)

model_phase1.DTT_MWP_con = pyo.Constraint(model_phase1.DTT_MWP_idx,
                                                      rule=DTT_MWP_rule)


# Lower bound on DTT vars based on minimum number of days worked per week

def DTT_TT_weeklyconservation_LB_rule(M, i, t, w):
    return sum(M.DailyTourType[i, t, d, w] for d in M.DAYS) >= M.TourType[i, t] * M.tt_min_dys_weeks[t, w]


model_phase1.DTT_TT_weeklyconservation_LB = \
  pyo.Constraint(model_phase1.DTT_TT_weeklyconservation_idx, rule=DTT_TT_weeklyconservation_LB_rule)


# Upper bound on DTT vars based on maximum number of days worked per week

def DTT_TT_weeklyconservation_UB_rule(M, i, t, w):
    return sum(M.DailyTourType[i, t, d, w] for d in M.DAYS) <= M.TourType[i, t] * M.tt_max_dys_weeks[t, w]


model_phase1.DTT_TT_weeklyconservation_UB = \
    pyo.Constraint(model_phase1.DTT_TT_weeklyconservation_idx, rule=DTT_TT_weeklyconservation_UB_rule)


# Cumulative (over weeks) versions of the lower and upper bound constraints immediately above

def DTT_TT_cumul_weeklyconservation_LB_rule(M, i, t, w):
    return sum(M.DailyTourType[i, t, d, z] for d in M.DAYS for z in pyo.sequence(w)) >= \
           M.TourType[i, t] * M.tt_min_cumul_dys_weeks[t, w]


model_phase1.DTT_TT_cumul_weeklyconservation_LB = \
    pyo.Constraint(model_phase1.DTT_TT_weeklyconservation_idx, rule=DTT_TT_cumul_weeklyconservation_LB_rule)


def DTT_TT_cumul_weeklyconservation_UB_rule(M, i, t, w):
    return sum(M.DailyTourType[i, t, d, z] for d in M.DAYS for z in pyo.sequence(w)) <= \
           M.TourType[i, t] * M.tt_max_cumul_dys_weeks[t, w]


model_phase1.DTT_TT_cumul_weeklyconservation_UB = \
    pyo.Constraint(model_phase1.DTT_TT_weeklyconservation_idx, rule=DTT_TT_cumul_weeklyconservation_UB_rule)


# The DSW_TT_weeklyconservation constraints feel redundant given that DSW summed over lengths is
# equal to DTT and we've already got DTT_TT constraints above.
# TODO - test for redundancy. Binary deactivation code already added.



def DSW_TT_weeklyconservation_idx_rule(M):
    return [(i,t,w) for (i,t) in M.okTourType for w in M.WEEKS]


model_phase1.DSW_TT_weeklyconservation_idx = pyo.Set(dimen=3, initialize=DSW_TT_weeklyconservation_idx_rule)


def DSW_TT_weeklyconservation_LB_rule(M, i, t, w):
    return sum(M.DailyShiftWorked[i, t, k, d, w] for d in M.DAYS for k in M.tt_length_x[t] if \
               (i, t, k, d, w) in M.DailyShiftWorked_idx) >= M.TourType[i, t] * M.tt_min_dys_weeks[t, w]


model_phase1.DSW_TT_weeklyconservation_LB = \
    pyo.Constraint(model_phase1.DSW_TT_weeklyconservation_idx, rule=DSW_TT_weeklyconservation_LB_rule)


def DSW_TT_weeklyconservation_UB_rule(M, i, t, w):
    return sum(M.DailyShiftWorked[i, t, k, d, w] for d in M.DAYS for k in M.tt_length_x[t] if \
               (i, t, k, d, w) in M.DailyShiftWorked_idx) <= M.TourType[i, t] * M.tt_max_dys_weeks[t, w]


model_phase1.DSW_TT_weeklyconservation_UB = \
    pyo.Constraint(model_phase1.DSW_TT_weeklyconservation_idx, rule=DSW_TT_weeklyconservation_UB_rule)


def DSW_TT_cumul_weeklyconservation_LB_rule(M, i, t, w):
    return sum(M.DailyShiftWorked[i, t, k, d, z] for d in M.DAYS for k in M.tt_length_x[t] for \
               z in pyo.sequence(w) if (i, t, k, d, z) in M.DailyShiftWorked_idx) >= \
           M.TourType[i, t] * M.tt_min_cumul_dys_weeks[t, w]


model_phase1.DSW_TT_cumul_weeklyconservation_LB = \
    pyo.Constraint(model_phase1.DSW_TT_weeklyconservation_idx, rule=DSW_TT_cumul_weeklyconservation_LB_rule)


def DSW_TT_cumul_weeklyconservation_UB_rule(M, i, t, w):
    return sum(M.DailyShiftWorked[i, t, k, d, z] for d in M.DAYS for k in M.tt_length_x[t] for \
               z in pyo.sequence(w) if (i, t, k, d, z) in M.DailyShiftWorked_idx) <= \
           M.TourType[i, t] * M.tt_max_cumul_dys_weeks[t, w]


model_phase1.DSW_TT_cumul_weeklyconservation_UB = \
    pyo.Constraint(model_phase1.DSW_TT_weeklyconservation_idx, rule=DSW_TT_cumul_weeklyconservation_UB_rule)


# Shiftlen versions of the DTT_TT constraints using the DSW vars.
# These are needed as they are shift length specific.

def DSW_TT_shiftlen_weeklyconservation_idx_rule(M):
    return [(i,t,k,w) for (i,t) in M.okTourType
                        for k in M.tt_length_x[t]
                        for w in M.WEEKS]


model_phase1.DSW_TT_shiftlen_weeklyconservation_idx = pyo.Set(dimen=4,initialize=DSW_TT_shiftlen_weeklyconservation_idx_rule) 


def DSW_TT_shiftlen_weeklyconservation_LB_rule(M,i,t,k,w):
    return sum(M.DailyShiftWorked[i,t,k,d,w] for d in M.DAYS)  >= M.TourType[i,t]*M.tt_shiftlen_min_dys_weeks[t,k,w]


model_phase1.DSW_TT_shiftlen_weeklyconservation_LB = \
  pyo.Constraint(model_phase1.DSW_TT_shiftlen_weeklyconservation_idx, rule=DSW_TT_shiftlen_weeklyconservation_LB_rule)


def DSW_TT_shiftlen_weeklyconservation_UB_rule(M,i,t,k,w):
    return sum(M.DailyShiftWorked[i,t,k,d,w] for d in M.DAYS) <= M.TourType[i,t]*M.tt_shiftlen_max_dys_weeks[t,k,w]


model_phase1.DSW_TT_shiftlen_weeklyconservation_UB_con = \
  pyo.Constraint(model_phase1.DSW_TT_shiftlen_weeklyconservation_idx, rule=DSW_TT_shiftlen_weeklyconservation_UB_rule)


def DSW_TT_shiftlen_cumul_weeklyconservation_LB_rule(M,i,t,k,w):
    return sum(M.DailyShiftWorked[i,t,k,d,z] for d in M.DAYS for z in pyo.sequence(w)) >= \
           M.TourType[i,t]*M.tt_shiftlen_min_cumul_dys_weeks[t,k,w]

model_phase1.DSW_TT_shiftlen_cumul_weeklyconservation_LB = \
  pyo.Constraint(model_phase1.DSW_TT_shiftlen_weeklyconservation_idx, rule=DSW_TT_shiftlen_cumul_weeklyconservation_LB_rule)


def DSW_TT_shiftlen_cumul_weeklyconservation_UB_rule(M,i,t,k,w):
    return sum(M.DailyShiftWorked[i,t,k,d,z] for d in M.DAYS for z in pyo.sequence(w)) <= \
           M.TourType[i,t]*M.tt_shiftlen_max_cumul_dys_weeks[t,k,w]


model_phase1.DSW_TT_shiftlen_cumul_weeklyconservation_UB = \
  pyo.Constraint(model_phase1.DSW_TT_shiftlen_weeklyconservation_idx, rule=DSW_TT_shiftlen_cumul_weeklyconservation_UB_rule)


def prds_worked_weekly_LB_rule(M,i,t,w):
    return sum(M.DailyShiftWorked[i,t,k,d,w]*M.lengths[k] for d in M.DAYS for k in M.tt_length_x[t] if (i,t,k,d,w) in M.DailyShiftWorked_idx) \
        >= M.TourType[i,t]*M.tt_min_prds_weeks[t,w]


model_phase1.prds_worked_weekly_LB = \
  pyo.Constraint(model_phase1.okTourType,model_phase1.WEEKS,rule=prds_worked_weekly_LB_rule)


def prds_worked_weekly_UB_rule(M,i,t,w):
    return sum(M.DailyShiftWorked[i,t,k,d,w]*M.lengths[k] for d in M.DAYS for k in M.tt_length_x[t] \
               if (i,t,k,d,w) in M.DailyShiftWorked_idx) <= M.TourType[i,t]*M.tt_max_prds_weeks[t,w]


model_phase1.prds_worked_weekly_UB = \
  pyo.Constraint(model_phase1.okTourType,model_phase1.WEEKS, rule=prds_worked_weekly_UB_rule)


# # Cumulative versions of the above 2 constraints
def prds_worked_cumul_weekly_LB_rule(M,i,t,w):
    return sum(M.DailyShiftWorked[i,t,k,d,z]*M.lengths[k] for d in M.DAYS for k in M.tt_length_x[t] \
               for z in pyo.sequence(w) if (i,t,k,d,z) in M.DailyShiftWorked_idx) >= \
           M.TourType[i,t]*M.tt_min_cumul_prds_weeks[t,w]


model_phase1.prds_worked_cumul_weekly_LB = \
  pyo.Constraint(model_phase1.okTourType,model_phase1.WEEKS, rule=prds_worked_cumul_weekly_LB_rule)


def prds_worked_cumul_weekly_UB_rule(M,i,t,w):
    return sum(M.DailyShiftWorked[i,t,k,d,z]*M.lengths[k] for d in M.DAYS for k in M.tt_length_x[t] \
               for z in pyo.sequence(w) if (i,t,k,d,z) in M.DailyShiftWorked_idx) <= \
           M.TourType[i,t]*M.tt_max_cumul_prds_weeks[t,w]


model_phase1.prds_worked_cumul_weekly_UB = \
  pyo.Constraint(model_phase1.okTourType,model_phase1.WEEKS, rule=prds_worked_cumul_weekly_UB_rule)


# These should be the shift length specific versions of the above 4 constraints
def prds_worked_shiflen_weekly_idx_rule(M):
    return [(i,t,k,w) for (i,t) in M.okTourType
                        for k in M.tt_length_x[t]
                        for w in M.WEEKS]


model_phase1.prds_worked_shiflen_weekly_idx = pyo.Set(dimen=4,initialize=prds_worked_shiflen_weekly_idx_rule) 


def prds_worked_shiflen_weekly_LB_rule(M,i,t,k,w):
    return sum(M.DailyShiftWorked[i,t,k,d,w]*M.lengths[k] for d in M.DAYS if (i,t,k,d,w) in M.DailyShiftWorked_idx) \
        >= M.TourType[i,t]*M.tt_shiftlen_min_prds_weeks[t,k,w]


model_phase1.prds_worked_shiflen_weekly_LB = \
  pyo.Constraint(model_phase1.prds_worked_shiflen_weekly_idx, rule=prds_worked_shiflen_weekly_LB_rule)


def prds_worked_shiflen_weekly_UB_rule(M,i,t,k,w):
    return sum(M.DailyShiftWorked[i,t,k,d,w]*M.lengths[k] for d in M.DAYS if (i,t,k,d,w) in M.DailyShiftWorked_idx) \
        <= M.TourType[i,t]*M.tt_shiftlen_max_prds_weeks[t,k,w]


model_phase1.prds_worked_shiflen_weekly_UB = \
  pyo.Constraint(model_phase1.prds_worked_shiflen_weekly_idx, rule=prds_worked_shiflen_weekly_UB_rule)


# Cumulative versions of the above 2 constraints
def prds_worked_cumul_shiflen_weekly_LB_rule(M,i,t,k,w):
    return sum(M.DailyShiftWorked[i,t,k,d,z]*M.lengths[k] for d in M.DAYS for z in pyo.sequence(w) \
               if (i,t,k,d,z) in M.DailyShiftWorked_idx) >= \
           M.TourType[i,t]*M.tt_shiftlen_min_cumul_prds_weeks[t,k,w]

model_phase1.prds_worked_cumul_shiflen_weekly_LB = \
  pyo.Constraint(model_phase1.prds_worked_shiflen_weekly_idx, rule=prds_worked_cumul_shiflen_weekly_LB_rule)

def prds_worked_cumul_shiflen_weekly_UB_rule(M,i,t,k,w):
    return sum(M.DailyShiftWorked[i,t,k,d,z]*M.lengths[k] for d in M.DAYS for z in pyo.sequence(w) \
               if (i,t,k,d,z) in M.DailyShiftWorked_idx) <= \
           M.TourType[i,t]*M.tt_shiftlen_max_cumul_prds_weeks[t,k,w]


model_phase1.prds_worked_cumul_shiflen_weekly_UB = \
  pyo.Constraint(model_phase1.prds_worked_shiflen_weekly_idx, rule=prds_worked_cumul_shiflen_weekly_UB_rule)



def prds_FSwkend_weeklyconservation_LB_rule(M,i,t,w):
    return sum(M.DailyShiftWorked[i,t,k,d,w]*M.lengths[k] for j in M.DAYS for k in M.LENGTHS if (i,t,k,d,w) in M.DailyShiftWorked_idx) \
                >= sum(M.WeekendDaysWorked[i,t,p]*M.FriSat_min_dys_weeks[t,p,w]*min(M.lengths[k] for k in M.tt_length_x[t]) for p in M.DAYS)


def prds_FSwkend_cumul_weeklyconservation_LB_rule(M,i,t,w):
    return sum(M.DailyShiftWorked[i,t,k,d,z]*M.lengths[k] for j in M.DAYS for k in M.LENGTHS for z in pyo.sequence(w) if (i,t,k,d,z) in M.DailyShiftWorked_idx) \
                >= sum(M.WeekendDaysWorked[i,t,p]*M.FriSat_min_dys_weeks[t,p,w]*min(M.lengths[k] for k in M.tt_length_x[t]) for p in M.DAYS)


def prds_FSwkend_weeklyconservation_LB_idx_rule(M):
    return [(i,t,w) for (i,t) in M.okTourType
                    for w in M.WEEKS
                    if 6 in M.weekend[i,t] and 7 in M.weekend[i,t]]


model_phase1.prds_FSwkend_weeklyconservation_LB_idx = \
    pyo.Set(dimen=3, initialize=prds_FSwkend_weeklyconservation_LB_idx_rule)


model_phase1.prds_FSwkend_weeklyconservation_LB = \
    pyo.Constraint(model_phase1.prds_FSwkend_weeklyconservation_LB_idx,rule=prds_FSwkend_weeklyconservation_LB_rule)


model_phase1.prds_FSwkend_cumul_weeklyconservation_LB = \
    pyo.Constraint(model_phase1.prds_FSwkend_weeklyconservation_LB_idx,rule=prds_FSwkend_cumul_weeklyconservation_LB_rule)


def DWT_FSwkend_weeklyconservation_LB_rule(M,i,t,w):
    return sum(M.DailyTourType[i,t,d,w] for d in M.DAYS if (i,t,d) in M.okDailyTourType) >= \
        sum(M.WeekendDaysWorked[i,t,p] * M.FriSat_min_dys_weeks[t,p,w] for \
            p in range(1,M.num_weekend_patterns[M.weekend_type[i,t].value,t].value+1))


def DWT_FSwkend_cumul_weeklyconservation_LB_rule(M,i,t,w):
    return sum(M.DailyTourType[i,t,d,z] for d in M.DAYS for z in pyo.sequence(w) if (i,t,d) in M.okDailyTourType) >= \
        sum(M.WeekendDaysWorked[i,t,p] * M.FriSat_min_cumul_dys_weeks[t,p,w] for \
            p in range(1,M.num_weekend_patterns[M.weekend_type[i,t].value,t].value+1))

                                         
def DWT_FSwkend_weeklyconservation_LB_idx_rule(M):
    return [(i,t,w) for (i,t) in M.okTourType
                    for w in M.WEEKS
                    if 6 in M.weekend[i,t] and 7 in M.weekend[i,t]]


model_phase1.DWT_FSwkend_weeklyconservation_LB_idx = pyo.Set(dimen=3,
                                                             initialize=DWT_FSwkend_weeklyconservation_LB_idx_rule)
model_phase1.DWT_FSwkend_weeklyconservation_LB = pyo.Constraint(model_phase1.DWT_FSwkend_weeklyconservation_LB_idx,
                                                                rule=DWT_FSwkend_weeklyconservation_LB_rule)
model_phase1.DWT_FSwkend_cumul_weeklyconservation_LB = pyo.Constraint(
    model_phase1.DWT_FSwkend_weeklyconservation_LB_idx, rule=DWT_FSwkend_cumul_weeklyconservation_LB_rule)


def max_ptfrac_rule(M):
    return sum(M.Shift[i, j, w, k, t] * M.lengths[k] for (i, j, w, k, t) in M.okShifts if M.tt_parttime[t] > 0) <= \
           M.max_parttime_frac.value * sum(M.Shift[i, j, w, k, t] * M.lengths[k] for (i, j, w, k, t) in M.okShifts)


model_phase1.max_ptfrac_con = pyo.Constraint(rule=max_ptfrac_rule)


# Chains - coordinates DWS and Shift within each chain for intra-tour start-time flexibility


# subject to chains_sweep_l{e in WEEKS, t in okTTYPES, k in tt_length_x[t], (b,j) in bchain[t,k],
#             p in period[b,j]..period[b,j],
#             w in 0..(numlinks[t,k,b,j]-1-(p-period[b,j])):
#   width>0 } :
#    sum{ (l,m) in linkspan[t,k,b,j,w+1]:
#     (l,m,k,t) in ok_shifts} Shift[l,m,e,k,t] >=
#     sum{u in p..p+w: (which_prd[u],which_day[u]) in okWindowBeginnings[t,k]
#     and sum{(l,m) in WindowWepochs[which_prd[u],which_day[u]]} allow_start[l,m,k,t]>0} 
#      DailyShiftWorked[which_prd[u],t,k,which_day[u],e] ; 


def chains_sweep_l_rule(M, t, k, b, j, w, p, v):
    # lhs = sum(M.Shift[l,m,n,k,t] for (l,m,n) in M.linkspan[t,k,b,j,w,v+1] if (l,m,n,k,t) in M.okShifts)

    return sum(
        M.Shift[l, m, n, k, t] for (l, m, n) in M.linkspan[t, k, b, j, w, v + 1] if (l, m, n, k, t) in M.okShifts) >= \
           sum(M.DailyShiftWorked[g_prd_to_tuple(M, u)[0], t, k, g_prd_to_tuple(M, u)[1], g_prd_to_tuple(M, u)[2]]
               for u in [vv for vv in range(p, p + M.g_start_window_width + 1)
                         if (g_prd_to_tuple(M, v)[0], g_prd_to_tuple(M, v)[1], g_prd_to_tuple(M, v)[2])
                         in M.okStartWindowRoots[t, k] and sum(M.allow_start[x, y, k, t] for (x, y, z)
                                                               in M.PotentialGlobalStartWindow[
                                                                   g_prd_to_tuple(M, v)[0], g_prd_to_tuple(M, v)[1],
                                                                   g_prd_to_tuple(M, v)[2]]) > 0])


def chains_sweep_l_idx_rule(M):
    index_list = []
    if M.g_start_window_width > 0:
        for t in M.activeTT:
            for k in [len for len in M.LENGTHS if len in M.tt_length_x[t]]:
                for b in M.PERIODS:
                    for j in M.DAYS:
                        for w in M.WEEKS:
                            for p in range(M.g_period[b, j, w].value, M.g_period[b, j, w] + 1):
                                if (b, j, w) in M.bchain[t, k]:
                                    for m in range(0,
                                                   M.n_links[t, k, b, j, w].value - 1 - (p - M.g_period[b, j, w]) + 1):
                                        index_list.append((t, k, b, j, w, p, m))
                                        
    return index_list


model_phase1.chains_sweep_l_idx = pyo.Set(dimen=7, ordered=True, initialize=chains_sweep_l_idx_rule)
model_phase1.chains_sweep_l_con = pyo.Constraint(model_phase1.chains_sweep_l_idx, rule=chains_sweep_l_rule)


# subject to chains_sweep_u{e in WEEKS,t in okTTYPES, k in tt_length_x[t],(b,j) in bchain[t,k],
#             w in 0..(numlinks[t,k,b,j]-1):
#   width>0 } :
#    sum{ i in period[b,j]..period[b,j]+w :
#     (which_prd[i],which_day[i],k,t) in ok_shifts} 
#      Shift[which_prd[i],which_day[i],e,k,t] <=
#     sum{u in period[b,j]..period[b,j]+w : (which_prd[u],which_day[u]) in okWindowBeginnings[t,k] and 
#      sum{(l,m) in WindowWepochs[which_prd[u],which_day[u]]} allow_start[l,m,k,t]>0} 
#      DailyShiftWorked[which_prd[u],t,k,which_day[u],e] ; 

def chains_sweep_u_rule(M, t, k, b, j, w, p, v):
    return sum(M.Shift[g_prd_to_tuple(M, i)[0], g_prd_to_tuple(M, i)[1], g_prd_to_tuple(M, i)[2], k, t]
               for i in range(p, p + v + 1) if
               (g_prd_to_tuple(M, i)[0], g_prd_to_tuple(M, i)[1], g_prd_to_tuple(M, i)[2]) in M.okShifts) <= sum(
        M.DailyShiftWorked[g_prd_to_tuple(M, u)[0], t, k, g_prd_to_tuple(M, u)[1], g_prd_to_tuple(M, u)[2]]
        for u in [vv for vv in range(p, p + M.g_start_window_width + 1)
                  if
                  (g_prd_to_tuple(M, v)[0], g_prd_to_tuple(M, v)[1], g_prd_to_tuple(M, v)[2]) in M.okStartWindowRoots[
                      t, k] \
                  and sum(M.allow_start[x, y, k, t] for (x, y, z) in M.PotentialGlobalStartWindow[
                      g_prd_to_tuple(M, v)[0], g_prd_to_tuple(M, v)[1], g_prd_to_tuple(M, v)[2]]) > 0])


def chains_sweep_u_idx_rule(M):
    index_list = []
    if M.g_start_window_width > 0:
        for t in M.activeTT:
            for k in [len for len in M.LENGTHS if len in M.tt_length_x[t]]:
                for b in M.PERIODS:
                    for j in M.DAYS:
                        for w in M.WEEKS:
                            for p in range(M.g_period[b, j, w].value, M.g_period[b, j, w] + 1):
                                if (b, j, w) in M.bchain[t, k]:
                                    for m in range(0, M.n_links[t, k, b, j, w] - 1):
                                        index_list.append((t, k, b, j, w, p, m))

    return index_list


model_phase1.chains_sweep_u_idx = pyo.Set(dimen=7, ordered=True, initialize=chains_sweep_u_idx_rule)
model_phase1.chains_sweep_u_con = pyo.Constraint(model_phase1.chains_sweep_l_idx, rule=chains_sweep_u_rule)
      
      


#subject to chains_tot{e in WEEKS, t in okTTYPES, k in tt_length_x[t],(i,j) in bchain[t,k]} :
# sum{(l,m) in linkspan[t,k,i,j,numlinks[t,k,i,j]]
#   :(l,m,k,t) in ok_shifts}
#    Shift[l,m,e,k,t] = sum{(n,o) in chain[i,j,t,k]:
#      sum{(p,d) in WindowWepochs[n,o]} allow_start[p,d,k,t]>0
#      } DailyShiftWorked[n,t,k,o,e];
#

def chains_tot_rule(M,t,k,b,j,w):
  # return  sum(M.Shift[l,m,n,k,t] for (l,m,n) in [(p,q,r) for (p,q,r) in M.linkspan[t,k,b,j,w,M.n_links[t,k,b,j,w]] if (p,q,r,k,t) in M.okShifts]) == \
#  return sum(M.Shift[l,m,n,k,t] for (l,m,n) in [(p,q,r) for (p,q,r) in M.linkspan[t,k,b,j,w,M.n_links[t,k,b,j,w]] ]) == \
#             sum(M.DailyShiftWorked[x,t,k,y,z] 
#                   for (x,y,z) in [(u,v,x) for (u,v,x) in M.chain[t,k,b,j,w]
#                     if sum(M.allow_start[p,q,k,t] for (p,q,r) in M.PotentialGlobalStartWindow[u,v,x])>0])
    return sum(M.Shift[l,m,n,k,t] for (l,m,n) in [(p,q,r) for (p,q,r) in M.linkspan[t,k,b,j,w,M.n_links[t,k,b,j,w]] if (p,q,r,k,t) in M.okShifts]) - \
             sum(M.DailyShiftWorked[x,t,k,y,z] 
                   for (x,y,z) in [(u,v,xx) for (u,v,xx) in M.chain[t,k,b,j,w]
                      if sum(M.allow_start[p,q,k,t] for (p,q,r) in M.PotentialGlobalStartWindow[u,v,x])>0]) == 0
                           
def chains_tot_idx_rule(M):
    return [(t,k,b,j,w) for t in M.activeTT
                          for k in [length for length in M.LENGTHS if length in M.tt_length_x[t]]
                          for b in M.PERIODS
                          for j in M.DAYS
                          for w in M.WEEKS
                          if M.g_start_window_width.value > 0 and (b,j,w) in M.bchain[t,k]]

                      
model_phase1.chains_tot_idx = pyo.Set(dimen=5,ordered=True,initialize=chains_tot_idx_rule)  
model_phase1.chains_tot_con = pyo.Constraint(model_phase1.chains_tot_idx,rule=chains_tot_rule)     


# TODO - the follow proxy constraints are only case of no intra-tour start time flexibility
# Coordinate Shift pyo.Variables with DailyShiftWorked for start window width=0.
# In this case, shift periods are same as shiftworked windows

def chains_tot_proxy1_rule(M, w, t, k, i, j):
    return M.Shift[i, j, w, k, t] == M.DailyShiftWorked[i, t, k, j, w]


def chains_tot_proxy1_idx_rule(M):
    return [(w, t, k, i, j) for w in M.WEEKS
            for t in M.activeTT
            for k in M.tt_length_x[t]
            for i in M.PERIODS
            for j in M.DAYS if (i, j, w, k, t) in M.okShifts]
                    

model_phase1.chains_tot_proxy1_idx = pyo.Set(dimen=5, initialize=chains_tot_proxy1_idx_rule)

model_phase1.chains_tot_proxy1_con = pyo.Constraint(model_phase1.chains_tot_proxy1_idx, rule=chains_tot_proxy1_rule)


def chains_tot_proxy2_rule(M, w, t, k, i, j):
    return M.DailyShiftWorked[i, t, k, j, w] == 0


def chains_tot_proxy2_idx_rule(M):
    return [(w, t, k, i, j) for w in M.WEEKS
            for t in M.activeTT
            for k in M.tt_length_x[t]
            for i in M.PERIODS
            for j in M.DAYS
            if (i, t, j) in M.okDailyTourType and M.allow_start[i, j, k, t] == 0]


model_phase1.chains_tot_proxy2_idx = pyo.Set(dimen=5, initialize=chains_tot_proxy2_idx_rule)

model_phase1.chains_tot_proxy2_con = pyo.Constraint(model_phase1.chains_tot_proxy2_idx, rule=chains_tot_proxy2_rule)

# Weekend subset constraints
# Need constraints to prevent cases such as the following. Consider two people with same tour type working
# 5 days per week. Assume that consecutive weekends are allowed and that the two weekend patterns for week 1
# are:  1 0 0 0 0 0 1 and 0 0 0 0 0 0 0. Now consider the following DailyTourType solution:
    
#  1 2 0 2 2 2 1

# The employee with the 1 0 0 0 0 0 1 pattern would be forced to work 6 days.

def weekend_subsets_5_4_idx_rule(M):
    index_list = []
    
    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_max_dys_weeks[t, w] == 5:
                    index_list.append((i, t, w, e, 2, 3, 4, 5))
                    index_list.append((i, t, w, e, 2, 3, 4, 6))
                    index_list.append((i, t, w, e, 2, 3, 5, 6))
                    index_list.append((i, t, w, e, 2, 4, 5, 6))
                    index_list.append((i, t, w, e, 3, 4, 5, 6))
        
    return index_list


model_phase1.weekend_subsets_5_4_idx = pyo.Set(dimen=8, initialize=weekend_subsets_5_4_idx_rule)


def weekend_subsets_5_4_rule(M, i, t, w, e, d1, d2, d3, d4):

    days_subset = [d1, d2, d3, d4]
    return sum(M.DailyTourType[i, t, d, w] for d in days_subset) <= \
           sum(M.MultiWeekPattern[i, t, p1, p2] * min(len(days_subset), (M.A_mwdw[t, p1, w] - M.A_num_wkend_days[p2, w, t, e])) for p1 in pyo.sequence(M.num_mwdw_patterns[t])
               for p2 in pyo.sequence(M.num_weekend_patterns[e, t]))


model_phase1.weekend_subsets_5_4_con2 = pyo.Constraint(model_phase1.weekend_subsets_5_4_idx,
                                                       rule=weekend_subsets_5_4_rule)

def weekend_subsets_4_3_idx_rule(M):
    index_list = []
    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] <= 4 <= M.tt_max_dys_weeks[t, w]:
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


model_phase1.weekend_subsets_4_3_idx = pyo.Set(dimen=7, initialize=weekend_subsets_4_3_idx_rule)


def weekend_subsets_4_3_rule(M, i, t, w, e, d1, d2, d3):

    days_subset = [d1, d2, d3]
    return sum(M.DailyTourType[i, t, d, w] for d in days_subset) <= \
           sum(M.MultiWeekPattern[i, t, p1, p2] * min(len(days_subset), (M.A_mwdw[t, p1, w] - M.A_num_wkend_days[p2, w, t, e])) for p1 in pyo.sequence(M.num_mwdw_patterns[t])
               for p2 in pyo.sequence(M.num_weekend_patterns[e, t]))


model_phase1.weekend_subsets_4_3_con2 = pyo.Constraint(model_phase1.weekend_subsets_4_3_idx,
                                                       rule=weekend_subsets_4_3_rule)


def weekend_subsets_3_2_idx_rule(M):
    index_list = []

    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] <= 3 <= M.tt_max_dys_weeks[t, w]:
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


model_phase1.weekend_subsets_3_2_idx = pyo.Set(dimen=6, initialize=weekend_subsets_3_2_idx_rule)


def weekend_subsets_3_2_rule(M, i, t, w, e, d1, d2):

    days_subset = [d1, d2]
    return sum(M.DailyTourType[i, t, d, w] for d in days_subset) <= \
           sum(M.MultiWeekPattern[i, t, p1, p2] * min(len(days_subset), (M.A_mwdw[t, p1, w] - M.A_num_wkend_days[p2, w, t, e])) for p1 in pyo.sequence(M.num_mwdw_patterns[t])
               for p2 in pyo.sequence(M.num_weekend_patterns[e, t]))


model_phase1.weekend_subsets_3_2_con2 = pyo.Constraint(model_phase1.weekend_subsets_3_2_idx,
                                                       rule=weekend_subsets_3_2_rule)


def weekend_subsets_2_1_idx_rule(M):
    index_list = []

    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] <= 2 <= M.tt_max_dys_weeks[t, w]:
                    index_list.append((i, t, w, e, 2))
                    index_list.append((i, t, w, e, 3))
                    index_list.append((i, t, w, e, 4))
                    index_list.append((i, t, w, e, 5))
                    index_list.append((i, t, w, e, 6))

    return index_list


model_phase1.weekend_subsets_2_1_idx = pyo.Set(dimen=5, initialize=weekend_subsets_2_1_idx_rule)


def weekend_subsets_2_1_rule(M, i, t, w, e, d1):
    days_subset = [d1]
    return sum(M.DailyTourType[i, t, d, w] for d in days_subset) <= \
           sum(M.MultiWeekPattern[i, t, p1, p2] * min(len(days_subset), (M.A_mwdw[t, p1, w] - M.A_num_wkend_days[p2, w, t, e])) for p1 in pyo.sequence(M.num_mwdw_patterns[t])
               for p2 in pyo.sequence(M.num_weekend_patterns[e, t]))


model_phase1.weekend_subsets_2_1_con2 = pyo.Constraint(model_phase1.weekend_subsets_2_1_idx,
                                                       rule=weekend_subsets_2_1_rule)


def main():
    pass


if __name__ == '__main__':
    main()
