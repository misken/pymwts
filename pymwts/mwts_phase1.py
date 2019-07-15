#-------------------------------------------------------------------------------
# Name:        mwts_phase1.py
# Purpose:     Phase 1 for implicit multi-week tour scheduling model
#
# Author:      isken
#
# Created:     05/31/2011
# Copyright:   (c) isken 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import pyomo.environ as pyo
#from pyomo.environ import *

#from pyutilib.misc import import_file
#from pymwts.mwts_utils import *

import itertools


# TODO Would be nice if Phase 1 and Phase 2 could share base pyo.Parameters
# model_phase1 = import_file('mwts_basepyo.Params.py').model


model_phase1 = pyo.AbstractModel()
model_phase1.name = "mwts_phase1"


# Constants 

infinity = float('inf')

# General temporal pyo.Parameters

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

model_phase1.n_prds_per_week = pyo.Param(within=pyo.PositiveIntegers, initialize=n_prds_per_week_init)


def n_prds_per_cycle_init(M):
    """
    Initialize convenience pyo.Parameter n_prds_per_cycle where cycle may include
    one or more weeks.
    """
    return M.n_weeks()*M.n_days_per_week()*M.n_prds_per_day()


model_phase1.n_prds_per_cycle = \
    pyo.Param(within=pyo.PositiveIntegers, initialize=n_prds_per_cycle_init)


def g_period_init(M, i, j, w):
    return ((w - 1) * M.n_days_per_week() * M.n_prds_per_day() +
            (j - 1) * M.n_prds_per_day() + i)


def g_tuple_to_period(i, j, w, n_prds_per_day, n_days_per_week):
    return ((w - 1) * n_days_per_week * n_prds_per_day +
            (j - 1) * n_prds_per_day + i)


def period_increment(M, i, j, w, incr):
    p = M.g_period[i, j, w]
    if (p + incr <= M.n_prds_per_cycle):
        return p + incr
    else:
        return p + incr - M.n_prds_per_cycle


def g_period_increment(M, p, incr):
    if p + incr <= M.n_prds_per_cycle:
        return p + incr
    else:
        return p + incr - M.n_prds_per_cycle


def g_period_difference(M, b_prd, e_prd):
    if e_prd >= b_prd:
        return e_prd - b_prd + 1
    else:
        return M.n_prds_per_cycle + e_prd - b_prd + 1


def g_prd_to_tuple(M, p):
    #    Param which_prd{p in 1..(n_days+1)*n_prds_per_day} :=
    #   p-n_prds_per_day*(ceil(p/n_prds_per_day-1));
    #
    # Param which_day{p in 1..(n_days+1)*n_prds_per_day} :=
    #   (if p>n_prds_per_day*n_days then 1 else 1+ceil(p/n_prds_per_day-1));

    n_week = ((p - 1) // M.n_prds_per_week.value) + 1
    prds_remainder = p - (n_week - 1) * M.n_prds_per_week
    if prds_remainder == 0:
        n_day = 1
    else:
        n_day = ((prds_remainder - 1) // M.n_prds_per_day.value) + 1

    prds_remainder = prds_remainder - (n_day - 1) * M.n_prds_per_day
    if prds_remainder == 0:
        n_period = 1
    else:
        n_period = prds_remainder

    return n_period, n_day, n_week

# For range sets, if start omitted, assumed range is 1..args[0]
model_phase1.PERIODS = pyo.RangeSet(1,model_phase1.n_prds_per_day)  # P
model_phase1.CYCLEPERIODS = pyo.RangeSet(1,model_phase1.n_prds_per_cycle)  # P
model_phase1.WINDOWS = pyo.RangeSet(1,model_phase1.n_prds_per_day)  
model_phase1.DAYS = pyo.RangeSet(1,model_phase1.n_days_per_week)    # D
model_phase1.WEEKS = pyo.RangeSet(1,model_phase1.n_weeks)           # W
model_phase1.WEEKENDS = pyo.RangeSet(1,2)

model_phase1.g_period = pyo.Param(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, initialize=g_period_init)

model_phase1.bins = model_phase1.PERIODS * model_phase1.DAYS * model_phase1.WEEKS  # B

def oneweek_bins_init(M):
    """
    Initialize convenience set oneweek_bins. (period,day) pairs. 
    
    Does not appear to be used.
    """
    return [(i,j) for i in M.PERIODS
                  for j in M.DAYS]
    
model_phase1.oneweek_bins = pyo.Set(dimen=2, ordered=True, initialize=oneweek_bins_init)


#### Tour type related pyo.Parameters

#-- Shift Lengths
model_phase1.n_lengths = pyo.Param(within=pyo.PositiveIntegers)    # Number of shift lengths
model_phase1.LENGTHS = pyo.RangeSet(1,model_phase1.n_lengths)
model_phase1.lengths = pyo.Param(model_phase1.LENGTHS)  # Vector of shift lengths

#-- Tour Types
model_phase1.n_tts = pyo.Param(within=pyo.PositiveIntegers)  # Number of different tour types
model_phase1.TTYPES = pyo.RangeSet(1,model_phase1.n_tts)
model_phase1.tt_length_x = pyo.Set(model_phase1.TTYPES,ordered=True,)  # Set of allowable length indices by tour type

### Bounds on tour type pyo.Variables
model_phase1.tt_lb =  pyo.Param(model_phase1.TTYPES)       # RHS from .MIX
model_phase1.tt_ub =  pyo.Param(model_phase1.TTYPES,  default=infinity)

#-- Weekend patterns

def maxwkend_init(M):
    
    maxpats = 2**(2*M.n_weeks.value)
    return maxpats

model_phase1.max_weekend_patterns = pyo.Param(initialize=maxwkend_init)

model_phase1.num_weekend_patterns = pyo.Param(model_phase1.WEEKENDS,model_phase1.TTYPES)       # Number of weekends worked patterns


## pyo.Param A[i,j,w,t,e] = 1 if weekend pattern i calls for work on day j of week k for tour type t having weekend type e and 0 otherwise

def A_idx_rule(M):
    return [(i,j,w,t,e) for i in pyo.sequence(M.max_weekend_patterns)
                 for j in M.DAYS
                 for w in M.WEEKS
                 for t in M.TTYPES
                 for e in pyo.sequence(2) if i <= M.num_weekend_patterns[e,t]]

model_phase1.A_idx = pyo.Set(dimen=5,ordered=True,initialize=A_idx_rule)
model_phase1.A = pyo.Param(model_phase1.A_idx, within=pyo.Boolean, default=0)

def A_wkend_days_idx_rule(M):
    return [(i,w,t,e) for i in pyo.sequence(M.max_weekend_patterns)
                 for w in M.WEEKS
                 for t in M.TTYPES
                 for e in pyo.sequence(2) if i <= M.num_weekend_patterns[e,t]]

model_phase1.A_wkend_days_idx = pyo.Set(dimen=4,ordered=True,initialize=A_wkend_days_idx_rule)


# -- Multiweek days worked patterns

def maxmwdw_init(M):
    maxmwdw = 4 ** M.n_weeks.value
    return maxmwdw


model_phase1.max_mwdw_patterns = pyo.Param(initialize=maxmwdw_init)

model_phase1.num_mwdw_patterns = pyo.Param(model_phase1.TTYPES)  # Number of mwdw worked patterns


def A_mwdw_idx_rule(M):
    return [(t, p, w) for t in M.TTYPES
            for p in pyo.sequence(M.max_mwdw_patterns)
            for w in M.WEEKS
            if p <= M.num_mwdw_patterns[t]]


model_phase1.A_mwdw_idx = pyo.Set(dimen=3, ordered=True, initialize=A_mwdw_idx_rule)
model_phase1.A_mwdw = pyo.Param(model_phase1.A_mwdw_idx, within=pyo.NonNegativeIntegers, default=0)


def activeTT_init(M):
    return [t for t in M.TTYPES if M.tt_ub[t] > 0]


model_phase1.activeTT = pyo.Set(dimen=1, ordered=True, initialize=activeTT_init)

def A_num_wkend_days_init(M,i,w,t,e):
    """
    Number of weekend days worked in the i'th pattern, for week k, ttype t,
    and weekend type e.
    """
    if e == 1:
        return M.A[i,1,w,t,e] + M.A[i,7,w,t,e]

            
    else:
        return M.A[i,6,w,t,e] + M.A[i,7,w,t,e]
        
    
model_phase1.A_num_wkend_days = pyo.Param(model_phase1.A_wkend_days_idx,initialize=A_num_wkend_days_init)


def A_tot_wkend_days_idx_rule(M):
    return [(i, t, e) for i in pyo.sequence(M.max_weekend_patterns)
                 for t in M.TTYPES
                 for e in pyo.sequence(2) if i <= M.num_weekend_patterns[e,t]]


model_phase1.A_tot_wkend_days_idx = pyo.Set(dimen=3, ordered=True, initialize=A_tot_wkend_days_idx_rule)


def A_tot_wkend_days_init(M, i, t, e):
    """
    Number of weekend days worked in the i'th pattern, ttype t,
    and weekend type e.
    """
    if e == 1:
        return sum(M.A[i, 1, w, t, e] + M.A[i, 7, w, t, e] for w in M.WEEKS)
    else:
        return sum(M.A[i, 6, w, t, e] + M.A[i, 7, w, t, e] for w in M.WEEKS)


model_phase1.A_tot_wkend_days = pyo.Param(model_phase1.A_tot_wkend_days_idx, initialize=A_tot_wkend_days_init)


def A_is_two_wkend_days_init(M,i,w,t,e):
    """
    Initialize indicator for each weekend worked pattern as to whether each week is
    sandwiched by consec weekend days worked.
    """
    if e == 1:
        if M.A[i,1,w,t,e] == 1 and M.A[i,7,w,t,e] == 1:
            return 1
        else:
            return 0
            
    else:
        if M.A[i,6,w,t,e] == 1 and M.A[i,7,w,t,e] == 1:
            return 1
        else:
            return 0
        
    
model_phase1.A_is_two_wkend_days = pyo.Param(model_phase1.A_wkend_days_idx,initialize=A_is_two_wkend_days_init)


def A_is_one_wkend_days_init(M,i,w,t,e):
    """
    Initialize indicator for each weekend worked pattern as to whether each week is
    sandwiched by consec weekend days worked.
    """
    if e == 1:
        if (M.A[i,1,w,t,e] + M.A[i,7,w,t,e]) == 1:
            return 1
        else:
            return 0
            
    else:
        if (M.A[i,6,w,t,e] + M.A[i,7,w,t,e]) == 1:
            return 1
        else:
            return 0
        
    
model_phase1.A_is_one_wkend_days = pyo.Param(model_phase1.A_wkend_days_idx,initialize=A_is_one_wkend_days_init)

def A_is_Sunday_init(M,i,w,t,e):
    """
    Initialize indicator for each weekend worked pattern as to whether it's
    a Sunday only pattern
    """
    if e == 1:
        if M.A[i,1,w,t,e] == 1 and M.A[i,7,w,t,e] == 0 :
            return 1
        else:
            return 0
            
    else:
        # TODO - FS weekend not implemented
        if M.A[i,1,w,t,e] == 1 and M.A[i,7,w,t,e] == 0 :
            return 1
        else:
            return 0
        
    
model_phase1.A_is_Sunday = pyo.Param(model_phase1.A_wkend_days_idx,initialize=A_is_Sunday_init)


def A_is_Saturday_init(M, i, w, t, e):
    """
    Initialize indicator for each weekend worked pattern as to whether it's
    a Saturday only pattern
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


model_phase1.A_is_Saturday = pyo.Param(model_phase1.A_wkend_days_idx, initialize=A_is_Saturday_init)


# Weekends consisting of a Fri and Sat imply updated lower bounds on some of the daily tour type pyo.Variables. 
# These were the key to modeling weekends worked patterns.

def FriSat_idx_rule(M):
    # return [(t,e,w) for t in M.activeTT for e in range(1,M.num_weekend_patterns[2,t] + 1) for w in model_phase1.WEEKS]
    index_list = []
    for t in M.activeTT:
        numpats = M.num_weekend_patterns[2, t]
        for p in pyo.sequence(numpats):
            for w in M.WEEKS:
                index_list.append((t, p, w))

    return index_list


model_phase1.FriSat_idx = pyo.Set(dimen=3, ordered=True, initialize=FriSat_idx_rule)

model_phase1.FriSat_min_dys_weeks = pyo.Param(model_phase1.FriSat_idx, default=0)
model_phase1.FriSat_min_cumul_dys_weeks = pyo.Param(model_phase1.FriSat_idx, default=0)

# Bounds on days and shifts worked over the week

# 1a. Min and max number of days worked by week by tour type
model_phase1.tt_min_dys_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.WEEKS, default=0.0)       
model_phase1.tt_max_dys_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.WEEKS, default=1e+6)         

# 1b. Min and max number of days worked by cumulative weeks by tour type
model_phase1.tt_min_cumul_dys_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.WEEKS, default=0.0) 
model_phase1.tt_max_cumul_dys_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.WEEKS, default=1e+6) 

# 2a. Min and max number of days worked by week by shiftlen by tour type
model_phase1.tt_shiftlen_min_dys_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.LENGTHS, model_phase1.WEEKS, default=0)        
model_phase1.tt_shiftlen_max_dys_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.LENGTHS, model_phase1.WEEKS, default=1e+6)    

# 2b. Min and max number of days worked by cumulative weeks by shiftlen by tour type
model_phase1.tt_shiftlen_min_cumul_dys_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.LENGTHS, model_phase1.WEEKS, default=0)         
model_phase1.tt_shiftlen_max_cumul_dys_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.LENGTHS, model_phase1.WEEKS, default=1e+6) 

# 3a. Min and max number of periods worked by week by tour type
model_phase1.tt_min_prds_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.WEEKS, default=0)        
model_phase1.tt_max_prds_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.WEEKS, default=1e+6) 

# 3b. Min and max number of periods worked by cumulative weeks by tour type
model_phase1.tt_min_cumul_prds_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.WEEKS, default=0)         
model_phase1.tt_max_cumul_prds_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.WEEKS, default=1e+6)   

# 4a. Min and max number of periods worked by week by tour type by shift length
model_phase1.tt_shiftlen_min_prds_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.LENGTHS, model_phase1.WEEKS, default=0)        
model_phase1.tt_shiftlen_max_prds_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.LENGTHS, model_phase1.WEEKS, default=1e+6)       

# 4b. Min and max number of periods worked by cumulative weeks by tour type by shift length
model_phase1.tt_shiftlen_min_cumul_prds_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.LENGTHS, model_phase1.WEEKS, default=0)         
model_phase1.tt_shiftlen_max_cumul_prds_weeks = pyo.Param(model_phase1.TTYPES, model_phase1.LENGTHS, model_phase1.WEEKS, default=1e+6)        


# To the above, we'll add pyo.Params and sets to allow direct modeling of side constraints
# of the form sum{subset of tour types} =, >=, <= some bound

model_phase1.tt_parttime = pyo.Param(model_phase1.TTYPES)     # 1 for part-time, 0 for full-time


# Allowable shift start times  - note that these are tour type specific
model_phase1.allow_start = pyo.Param(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.LENGTHS, model_phase1.TTYPES,
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
# COST RELATED pyo.ParamETERS
# -----------------------------------------------------------------------

model_phase1.tt_cost_multiplier = pyo.Param(model_phase1.TTYPES)           # Tour type differential

model_phase1.cu1 = pyo.Param()
model_phase1.cu2 = pyo.Param()        
model_phase1.usb = pyo.Param()

# -----------------------------------------------------------------------
# Weekend RELATED pyo.ParamETERS
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

model_phase1.b_window_epoch = pyo.Param(model_phase1.bins,initialize=b_window_epoch_init)


def e_window_epoch_init(M,i,j,w):
    epoch = M.n_days_per_week.value * M.n_prds_per_day.value*(w-1) + \
            M.n_prds_per_day.value * (j-1) + i + M.g_start_window_width.value

    if epoch <= M.n_prds_per_cycle.value:
        e_window = epoch
    else:
        e_window = epoch-M.n_prds_per_cycle.value

    return e_window

model_phase1.e_window_epoch = pyo.Param(model_phase1.bins,initialize=e_window_epoch_init)



###### PotentialGlobalStartWindow ######

# Index: PotentialGlobalStartWindow is defined for all (period,day,week) bins.
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
    # Compute global period of (i,j,w)
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
    # to convert it back to a tuple from a global period number.


    window_list.append(prd)
    return window_list


model_phase1.echain = pyo.Set(model_phase1.chain_idx,ordered=True,dimen=3,initialize=echain_init)


def n_links_init(M,t,k,i,j,w):

    # Compute global period of (i,j,w)
    b_g_prd = M.g_period[i, j, w]
    e_prd = M.echain[t, k, i, j, w][1]
    e_g_prd = M.g_period[e_prd[0], e_prd[1], e_prd[2]]

    return g_period_difference(M,b_g_prd,e_g_prd)


model_phase1.n_links = pyo.Param(model_phase1.chain_idx,initialize=n_links_init)


def chain_init(M,t,k,i,j,w):
    window_list =[(i,j,w)]
    # Compute global period of (i,j,w)
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




##### Shift pyo.Variables

#model_phase1.Shift_cost = pyo.Param(model_phase1.ok_shifts,  initialize=Shift_cost_init)



model_phase1.Shift = pyo.Var(model_phase1.okShifts, within=pyo.NonNegativeIntegers)

# Shift[i,j,w,k,t] = Number of shifts of length k starting in period i
# of day j in week w for a tour of type t

##### Windowed tour type pyo.Variables
# TourType[i,j] Number of employees working tour type j starting in window i  


#/* set okTourType := setof{i in WINDOWS,t in okTTYPES :
#    sum{j in DAYS} (if (i,j) in ok_window_beginnings[t] then 1 else 0) >= tt_min_dys_week[t]} (i,t);*/
#    
#set okTourType := setof{i in WINDOWS,t in okTTYPES} (i,t);

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
    for (i,t) in M.okTourType:
        for p in pyo.sequence(M.max_weekend_patterns):
            weekendtype = M.weekend_type[i,t]
            if p <= M.num_weekend_patterns[weekendtype,t]:
                index_list.append((i, t, p))
    
    return index_list


model_phase1.weekenddaysworked_idx = pyo.Set(dimen=3, initialize=weekenddaysworked_idx_rule)

model_phase1.WeekendDaysWorked = pyo.Var(model_phase1.weekenddaysworked_idx, within=pyo.NonNegativeIntegers)
    
    # WeekendDaysWorked[d,i,t] = Number of employees working days-off patterns d
    # in start window i and of tour type t


def multiweekdaysworked_idx_rule(M):
    index_list = []
    for (i,t) in M.okTourType:
        for p in pyo.sequence(M.max_mwdw_patterns):
            if p <= M.num_mwdw_patterns[t]:
                index_list.append((i, t, p))

    return index_list


model_phase1.multiweekdaysworked_idx = pyo.Set(dimen=3, initialize=multiweekdaysworked_idx_rule)

model_phase1.MultiWeekDaysWorked = pyo.Var(model_phase1.multiweekdaysworked_idx, within=pyo.NonNegativeIntegers)

# MultiWeekDaysWorked[t, p] = Number of employees of tour type t working mwdw pattern p

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
    obj3 = sum(M.WeekendDaysWorked[i, t, p] * M.A_tot_wkend_days[p, t, e] ** 2 for (i, t, p) in M.weekenddaysworked_idx for e in M.WEEKENDS if M.weekend_type[i, t] == e)
    return obj1 + obj2 + obj3


model_phase1.total_cost = pyo.Objective(rule=objective_rule, sense=pyo.minimize)


##### Budget constraints

def max_labor_budget_rule(M): 
    return sum(M.Shift[i,j,w,k,t] * M.lengths[k] * M.tt_cost_multiplier[t] for (i, j, w, k, t) in M.okShifts) <= M.labor_budget

model_phase1.max_labor_budget = pyo.Constraint(rule=max_labor_budget_rule)

##### Coverage constraints

# Breaking them up into four different constraints, one for each case in terms of handling end of day horizon wrapping
# Breaking them up into four different constraints, one for each case in terms of handling end of day horizon wrapping

# subject to coverage1{i in PERIODS,j in DAYS,w in WEEKS} :
#  cov1[i,j,w] = sum{l in LENGTHS,t in okTTYPES,p in 0..lengths[l]-1: (i-p)>0} if allow_start[(i-p),j,l,t]>0 then
#                      Shift[(i-p),j,w,l,t] else 0;

def coverage1_rule(M,i,j,w):
    return sum(M.Shift[(i-p),j,w,l,t] for l in M.LENGTHS
                             for t in M.activeTT
                             for p in range(0,M.lengths[l])
                             if (i-p)>0 and M.allow_start[(i-p),j,l,t]>0) == M.cov1[i,j,w]

model_phase1.coverage1 = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, rule=coverage1_rule)

# subject to coverage2{i in PERIODS,j in DAYS,w in WEEKS} :
#  cov2[i,j,w] = sum{l in LENGTHS,t in okTTYPES,p in 0..lengths[l]-1: (i-p)<=0 and j>1}
#       if allow_start[n_prds_per_day+(i-p),j-1,l,t]>0 then
#         Shift[n_prds_per_day+(i-p),j-1,w,l,t] else 0;

def coverage2_rule(M,i,j,w):
    return sum(M.Shift[M.n_prds_per_day.value+(i-p),j-1,w,l,t] for l in M.LENGTHS
                             for t in M.activeTT
                             for p in range(0,M.lengths[l])
                             if (i-p)<=0 and j>1 and
               M.allow_start[M.n_prds_per_day.value+(i-p),j-1,l,t]>0) == M.cov2[i,j,w]

model_phase1.coverage2 = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, rule=coverage2_rule)

# subject to coverage3{i in PERIODS,j in DAYS,w in WEEKS} :
#  cov3[i,j,w] = sum{l in LENGTHS,t in okTTYPES,p in 0..lengths[l]-1: (i-p)<=0 and j=1 and w>1}
#       if allow_start[n_prds_per_day+(i-p),n_days_per_week,l,t]>0 then
#     Shift[n_prds_per_day+(i-p),n_days_per_week,w-1,l,t] else 0;

def coverage3_rule(M,i,j,w):
    return sum(M.Shift[M.n_prds_per_day.value+(i-p),M.n_days_per_week.value,w-1,l,t] for l in M.LENGTHS
                             for t in M.activeTT
                             for p in range(0,M.lengths[l])
                             if (i-p)<=0 and j==1 and w>1 and
               M.allow_start[M.n_prds_per_day.value+(i-p),M.n_days_per_week.value,l,t]>0) == M.cov3[i,j,w]

model_phase1.coverage3 = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, rule=coverage3_rule)

# subject to coverage4{i in PERIODS,j in DAYS,w in WEEKS} :
#  cov4[i,j,w] = sum{l in LENGTHS,t in okTTYPES,p in 0..lengths[l]-1: (i-p)<=0 and j=1 and w=1}
#      if allow_start[n_prds_per_day+(i-p),n_days_per_week,l,t]>0 then
#    Shift[n_prds_per_day+(i-p),n_days_per_week,n_weeks,l,t] else 0;

def coverage4_rule(M,i,j,w):
    return sum(M.Shift[M.n_prds_per_day.value+(i-p),M.n_days_per_week.value,M.n_weeks.value,l,t] for l in M.LENGTHS
                             for t in M.activeTT
                             for p in range(0,M.lengths[l])
                             if (i-p)<=0 and j==1 and w==1 and
               M.allow_start[M.n_prds_per_day.value+(i-p),M.n_days_per_week.value,l,t]>0) == M.cov4[i,j,w]

model_phase1.coverage4 = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, rule=coverage4_rule)


def tot_coverage_rule(M,i,j,w):
    return M.cov1[i,j,w] + M.cov2[i,j,w] + M.cov3[i,j,w] + M.cov4[i,j,w] == M.cov[i,j,w]


model_phase1.tot_coverage = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, rule=tot_coverage_rule)


def coverage_rule(M,i,j,w):
    return M.cov[i,j,w] + M.under1[i,j,w] + M.under2[i,j,w] >= M.dmd_staff[i,j,w]


model_phase1.coverage = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, rule=coverage_rule)


def minstaff_rule(M,i,j,w):
    return M.cov[i,j,w] >= M.min_staff[i,j,w]


model_phase1.minstaff = pyo.Constraint(model_phase1.PERIODS, model_phase1.DAYS, model_phase1.WEEKS, rule=minstaff_rule)

#### WT_bounds - bounds from .MIX file
#
def TourType_LB_rule(M,t):
    return sum(M.TourType[i,t] for (i,s) in M.okTourType if s == t) >= M.tt_lb[t]

model_phase1.TourType_LB_con = pyo.Constraint(model_phase1.activeTT, rule=TourType_LB_rule)  

def TourType_UB_rule(M,t):
    return sum(M.TourType[i,t] for (i,s) in M.okTourType if s == t) <= M.tt_ub[t]

model_phase1.TourType_UB_con = pyo.Constraint(model_phase1.activeTT, rule=TourType_UB_rule)  

# Each tour pyo.Variable must get a WeekendDaysWorked pattern assigned to it.
#  subject to weekend_total{(i,t) in okTourType} :
#    (sum{p in 1..num_weekend_patterns[weekend_type[i,t],t]} WeekendDaysWorked[p,i,t]) - TourType[i,t]  = 0;


def weekend_total_rule(M,i,t):
    """
    Each tour pyo.Variable must get a WeekendDaysWorked pattern assigned to it
    """
    return sum(M.WeekendDaysWorked[i,t,p] for p in pyo.sequence(M.num_weekend_patterns[M.weekend_type[i,t],t])) \
                == M.TourType[i,t]
                           

def weekend_total_idx_rule(M):
    """
    The index for the weekend_total constraints is (window, ttype) tuples
    in the set okTourType.
    """
    return [(i,t) for (i,t) in M.okTourType]
                        
model_phase1.weekend_total_idx = pyo.Set(dimen=2, initialize=weekend_total_idx_rule)

model_phase1.weekend_total_con = pyo.Constraint(model_phase1.weekend_total_idx, rule=weekend_total_rule)


# Integrate days-off and shift scheduling sub-problems

"""
The number of people working on Sun or Sat as part of a type 1 weekend worked pattern must be equal to the
number of people working on Sun or Sat as specified by DailyTourType pyo.Variables.
"""

# subject to weekend_integration_1{j in DAYS,w in WEEKS,i in 1..n_windows,t in okTTYPES :
#    j in weekend[i,t] and 1 in weekend[i,t] and 7 in weekend[i,t] and (i,t,j) in okDailyTourType} :
#      sum{p in 1..num_weekend_patterns[weekend_type[i,t],t]} A[p,j,w,t,1]*WeekendDaysWorked[p,i,t] =
#        DailyTourType[i,t,j,w];
      

def weekend_integration_1_SS_idx_rule(M):
    return [(j,w,i,t) for j in M.DAYS
                    for w in M.WEEKS
                    for i in M.WINDOWS                        
                    for t in M.activeTT
                    if j in M.weekend[i,t] and 1 in M.weekend[i,t] and 7 in M.weekend[i,t] \
                      and (i,t,j) in M.okDailyTourType]
                        
model_phase1.weekend_integration_1_SS_idx = pyo.Set(dimen=4, initialize=weekend_integration_1_SS_idx_rule)


def weekend_integration_1_SS_rule(M,j,w,i,t):
    """
    The number of people working on Sun or Sat as part of a type 1 weekend worked pattern must be equal to the
    number of people working on Sun or Sat as specified by DailyTourType pyo.Variables.
    """
    return sum(M.A[p,j,w,t,1] * M.WeekendDaysWorked[i,t,p] for p in pyo.sequence(M.num_weekend_patterns[1,t])) \
                == M.DailyTourType[i,t,j,w]
                

model_phase1.weekend_integration_1_SS_con = pyo.Constraint(model_phase1.weekend_integration_1_SS_idx,
                                                       rule=weekend_integration_1_SS_rule)


# subject to weekend_integration_2{j in DAYS,w in WEEKS,i in 1..n_windows,
#    t in okTTYPES : j in weekend[i,t] and 6 in weekend[i,t] and 7 in weekend[i,t] and (i,t,j) in okDailyTourType} :
#     sum{p in 1..num_weekend_patterns[weekend_type[i,t],t]} A[p,j,w,t,2]*WeekendDaysWorked[p,i,t] =
#         DailyTourType[i,t,j,w];


def weekend_integration_2_FS_rule(M,j,w,i,t):
    """
    The number of people working on Fri or Sat as part of a type 2 weekend worked pattern must be equal to the
    number of people working on Sun or Sat as specified by DailyTourType pyo.Variables.
    """
    return sum(M.A[p,j,w,t,2] * M.WeekendDaysWorked[i,t,p] for p in pyo.sequence(M.num_weekend_patterns[2,t])) \
                == M.DailyTourType[i,t,j,w]


def weekend_integration_2_FS_idx_rule(M):
    return [(j,w,i,t) for j in M.DAYS
                    for w in M.WEEKS
                    for i in M.WINDOWS                        
                    for t in M.activeTT
                    if j in M.weekend[i,t] and 6 in M.weekend[i,t] and 7 in M.weekend[i,t] \
                      and (i,t,j) in M.okDailyTourType]
                        
model_phase1.weekend_integration_2_FS_idx = pyo.Set(dimen=4,initialize=weekend_integration_2_FS_idx_rule)  

model_phase1.weekend_integration_2_FS_con = pyo.Constraint(model_phase1.weekend_integration_2_FS_idx,rule=weekend_integration_2_FS_rule)

# Integrate shift, days worked, and tour type pyo.Variables

# Make sure number of shifts scheduled each day for each tour type doesn't exceed the number of
# tour types scheduled to work that day

# Make sure number of shifts scheduled each week for each tour is equal to the number of days worked each week

# TODO: I think the following is redundant given DTT_DSW constraints
#
# subject to shift_DTT_dailyconservation{i in 1..n_windows, j in DAYS, w in WEEKS, t in TTYPES : (i,t,j) in okDailyTourType } :
#    sum{k in LENGTHS: k in tt_length_x[t] and (i,j,k,t) in ok_shifts} Shift[i,j,w,k,t] = DailyTourType[i,t,j,w];


# TODO The following are only correct for start window = 0
# def shift_DTT_dailyconservation_rule(M,i,j,w,t):
#    return sum(M.Shift[i,j,w,k,t] for k in M.tt_length_x[t] if (i,j,w,k,t) in M.okShifts) \
#                                  == M.DailyTourType[i,t,j,w]


# The index set gets used in the DTT_DSW_con and DTT_TT_UB

def shift_DTT_dailyconservation_idx_rule(M):    
    return [(i,j,w,t) for i in M.WINDOWS 
                      for j in M.DAYS   
                      for w in M.WEEKS                    
                      for t in M.activeTT
                      if (i,t,j) in M.okDailyTourType]
                         

model_phase1.shift_DTT_dailyconservation_idx = pyo.Set(dimen=4, initialize=shift_DTT_dailyconservation_idx_rule)  

# model_phase1.shift_DTT_dailyconservation_con =
#    Constraint(model_phase1.shift_DTT_dailyconservation_idx,rule=shift_DTT_dailyconservation_rule)


##
#### ---- Make sure number of people working on any given day <= number of employees
###subject to DWT_WT_UB{j in DAYS,i in 1..n_windows, t in okTTYPES, w in WEEKS : (i,t,j) in okDailyTourType} :
###    DailyTourType[i,t,j,w] <= TourType[i,t];
###    
###    
###subject to DWT_DSW{j in DAYS,i in 1..n_windows, t in okTTYPES, w in WEEKS : (i,t,j) in okDailyTourType} :
###    sum{k in tt_length_x[t]}DailyShiftWorked[i,t,k,j,w] = DailyTourType[i,t,j,w];
#### ------------------------------------------------------------------------------------------------    
##### ---- Make sure number of people working on any given day <= number of employees


def DTT_TT_UB_rule(M,i,j,w,t):
    """
    Every day of every week for each (window,ttype), there can be no more people scheduled (DailyTourType) than
    number of people assigned to TourType[window,ttype]
    """
    return M.DailyTourType[i,t,j,w] <= M.TourType[i,t]
    
 
model_phase1.DTT_TT_UB = \
    pyo.Constraint(model_phase1.shift_DTT_dailyconservation_idx, rule=DTT_TT_UB_rule)


# For case where two weekend days worked and min days worked per week = 2, we do a heuristic
# adjustment to the max of DailyTourType each day to avoid infeasibility due to forcing
# a two weekend day person to work a third day and perhaps leading to some other tour not
# getting >= 2 shifts over the week


def dailyconservation_wkendadj_idx_rule(M):
    return [(i,j,w,t) for i in M.WINDOWS 
                      for j in M.DAYS   
                      for w in M.WEEKS                    
                      for t in M.activeTT
                      if (i,t,j) in M.okDailyTourType and M.tt_min_dys_weeks[t,w] == 2]
                         

model_phase1.dailyconservation_wkendadj_idx = pyo.Set(dimen=4, initialize=dailyconservation_wkendadj_idx_rule)


def DTT_TT_fullwkendadj_UB_rule(M,i,j,w,t):
    return M.DailyTourType[i,t,j,w] <= M.TourType[i,t] - sum(M.WeekendDaysWorked[i,t,p]
                                                             for p in M.two_wkend_days[w,t,M.weekend_type[i,t]])
    
model_phase1.DTT_TT_fullwkendadj_UB = pyo.Constraint(model_phase1.dailyconservation_wkendadj_idx,
                                                 rule=DTT_TT_fullwkendadj_UB_rule)

# def DTT_TT_halfwkendadj_UB_rule(M,i,j,w,t):
#    return M.DailyTourType[i,t,j,w] <= M.TourType[i,t] - sum(M.WeekendDaysWorked[i,t,p] for p in M.two_wkend_days[w,t,M.weekend_type[i,t]])
#    
# model_phase1.DTT_TT_fullwkendadj_UB = Constraint(model_phase1.dailyconservation_wkendadj_idx,rule=DTT_TT_fullwkendadj_UB_rule)


# Coordinate DailyShiftWorked and DailyTourType pyo.Variables 
# For each day of week in each week


def DTT_DSW_rule(M,i,j,w,t): 
    return sum(M.DailyShiftWorked[i,t,k,j,w] for k in M.tt_length_x[t]) == M.DailyTourType[i,t,j,w]  


model_phase1.DTT_DSW_con = pyo.Constraint(model_phase1.shift_DTT_dailyconservation_idx,rule=DTT_DSW_rule)

# ---- Min and max bounds on days worked each week

# subject to DWT_WT_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, w in WEEKS } :
#    sum{d in DAYS: (i,t,d) in okDailyTourType} DailyTourType[i,t,d,w] >= TourType[i,t]*tt_min_dys_weeks[t,w];

# subject to DWT_WT_weeklyconservation_UB{i in 1..n_windows, t in okTTYPES, w in WEEKS } :
#    sum{d in DAYS: (i,t,d) in okDailyTourType} DailyTourType[i,t,d,w] <= TourType[i,t]*tt_max_dys_weeks[t,w];


# Index for both lower and upper bound versions of these constraints are the same

def DTT_TT_weeklyconservation_idx_rule(M):
    return [(i,t,w) for (i,t) in M.okTourType
                        for w in M.WEEKS]


model_phase1.DTT_TT_weeklyconservation_idx = pyo.Set(dimen=3, initialize=DTT_TT_weeklyconservation_idx_rule)


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

def DTT_TT_cumul_weeklyconservation_LB_rule(M,i,t,w):
    return sum(M.DailyTourType[i,t,d,z] for d in M.DAYS for z in pyo.sequence(w)) >= \
           M.TourType[i,t]*M.tt_min_cumul_dys_weeks[t,w]


model_phase1.DTT_TT_cumul_weeklyconservation_LB = \
  pyo.Constraint(model_phase1.DTT_TT_weeklyconservation_idx, rule=DTT_TT_cumul_weeklyconservation_LB_rule)


def DTT_TT_cumul_weeklyconservation_UB_rule(M, i, t, w):
    return sum(M.DailyTourType[i, t, d, z] for d in M.DAYS for z in pyo.sequence(w)) <= \
           M.TourType[i, t] * M.tt_max_cumul_dys_weeks[t, w]


model_phase1.DTT_TT_cumul_weeklyconservation_UB = \
    pyo.Constraint(model_phase1.DTT_TT_weeklyconservation_idx, rule=DTT_TT_cumul_weeklyconservation_UB_rule)



def DSW_TT_weeklyconservation_idx_rule(M):
    return [(i,t,w) for (i,t) in M.okTourType for w in M.WEEKS]


model_phase1.DSW_TT_weeklyconservation_idx = pyo.Set(dimen=3,initialize=DSW_TT_weeklyconservation_idx_rule) 

                        
def DSW_TT_weeklyconservation_LB_rule(M,i,t,w):
    return sum(M.DailyShiftWorked[i,t,k,d,w] for d in M.DAYS for k in M.tt_length_x[t] if \
               (i,t,k,d,w) in M.DailyShiftWorked_idx) >= M.TourType[i,t]*M.tt_min_dys_weeks[t,w]


model_phase1.DSW_TT_weeklyconservation_LB = \
  pyo.Constraint(model_phase1.DSW_TT_weeklyconservation_idx, rule=DSW_TT_weeklyconservation_LB_rule)


def DSW_TT_weeklyconservation_UB_rule(M,i,t,w):
    return sum(M.DailyShiftWorked[i,t,k,d,w] for d in M.DAYS for k in M.tt_length_x[t] if \
               (i,t,k,d,w) in M.DailyShiftWorked_idx) <= M.TourType[i,t]*M.tt_max_dys_weeks[t,w]


model_phase1.DSW_TT_weeklyconservation_UB = \
  pyo.Constraint(model_phase1.DSW_TT_weeklyconservation_idx, rule=DSW_TT_weeklyconservation_UB_rule)


def DSW_TT_cumul_weeklyconservation_LB_rule(M,i,t,w):
    return sum(M.DailyShiftWorked[i,t,k,d,z] for d in M.DAYS for k in M.tt_length_x[t] for \
               z in pyo.sequence(w) if (i,t,k,d,z) in M.DailyShiftWorked_idx) >= \
           M.TourType[i,t]*M.tt_min_cumul_dys_weeks[t,w]

model_phase1.DSW_TT_cumul_weeklyconservation_LB = \
  pyo.Constraint(model_phase1.DSW_TT_weeklyconservation_idx, rule=DSW_TT_cumul_weeklyconservation_LB_rule)


def DSW_TT_cumul_weeklyconservation_UB_rule(M,i,t,w):
    return sum(M.DailyShiftWorked[i,t,k,d,z] for d in M.DAYS for k in M.tt_length_x[t] for \
               z in pyo.sequence(w) if (i,t,k,d,z) in M.DailyShiftWorked_idx) <= \
           M.TourType[i,t]*M.tt_max_cumul_dys_weeks[t,w]

model_phase1.DSW_TT_cumul_weeklyconservation_UB = \
  pyo.Constraint(model_phase1.DSW_TT_weeklyconservation_idx, rule=DSW_TT_cumul_weeklyconservation_UB_rule)


# Shiftlen versions of the above.

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
#  
#  
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


##### ---- OK - Min and max bounds on prds worked each week
##Shift[i,j,w,k,t]*lengths[k] >= TourType[i,t]*tt_min_prds_weeks[t,w]
#def prds_worked_weekly_LB_rule(M,i,t,w):
#    return sum(M.Shift[i,d,w,k,t]*M.lengths[k] for d in M.DAYS for k in M.tt_length_x[t] if (i,d,w,k,t) in M.okShifts) \
#        >= M.TourType[i,t]*M.tt_min_prds_weeks[t,w]


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


## Cumulative versions of the above 2 constraints
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


# The following were commented out for some reason though they definitely appear in GMPL model

# Shiftlen version of the above (reused model_phase1.prds_worked_shiflen_weekly_idx from above)

# ---- OK - Min and max bounds on shifts worked each week
#subject to shift_WT_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, k in LENGTHS, w in WEEKS } :
#    sum{d in DAYS: (i,d,k,t) in ok_shifts} Shift[i,d,w,k,t] >= TourType[i,t]*tt_shiftlen_min_dys_weeks[t,k,w];
#
#subject to shift_WT_weeklyconservation_UB{i in 1..n_windows, t in okTTYPES, k in LENGTHS, w in WEEKS } :
#    sum{d in DAYS: (i,d,k,t) in ok_shifts} Shift[i,d,w,k,t] <= TourType[i,t]*tt_shiftlen_max_dys_weeks[t,k,w];
#
#
## ---- OK - Min and max bounds on shifts worked over cumulative number of weeks
# subject to shift_WT_cumul_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, k in LENGTHS, w in WEEKS } :
#    sum{d in DAYS,z in 1..w: (i,d,k,t) in ok_shifts} Shift[i,d,z,k,t] >= TourType[i,t]*tt_shiftlen_min_cumul_dys_weeks[t,k,w];
#
# subject to shift_WT_cumul_weeklyconservation_UB{i in 1..n_windows, t in okTTYPES, k in LENGTHS, w in WEEKS } :
#    sum{d in DAYS,z in 1..w: (i,d,k,t) in ok_shifts} Shift[i,d,z,k,t] <= TourType[i,t]*tt_shiftlen_max_cumul_dys_weeks[t,k,w];

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


#
#####subject to DWT_ID_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, w in WEEKS:
##### 6 in weekend[i,t] and 7 in weekend[i,t]} :
#####    sum{d in DAYS: (i,t,d) in okDailyTourType} DailyTourType[i,t,d,w] >= 
#####    sum{p in 1..num_weekend_patterns[t]}WeekendDaysWorked[p,i,t]*FriSat_min_dys_weeks[t,p,w];
#####    
#####subject to DWT_ID_cumul_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, w in WEEKS:
##### 6 in weekend[i,t] and 7 in weekend[i,t]} :
#####    sum{d in DAYS, z in 1..w: (i,t,d) in okDailyTourType} DailyTourType[i,t,d,z] >= 
#####    sum{p in 1..num_weekend_patterns[t]}WeekendDaysWorked[p,i,t]*FriSat_min_cumul_dys_weeks[t,p,w];
##
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
                        

model_phase1.DWT_FSwkend_weeklyconservation_LB_idx = pyo.Set(dimen=3,initialize=DWT_FSwkend_weeklyconservation_LB_idx_rule)
model_phase1.DWT_FSwkend_weeklyconservation_LB = pyo.Constraint(model_phase1.DWT_FSwkend_weeklyconservation_LB_idx,rule=DWT_FSwkend_weeklyconservation_LB_rule)
model_phase1.DWT_FSwkend_cumul_weeklyconservation_LB = pyo.Constraint(model_phase1.DWT_FSwkend_weeklyconservation_LB_idx,rule=DWT_FSwkend_cumul_weeklyconservation_LB_rule)

# TODO Reformulate the following max_ptfrac constraint as I think it's causing problems
##### Part-time fraction constraints
####
####
#####subject to max_ptfrac :
#####   (  sum{i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS,t in TTYPES: ptfrac_tog > 0 and (i,j,k,t) in ok_shifts} 
#####       ( if tt_parttime[t] > 0 then 
#####        Shift[i,j,w,k,t]*lengths[k]*(1.0-1.0/ptfrac) else
#####            Shift[i,j,w,k,t]*lengths[k] ) ) >= 0.0 ;
##

def max_ptfrac_rule(M):
    return sum(M.Shift[i,j,w,k,t] * M.lengths[k] for (i,j,w,k,t) in M.okShifts if M.tt_parttime[t] > 0) <= \
        M.max_parttime_frac.value * sum(M.Shift[i,j,w,k,t] * M.lengths[k] for (i,j,w,k,t) in M.okShifts)
    
    
model_phase1.max_ptfrac_con = pyo.Constraint(rule=max_ptfrac_rule)


# Chains - coordinates DWT and ix within each chain

# How to deal with multiple weeks?
# How to deal with fact that single tour type can have multiple shift lengths? It's a problem since we CANNOT equate sums
# of shift pyo.Variables for a specific shift length with DailyTourType pyo.Variables.

#subject to chains_sweep_l{e in WEEKS, t in okTTYPES, k in tt_length_x[t], (b,j) in bchain[t,k],
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
           sum(M.DailyShiftWorked[g_prd_to_tuple(M, u)[0], t, k, g_prd_to_tuple(M, u)[1], g_prd_to_tuple(M, u)[2]] \
               for u in [vv for vv in range(p, p + M.g_start_window_width + 1)
                         if (g_prd_to_tuple(M, v)[0], g_prd_to_tuple(M, v)[1], g_prd_to_tuple(M, v)[2]) \
                         in M.okStartWindowRoots[t, k] and sum(M.allow_start[x, y, k, t] for (x, y, z) \
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
                                           
    
#    return [(t,k,b,j,w,m) for t in M.activeTT
#                          for k in M.LENGTHS
#                          for b in M.PERIODS
#                          for j in M.DAYS
#                          for w in M.WEEKS
#                          for p in range(M.g_period[b,j,w].value,M.g_period[b,j,w].value+1)
#                          for m in range(0,M.n_links[t,k,b,j,w].value-1-(p-M.g_period[b,j,w].value)+1)
#                          if M.g_start_window_width>0 and (b,j,w) in M.bchain[t,k] and k in M.tt_length_x[t]]


model_phase1.chains_sweep_l_idx = pyo.Set(dimen=7, ordered=True, initialize=chains_sweep_l_idx_rule)
model_phase1.chains_sweep_l_con = pyo.Constraint(model_phase1.chains_sweep_l_idx, rule=chains_sweep_l_rule)


#subject to chains_sweep_u{e in WEEKS,t in okTTYPES, k in tt_length_x[t],(b,j) in bchain[t,k],
#             w in 0..(numlinks[t,k,b,j]-1):
#   width>0 } :
#    sum{ i in period[b,j]..period[b,j]+w :
#     (which_prd[i],which_day[i],k,t) in ok_shifts} 
#      Shift[which_prd[i],which_day[i],e,k,t] <=
#     sum{u in period[b,j]..period[b,j]+w : (which_prd[u],which_day[u]) in okWindowBeginnings[t,k] and 
#      sum{(l,m) in WindowWepochs[which_prd[u],which_day[u]]} allow_start[l,m,k,t]>0} 
#      DailyShiftWorked[which_prd[u],t,k,which_day[u],e] ; 
      
def chains_sweep_u_rule(M,t,k,b,j,w,p,v):
    
    return  sum(M.Shift[g_prd_to_tuple(M, i)[0],g_prd_to_tuple(M, i)[1],g_prd_to_tuple(M, i)[2],k,t] 
              for i in range(p,p+v+1) if (g_prd_to_tuple(M, i)[0],g_prd_to_tuple(M, i)[1],g_prd_to_tuple(M, i)[2]) in M.okShifts) <= sum(M.DailyShiftWorked[g_prd_to_tuple(M, u)[0],t,k,g_prd_to_tuple(M, u)[1],g_prd_to_tuple(M, u)[2]] 
                   for u in [vv for vv in range(p,p + M.g_start_window_width + 1)
                if (g_prd_to_tuple(M, v)[0],g_prd_to_tuple(M, v)[1],g_prd_to_tuple(M, v)[2]) in M.okStartWindowRoots[t,k] \
                  and sum(M.allow_start[x,y,k,t] for (x,y,z) in M.PotentialGlobalStartWindow[g_prd_to_tuple(M, v)[0],g_prd_to_tuple(M, v)[1],g_prd_to_tuple(M, v)[2]])>0])
                           
def chains_sweep_u_idx_rule(M):
    index_list = []
    if M.g_start_window_width>0:
        for t in M.activeTT:
            for k in [len for len in M.LENGTHS if len in M.tt_length_x[t]]:
                for b in M.PERIODS:
                    for j in M.DAYS:
                        for w in M.WEEKS:
                            for p in range(M.g_period[b,j,w].value,M.g_period[b,j,w]+1):
                                if (b,j,w) in M.bchain[t,k]:
                                    for m in range(0,M.n_links[t,k,b,j,w]-1):
                                        index_list.append((t,k,b,j,w,p,m))
                                        
    return index_list

                        
model_phase1.chains_sweep_u_idx = pyo.Set(dimen=7,ordered=True,initialize=chains_sweep_u_idx_rule)  
model_phase1.chains_sweep_u_con = pyo.Constraint(model_phase1.chains_sweep_l_idx,rule=chains_sweep_u_rule)     
      
      


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





# TODO - the follow proxy constraints are only for PotentialGlobalStartWindow = 0
# Coordinate Shift pyo.Variables with DailyShiftWorked for start window width=0.
# In this case, shift periods are same as shiftworked windows

def chains_tot_proxy1_rule(M,w,t,k,i,j):
    return M.Shift[i,j,w,k,t] == M.DailyShiftWorked[i,t,k,j,w]


def chains_tot_proxy1_idx_rule(M):
    return [(w,t,k,i,j) for w in M.WEEKS
                        for t in M.activeTT
                        for k in M.tt_length_x[t]
                        for i in M.PERIODS
                        for j in M.DAYS if (i,j,w,k,t) in M.okShifts ]                      
                    

model_phase1.chains_tot_proxy1_idx = pyo.Set(dimen=5, initialize=chains_tot_proxy1_idx_rule)

model_phase1.chains_tot_proxy1_con = pyo.Constraint(model_phase1.chains_tot_proxy1_idx, rule=chains_tot_proxy1_rule)


def chains_tot_proxy2_rule(M,w,t,k,i,j):
    return M.DailyShiftWorked[i,t,k,j,w] == 0


def chains_tot_proxy2_idx_rule(M):
    return [(w,t,k,i,j) for w in M.WEEKS
                        for t in M.activeTT
                        for k in M.tt_length_x[t]
                        for i in M.PERIODS
                        for j in M.DAYS
                        if (i,t,j) in M.okDailyTourType and M.allow_start[i,j,k,t] == 0] 
                    

model_phase1.chains_tot_proxy2_idx = pyo.Set(dimen=5,initialize=chains_tot_proxy2_idx_rule)  

model_phase1.chains_tot_proxy2_con = pyo.Constraint(model_phase1.chains_tot_proxy2_idx,rule=chains_tot_proxy2_rule)

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
    
    # total days in subset worked by all wkend patterns -  days worked by those with < 2 weekend days <= 4x-x or 3x
    # where x is number of weekend patterns with 2 wkend days

    days = [d1, d2, d3, d4]
    return sum(M.DailyTourType[i, t, d, w] for d in days) <= \
        (len(days)) * sum(M.WeekendDaysWorked[i, t, p] for p in M.oneorzero_wkend_day[w, t, e]) \
        + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])


model_phase1.weekend_subsets_5_4_con = pyo.Constraint(model_phase1.weekend_subsets_5_4_idx,
                                                      rule=weekend_subsets_5_4_rule)
    

# weekend_subsets_5_5 - put UB on number of M-F shifts based on number of weekend patterns used
# with 0, 1, 2 weekend days worked.
# On second thought, this seems totally redundant in that the total number of shifts per week
# bounds in conjunction with the
# integration of weekends off pyo.Variables and daily tour type pyo.Variables should totally
# determine the number of M-F shifts.

def weekend_subsets_5_5_idx_rule(M):
    index_list = []

    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] == 5:
                    index_list.append((i, t, w, e, 2, 3, 4, 5, 6))

    return index_list


model_phase1.weekend_subsets_5_5_idx = pyo.Set(dimen=9, initialize=weekend_subsets_5_5_idx_rule)


def weekend_subsets_5_5_rule(M, i, t, w, e, d1, d2, d3, d4, d5):
    days = [d1, d2, d3, d4, d5]
    return sum(M.DailyTourType[i, t, d, w] for d in days) <= (len(days)) * sum(
        M.WeekendDaysWorked[i, t, p] for p in M.zero_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.one_wkend_day[w, t, e]) \
           + (len(days) - 2) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])


model_phase1.weekend_subsets_5_5_con = pyo.Constraint(model_phase1.weekend_subsets_5_5_idx,
                                                      rule=weekend_subsets_5_5_rule)


def weekend_subsets_5_5lb_idx_rule(M):
    index_list = []

    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] == 5:
                    index_list.append((i, t, w, e, 1, 2, 3, 4, 5, 6))
                    index_list.append((i, t, w, e, 2, 3, 4, 5, 6, 7))

    return index_list


model_phase1.weekend_subsets_5_5lb_idx = pyo.Set(dimen=10, initialize=weekend_subsets_5_5lb_idx_rule)


def weekend_subsets_5_5lb_rule(M, i, t, w, e, d1, d2, d3, d4, d5, d6):
    days = [d1, d2, d3, d4, d5, d6]
    if 1 in days:
        return sum(M.DailyTourType[i, t, d, w] for d in days) >= (len(days) - 1) * sum(
            M.WeekendDaysWorked[i, t, p] for p in M.Sun_wkend_day[w, t, e])
    else:
        return sum(M.DailyTourType[i, t, d, w] for d in days) >= (len(days) - 1) * sum(
            M.WeekendDaysWorked[i, t, p] for p in M.Sat_wkend_day[w, t, e])


model_phase1.weekend_subsets_5_5lb_con = pyo.Constraint(model_phase1.weekend_subsets_5_5lb_idx,
                                                        rule=weekend_subsets_5_5lb_rule)


#def weekend_subsets_5_4lb_idx_rule(M):
#    index_list = []
#    
#    for (i,t) in M.okTourType:
#        for w in M.WEEKS:
#            for e in M.WEEKENDS:
#                if M.tt_max_dys_weeks[t,w].value ==  5:
#                    index_list.append((i,t,w,e,2,3,4,5,6))
#                    
#                    
#        
#    return index_list
#        
#model_phase1.weekend_subsets_5_4lb_idx = Set(dimen=9,initialize=weekend_subsets_5_4lb_idx_rule)    
#    
#def weekend_subsets_5_4lb_rule(M,i,t,w,e,d1,d2,d3,d4,d5):
#    
#    # total days in subset worked by all wkend patterns -  days worked by those with < 2 weekend days <= 4x-x or 3x
#    # where x is number of weekend patterns with 2 wkend days
#    # M.WeekendDaysWorked[i,t,p]
#    days = [d1,d2,d3,d4,d5]
#    return sum(M.DailyTourType[i,t,d,w] for d in days) >= (len(days))*sum(M.WeekendDaysWorked[i,t,p] for p in M.zero_wkend_day[w,t,e]) \
#                                                          + (len(days)-1)*sum(M.WeekendDaysWorked[i,t,p] for p in M.one_wkend_day[w,t,e]) \
#                                                          + (len(days)-2)*sum(M.WeekendDaysWorked[i,t,p] for p in M.two_wkend_days[w,t,e])
#    
#model_phase1.weekend_subsets_5_4lb_con = Constraint(model_phase1.weekend_subsets_5_4lb_idx,rule=weekend_subsets_5_4lb_rule) 


def weekend_subsets_5_5sun_idx_rule(M):
    index_list = []

    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] == 5:
                    # don't think I need this one since in 5_5 - index_list.append((i,t,w,e,2,3,4,5,6))
                    index_list.append((i, t, w, e, 2, 3, 4, 5, 7))
                    index_list.append((i, t, w, e, 2, 3, 4, 6, 7))
                    index_list.append((i, t, w, e, 2, 3, 5, 6, 7))
                    index_list.append((i, t, w, e, 2, 4, 5, 6, 7))
                    index_list.append((i, t, w, e, 3, 4, 5, 6, 7))

    return index_list


model_phase1.weekend_subsets_5_5sun_idx = pyo.Set(dimen=9, initialize=weekend_subsets_5_5sun_idx_rule)


def weekend_subsets_5_5sun_rule(M, i, t, w, e, d1, d2, d3, d4, d5):
    # total days in subset worked by all wkend patterns -  days worked by those with < 2 weekend days <= 4x-x or 3x
    # where x is number of weekend patterns with 2 wkend days
    # M.WeekendDaysWorked[i,t,p]
    days = [d1, d2, d3, d4, d5]
    return sum(M.DailyTourType[i, t, d, w] for d in days) <= (len(days)) * sum(
        M.WeekendDaysWorked[i, t, p] for p in M.Sat_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.Sun_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.zero_wkend_day[w, t, e])


model_phase1.weekend_subsets_5_5sun_con = pyo.Constraint(model_phase1.weekend_subsets_5_5sun_idx,
                                                         rule=weekend_subsets_5_5sun_rule)


def weekend_subsets_5_5sat_idx_rule(M):
    index_list = []

    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] == 5:
                    # don't think I need this one since in 5_5 - index_list.append((i,t,w,e,2,3,4,5,6))
                    index_list.append((i, t, w, e, 2, 3, 4, 5, 1))
                    index_list.append((i, t, w, e, 2, 3, 4, 6, 1))
                    index_list.append((i, t, w, e, 2, 3, 5, 6, 1))
                    index_list.append((i, t, w, e, 2, 4, 5, 6, 1))
                    index_list.append((i, t, w, e, 3, 4, 5, 6, 1))

    return index_list


model_phase1.weekend_subsets_5_5sat_idx = pyo.Set(dimen=9, initialize=weekend_subsets_5_5sat_idx_rule)


def weekend_subsets_5_5sat_rule(M, i, t, w, e, d1, d2, d3, d4, d5):
    # total days in subset worked by all wkend patterns -  days worked by those with < 2 weekend days <= 4x-x or 3x
    # where x is number of weekend patterns with 2 wkend days
    # M.WeekendDaysWorked[i,t,p]
    days = [d1, d2, d3, d4, d5]
    return sum(M.DailyTourType[i, t, d, w] for d in days) <= (len(days)) * sum(
        M.WeekendDaysWorked[i, t, p] for p in M.Sun_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.Sat_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.zero_wkend_day[w, t, e])


model_phase1.weekend_subsets_5_5sat_con = pyo.Constraint(model_phase1.weekend_subsets_5_5sat_idx,
                                                         rule=weekend_subsets_5_5sat_rule)


def weekend_subsets_5_4sat_idx_rule(M):
    index_list = []

    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] == 5:
                    # don't think I need this one since in 5_5 - index_list.append((i,t,w,e,2,3,4,5,6))
                    index_list.append((i, t, w, e, 1, 2, 3, 4))
                    index_list.append((i, t, w, e, 1, 2, 3, 5))
                    index_list.append((i, t, w, e, 1, 2, 3, 6))

                    index_list.append((i, t, w, e, 1, 2, 4, 5))
                    index_list.append((i, t, w, e, 1, 2, 4, 6))

                    index_list.append((i, t, w, e, 1, 3, 4, 5))
                    index_list.append((i, t, w, e, 1, 3, 4, 6))

    return index_list


model_phase1.weekend_subsets_5_4sat_idx = pyo.Set(dimen=8, initialize=weekend_subsets_5_4sat_idx_rule)


def weekend_subsets_5_4sat_rule(M, i, t, w, e, d1, d2, d3, d4):
    # total days in subset worked by all wkend patterns -  days worked by those with < 2 weekend days <= 4x-x or 3x
    # where x is number of weekend patterns with 2 wkend days
    # M.WeekendDaysWorked[i,t,p]
    days = [d1, d2, d3, d4]
    return sum(M.DailyTourType[i, t, d, w] for d in days) <= (len(days)) * sum(
        M.WeekendDaysWorked[i, t, p] for p in M.Sun_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.Sat_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.zero_wkend_day[w, t, e])


model_phase1.weekend_subsets_5_4sat_con = pyo.Constraint(model_phase1.weekend_subsets_5_4sat_idx,
                                                         rule=weekend_subsets_5_4sat_rule)


def weekend_subsets_5_4sun_idx_rule(M):
    index_list = []

    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_min_dys_weeks[t, w] == 5:
                    # don't think I need this one since in 5_5 - index_list.append((i,t,w,e,2,3,4,5,6))
                    index_list.append((i, t, w, e, 7, 2, 3, 4))
                    index_list.append((i, t, w, e, 7, 2, 3, 5))
                    index_list.append((i, t, w, e, 7, 2, 3, 6))

                    index_list.append((i, t, w, e, 7, 2, 4, 5))
                    index_list.append((i, t, w, e, 7, 2, 4, 6))

                    index_list.append((i, t, w, e, 7, 3, 4, 5))
                    index_list.append((i, t, w, e, 7, 3, 4, 6))

    return index_list


model_phase1.weekend_subsets_5_4sun_idx = pyo.Set(dimen=8, initialize=weekend_subsets_5_4sun_idx_rule)


def weekend_subsets_5_4sun_rule(M, i, t, w, e, d1, d2, d3, d4):
    # total days in subset worked by all wkend patterns -  days worked by those with < 2 weekend days <= 4x-x or 3x
    # where x is number of weekend patterns with 2 wkend days
    # M.WeekendDaysWorked[i,t,p]
    days = [d1, d2, d3, d4]
    return sum(M.DailyTourType[i, t, d, w] for d in days) <= (len(days)) * sum(
        M.WeekendDaysWorked[i, t, p] for p in M.Sat_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.Sun_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.zero_wkend_day[w, t, e])


model_phase1.weekend_subsets_5_4sun_con = pyo.Constraint(model_phase1.weekend_subsets_5_4sun_idx,
                                                         rule=weekend_subsets_5_4sun_rule)


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
    # total days in subset worked by all wkend patterns -  days worked by those with < 2 weekend days <= 4x-x or 3x
    # where x is number of weekend patterns with 2 wkend days
    # M.WeekendDaysWorked[i,t,p]
    days = [d1, d2, d3]
    return sum(M.DailyTourType[i, t, d, w] for d in days) <= (len(days)) * sum(
        M.WeekendDaysWorked[i, t, p] for p in M.oneorzero_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])


model_phase1.weekend_subsets_4_3_con = pyo.Constraint(model_phase1.weekend_subsets_4_3_idx,
                                                      rule=weekend_subsets_4_3_rule)


# July 2019 - Not sure the weekend subset constraints as originally conceived are correct.
# Attempting new versions. Starting with this one to deal with tt8 (3433 12 hr) infeasibilities.

def weekend_subsets_4_3_rule2(M, i, t, w, e, d1, d2, d3):

    days_subset = [d1, d2, d3]
    return sum(M.DailyTourType[i, t, d, w] for d in days_subset) <= \
           sum(M.DailyTourType[i, t, d, w] for d in M.DAYS) \
           - sum(M.WeekendDaysWorked[i, t, p] for p in M.one_wkend_day[w, t, e]) \
           - 2 * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e]) \
           - sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])

model_phase1.weekend_subsets_4_3_con2 = pyo.Constraint(model_phase1.weekend_subsets_4_3_idx,
                                                       rule=weekend_subsets_4_3_rule2)


def weekend_subsets_5_4_rule2(M, i, t, w, e, d1, d2, d3, d4):

    days_subset = [d1, d2, d3, d4]
    return sum(M.DailyTourType[i, t, d, w] for d in days_subset) <= \
           sum(M.DailyTourType[i, t, d, w] for d in M.DAYS) \
           - sum(M.WeekendDaysWorked[i, t, p] for p in M.one_wkend_day[w, t, e]) \
           - 2 * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e]) \
           - sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])


model_phase1.weekend_subsets_5_4_con2 = pyo.Constraint(model_phase1.weekend_subsets_5_4_idx,
                                                       rule=weekend_subsets_5_4_rule2)



def weekend_subsets_4_4_idx_rule(M):
    index_list = []
    for (i, t) in M.okTourType:
        for w in M.WEEKS:
            for e in M.WEEKENDS:
                if M.tt_max_dys_weeks[t, w] == 4:
                    index_list.append((i, t, w, e, 2, 3, 4, 5))
                    index_list.append((i, t, w, e, 2, 3, 4, 6))
                    index_list.append((i, t, w, e, 2, 3, 5, 6))
                    index_list.append((i, t, w, e, 2, 4, 5, 6))
                    index_list.append((i, t, w, e, 3, 4, 5, 6))

    return index_list


model_phase1.weekend_subsets_4_4_idx = pyo.Set(dimen=8, initialize=weekend_subsets_4_4_idx_rule)


def weekend_subsets_4_4_rule(M, i, t, w, e, d1, d2, d3, d4):
    # total days in subset worked by all wkend patterns -  days worked by those with < 2 weekend days <= 4x-x or 3x
    # where x is number of weekend patterns with 2 wkend days
    # M.WeekendDaysWorked[i,t,p]
    days = [d1, d2, d3, d4]
    return sum(M.DailyTourType[i, t, d, w] for d in days) <= (len(days)) * sum(
        M.WeekendDaysWorked[i, t, p] for p in M.zero_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.one_wkend_day[w, t, e]) \
           + (len(days) - 2) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])


model_phase1.weekend_subsets_4_4_con = pyo.Constraint(model_phase1.weekend_subsets_4_4_idx,
                                                      rule=weekend_subsets_4_4_rule)


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
    # total days in subset worked by all wkend patterns -  days worked by those with < 2 weekend days <= 4x-x or 3x
    # where x is number of weekend patterns with 2 wkend days
    # M.WeekendDaysWorked[i,t,p]
    days = [d1, d2]
    return sum(M.DailyTourType[i, t, d, w] for d in days) <= (len(days)) * sum(
        M.WeekendDaysWorked[i, t, p] for p in M.oneorzero_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])


model_phase1.weekend_subsets_3_2_con = pyo.Constraint(model_phase1.weekend_subsets_3_2_idx,
                                                      rule=weekend_subsets_3_2_rule)

def weekend_subsets_3_2_rule2(M, i, t, w, e, d1, d2):

    days_subset = [d1, d2]
    return sum(M.DailyTourType[i, t, d, w] for d in days_subset) <= \
        sum(M.DailyTourType[i, t, d, w] for d in M.DAYS)  \
        - sum(M.WeekendDaysWorked[i, t, p] for p in M.one_wkend_day[w, t, e]) \
        - 2 * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e]) \
        - sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])


model_phase1.weekend_subsets_3_2_con2 = pyo.Constraint(model_phase1.weekend_subsets_3_2_idx,
                                                       rule=weekend_subsets_3_2_rule2)


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
    # total days in subset worked by all wkend patterns -  days worked by those with < 2 weekend days <= 4x-x or 3x
    # where x is number of weekend patterns with 2 wkend days
    # M.WeekendDaysWorked[i,t,p]
    days = [d1]
    return sum(M.DailyTourType[i, t, d, w] for d in days) <= (len(days)) * sum(
        M.WeekendDaysWorked[i, t, p] for p in M.oneorzero_wkend_day[w, t, e]) \
           + (len(days) - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])


model_phase1.weekend_subsets_2_1_con = pyo.Constraint(model_phase1.weekend_subsets_2_1_idx,
                                                      rule=weekend_subsets_2_1_rule)

def weekend_subsets_2_1_rule2(M, i, t, w, e, d1):

    days_subset = [d1]
    return sum(M.DailyTourType[i, t, d, w] for d in days_subset) <= \
        sum(M.DailyTourType[i, t, d, w] for d in M.DAYS)  \
        - sum(M.WeekendDaysWorked[i, t, p] for p in M.one_wkend_day[w, t, e]) \
        - 2 * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e]) \
        - sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])


model_phase1.weekend_subsets_2_1_con2 = pyo.Constraint(model_phase1.weekend_subsets_2_1_idx,
                                                       rule=weekend_subsets_2_1_rule2)


def DTT_mwdw_idx_rule(M):
    index_list = []
    for t in M.activeTT:
        numpats = M.num_mwdw_patterns[t]
        if numpats > 0:
            for i in M.WINDOWS:
                if (i, t) in M.okTourType:
                    for w in M.WEEKS:
                        index_list.append((i, t, w))
    return index_list


def DTT_mwdw_rule(M, i, t, w):
    return sum(M.DailyTourType[i,t,j,w] for j in M.DAYS) == sum(M.MultiWeekDaysWorked[i, t, p] * M.A_mwdw[t, p, w] for p in pyo.sequence(M.num_mwdw_patterns[t]))


model_phase1.DTT_mwdw_idx = pyo.Set(dimen=3, initialize=DTT_mwdw_idx_rule)

model_phase1.DTT_mwdw_con = pyo.Constraint(model_phase1.DTT_mwdw_idx,
                                                      rule=DTT_mwdw_rule)

def TT_mwdw_idx_rule(M):
    index_list = []
    for t in M.activeTT:
        numpats = M.num_mwdw_patterns[t]
        if numpats > 0:
            for i in M.WINDOWS:
                if (i, t) in M.okTourType:
                    index_list.append((i, t))
    return index_list


def TT_mwdw_rule(M, i, t):
    return M.TourType[i,t] == sum(M.MultiWeekDaysWorked[i, t, p] for p in pyo.sequence(M.num_mwdw_patterns[t]))


model_phase1.TT_mwdw_idx = pyo.Set(dimen=2, initialize=TT_mwdw_idx_rule)

model_phase1.TT_mwdw_con = pyo.Constraint(model_phase1.TT_mwdw_idx,
                                                       rule=TT_mwdw_rule)

# The following are various ad-hoc constraints that I was trying for fixing Phase 2 infeasibility problems.
# Any such constraints should include deactivation switches in solvemwts.py.

def ad_hoc_weekend_subsets_ttype7_idx_rule(M):
    index_list = []

    x = itertools.combinations(list(range(2, 7)), 2)
    y = itertools.combinations(list(range(9, 14)), 2)
    z = itertools.product(x, y)
    daypairs = [list(p) for p in z]

    for (i, t) in M.okTourType:
        if t == 6 or t == 7:
            for w in [1, 3]:
                for e in M.WEEKENDS:
                    # if M.tt_max_dys_weeks[t,w] ==  5:
                    for d in daypairs:
                        d1 = d[0][0]
                        d2 = d[0][1]
                        d3 = d[1][0] - 7
                        d4 = d[1][1] - 7
                        index_list.append((i, t, w, e, d1, d2, d3, d4))

    return index_list


model_phase1.ad_hoc_weekend_subsets_ttype7_idx = pyo.Set(dimen=8, initialize=ad_hoc_weekend_subsets_ttype7_idx_rule)


def ad_hoc_weekend_subsets_ttype7_rule(M, i, t, w, e, w1d1, w1d2, w2d1, w2d2):
    # ttype 7 is a half-time 8 (3-2 or 2-3)
    # Trying to avoid a forced 3-3 pattern due to weekend pattern
    # .

    wk1_days = [w1d1, w1d2]
    wk2_days = [w2d1, w2d2]

    subset_len = len(wk1_days) + len(wk2_days)

    return sum(M.DailyTourType[i, t, d, w] for d in wk1_days) + sum(
        M.DailyTourType[i, t, d, w + 1] for d in wk2_days) <= \
           subset_len * sum(M.WeekendDaysWorked[i, t, p] for p in M.zero_wkend_day[w, t, e]) \
           + (subset_len - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.one_wkend_day[w, t, e]) \
           + (subset_len - 2) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])


model_phase1.ad_hoc_weekend_subsets_ttype7_con = pyo.Constraint(model_phase1.ad_hoc_weekend_subsets_ttype7_idx,
                                                                rule=ad_hoc_weekend_subsets_ttype7_rule)


######################

def ad_hoc_weekend_subsets_ttype8_idx_rule(M):
    index_list = []

    w1 = itertools.combinations(list(range(2, 7)), 3)
    w2 = itertools.combinations(list(range(9, 14)), 3)
    w3 = itertools.combinations(list(range(16, 21)), 3)
    w4 = itertools.combinations(list(range(23, 28)), 3)
    z = itertools.product(w1, w2, w3, w4)
    daytrips = [list(p) for p in z]

    for (i, t) in M.okTourType:
        if t == 8:
            for w in [1]:
                for e in M.WEEKENDS:
                    # if M.tt_max_dys_weeks[t,w] ==  5:
                    for d in daytrips:
                        d1 = d[0][0]
                        d2 = d[0][1]
                        d3 = d[0][2]
                        d4 = d[1][0] - 7
                        d5 = d[1][1] - 7
                        d6 = d[1][2] - 7
                        d7 = d[2][0] - 14
                        d8 = d[2][1] - 14
                        d9 = d[2][2] - 14
                        d10 = d[3][0] - 21
                        d11 = d[3][1] - 21
                        d12 = d[3][2] - 21
                        index_list.append((i, t, w, e, d1, d2, d3,
                                           d4, d5, d6,
                                           d7, d8, d9,
                                           d10, d11, d12))

    return index_list


model_phase1.ad_hoc_weekend_subsets_ttype8_idx = pyo.Set(dimen=16, initialize=ad_hoc_weekend_subsets_ttype8_idx_rule)


def ad_hoc_weekend_subsets_ttype8_rule(M, i, t, w, e, w1d1, w1d2, w1d3, w2d1, w2d2, w2d3, w3d1, w3d2, w3d3, w4d1, w4d2,
                                       w4d3):
    # ttype 7 is a half-time 8 (3-2 or 2-3)
    # Trying to avoid a forced 3-3 pattern due to weekend pattern
    # .

    wk1_days = [w1d1, w1d2, w1d3]
    wk2_days = [w2d1, w2d2, w2d3]
    wk3_days = [w3d1, w3d2, w3d3]
    wk4_days = [w4d1, w4d2, w4d3]

    subset_len = len(wk1_days) + len(wk2_days) + len(wk3_days) + len(wk4_days)

    return sum(M.DailyTourType[i, t, d, w] for d in wk1_days) + sum(M.DailyTourType[i, t, d, w + 1] for d in wk2_days) \
           + sum(M.DailyTourType[i, t, d, w + 2] for d in wk3_days) + sum(
        M.DailyTourType[i, t, d, w + 3] for d in wk4_days) <= \
           subset_len * sum(M.WeekendDaysWorked[i, t, p] for p in M.zero_wkend_day[w, t, e]) \
           + (subset_len - 1) * sum(M.WeekendDaysWorked[i, t, p] for p in M.one_wkend_day[w, t, e]) \
           + (subset_len - 2) * sum(M.WeekendDaysWorked[i, t, p] for p in M.two_wkend_days[w, t, e])


model_phase1.ad_hoc_weekend_subsets_ttype8_con = pyo.Constraint(model_phase1.ad_hoc_weekend_subsets_ttype8_idx,
                                                                rule=ad_hoc_weekend_subsets_ttype8_rule)


def main():
    pass


if __name__ == '__main__':
    main()
