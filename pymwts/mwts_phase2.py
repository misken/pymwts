"""
Phase 2 for implicit multi-week tour scheduling model
"""

# Author: misken
# License: TBD

import sys

import pyomo.environ as pyo

import mwts_shared
from mwts_shared import epoch_to_tuple, epoch_increment

# from pyutilib.misc import import_file


# model = import_file('mwts_baseparams.py').model
model = pyo.AbstractModel()
model.name = "mwts_phase2"

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
    pyo.Param(within=pyo.PositiveIntegers, initialize=mwts_shared.n_prds_per_cycle_init)


# Temporal ranges used for various index sets
model.PERIODS = pyo.RangeSet(1, model.n_prds_per_day)
model.EPOCHS = pyo.RangeSet(1, model.n_prds_per_cycle)
model.WINDOWS = pyo.RangeSet(1, model.n_prds_per_day)
model.DAYS = pyo.RangeSet(1, model.n_days_per_week)
model.WEEKS = pyo.RangeSet(1, model.n_weeks)
model.WEEKENDS = pyo.RangeSet(1, 2)

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
# TODO: Some of these may be able to go away after I decide on possibly
# redundant constraints.

# 1a. Min and max number of days worked by week by tour type
model.tt_min_dys_weeks = pyo.Param(model.TTYPES,
                                          model.WEEKS, default=0.0)

model.tt_max_dys_weeks = pyo.Param(model.TTYPES,
                                          model.WEEKS, default=1e+6)

# 1b. Min and max number of days worked by cumulative weeks by tour type
model.tt_min_cumul_dys_weeks = pyo.Param(model.TTYPES,
                                                model.WEEKS, default=0.0)

model.tt_max_cumul_dys_weeks = pyo.Param(model.TTYPES,
                                                model.WEEKS, default=1e+6)

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
                                       model.TTYPES)

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







def activeWIN_init(M):
    return [i for i in M.WINDOWS]


model.activeWIN = pyo.Set(dimen=1, ordered=True, initialize=activeWIN_init)






# TODO: start windows ---------------------------------------------------------
# START WINDOWS - should these be tour type specific since allow start is tour type specific?
#
# To start with,I made them week specific but not sure if they should be. Seems they should be 
# period, day, tour type.
#
# Maybe best strategy is to do the variables and constraints first and then work backwards to
# define windows appropriately.
# -----------------------------------------------------------------------

model.g_start_window_width = pyo.Param()  # Width of start-time windows


# #/**** Beginning of each start window (in total periods from Sunday @ midnight)****/
# #pyo.Param b_window_wepoch{i in PERIODS,j in DAYS} := n_prds_per_day*(j-1)+i;
# #
# #/**** End of each start window (in total periods from Sunday @ midnight) ****/
# #pyo.Param e_window_wepoch{i in PERIODS,j in DAYS} :=
# # ( if n_prds_per_day*(j-1)+i+width <= n_prds_per_day*n_days then
# #    n_prds_per_day*(j-1)+i+width
# #  else
# #    (n_prds_per_day*(j-1)+i+width )-n_prds_per_day*n_days);


# if model.g_start_window_width > 0:
# This version spans multiple weeks
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

# Index: okStartWindowRoots is defined for all (tour type,shift length) pairs such that the shift length
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
                                          dimen=3, ordered=True,
                                          initialize=okStartWindowRoots_init)


# #
# #
# ## CHAINS - NON-OVERLAPPING SEQUENCES OF (i,j) period,day PAIRS THAT
# ##          CAN BE ISOLATED FOR COORDINATING ix AND DWT pyo.VarIABLES.
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
            for (p, q, r) in (M.okStartWindowRoots[t, k] - pyo.Set(initialize=[(i, j, w)])):
                if (i, j, w) not in M.PotentialStartWindow[p, q, r, k, t]:
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
            window_list.append((a, b, c))
        if m > 1:
            for (a, b, c) in M.linkspan[t, k, i, j, w, m - 1]:
                window_list.append((a, b, c))

    return window_list


model.linkspan = pyo.Set(model.link_idx, dimen=3,
                                ordered=True, initialize=linkspan_init)


# Phase 1 Shift variables --> Phase 2 parameters ------------------------------

model.Shift = pyo.Param(model.okShifts, default=0)

# Shift[i,j,w,k,t] = Number of shifts of length k starting in period i
# of day j in week w for a tour of type t


model.okTourType = model.WINDOWS * model.TTYPES


def TourType_idx_rule(M):
    return [(i, t) for i in M.WINDOWS
            for t in M.TTYPES
            if (i, t) in M.okTourType]


model.TourType_idx = pyo.Set(dimen=2, initialize=TourType_idx_rule)
model.TourType = pyo.Param(model.TourType_idx, default=0)


model.okTourTypeDay = model.WINDOWS * model.TTYPES * model.DAYS


def TourTypeDay_idx_rule(M):
    return [(i, t, j, w) for i in M.WINDOWS
            for t in M.TTYPES
            for j in M.DAYS
            for w in M.WEEKS
            if (i, t, j) in M.okTourTypeDay]


model.TourTypeDay_idx = pyo.Set(dimen=4, initialize=TourTypeDay_idx_rule)
model.TourTypeDay = pyo.Param(model.TourTypeDay_idx)


# TODO - modify this index when get width>0 working. These are windows, not period start times.
# Just using okShifts for now for w=0 case.  

def TourTypeDayShift_idx_rule(M):
    return [(i, t, k, j, w) for i in M.WINDOWS
            for t in M.TTYPES
            for k in M.tt_length_x[t]
            for j in M.DAYS
            for w in M.WEEKS
            if (i, t, j) in M.okTourTypeDay]


model.TourTypeDayShift_idx = pyo.Set(dimen=5, initialize=TourTypeDayShift_idx_rule)
model.TourTypeDayShift = pyo.Param(model.TourTypeDayShift_idx, default=0)


# #### Weekend Days off variables

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

model.WeekendDaysWorked = pyo.Param(model.weekend_days_worked_idx, default=0)


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
    pyo.Param(model.multiweekdaysworked_idx,
              within=pyo.NonNegativeIntegers)


model.n_tours = pyo.Param(within=pyo.PositiveIntegers)
model.TOURS = pyo.RangeSet(1, model.n_tours)

# The following two parameters are created after solving the Phase 1 model
# and stored in the Phase 2 dat file. They are created
# by mwts_utils.tour_WIN_TT_to_param.


model.WIN_x = pyo.Param(model.TOURS)
model.TT_x = pyo.Param(model.TOURS)


# Decision Variables ----------------------------------------------------------

def TourShifts_idx_rule(M):
    index_list = []
    for s in M.TOURS:
        for i in M.PERIODS:
            for j in M.DAYS:
                for w in M.WEEKS:
                    for t in [a for a in M.activeTT if a == M.TT_x[s]]:
                        for k in M.tt_length_x[t]:
                            for (p, d, q) in [(x, y, z) for (x, y, z) in M.okStartWindowRoots[t, k] if
                                              x == M.WIN_x[s] and x in M.activeWIN]:
                                if (i, j, w) in M.PotentialStartWindow[p, d, q, k, t]:
                                    index_list.append((s, i, j, w, k, t))

    # print("TourShift_idx", index_list)
    return index_list


model.TourShifts_idx = pyo.Set(dimen=6, initialize=TourShifts_idx_rule)
model.TourShifts = pyo.Var(model.TourShifts_idx, within=pyo.Boolean)


###
###var tourdof{l in 1..n_tours,d in 1..max_weekend_patterns,i in WINDOWS,t in okTTYPES :
###       t=TT_x[l] and i=WIN_x[l] and (i,t) in okTourType and d <= num_weekend_patterns[weekend_type[i,t],t] } binary ;
###    
###    # tourdof[l,d,i,t] = 1 if an ID[d,i,t] days off variable is assigned to tour l
###    #            = 0 otherwise
###    # tourdof ==> "binary days-off"  

def TourWkendDof_idx_rule(M):
    index_list = []
    for s in M.TOURS:
        for i in [a for a in M.activeWIN if a == M.WIN_x[s]]:
            for t in [a for a in M.activeTT if a == M.TT_x[s]]:
                for pattern in pyo.sequence(M.num_weekend_patterns[M.weekend_type[i, t], t]):
                    if (i, t) in M.okTourType:
                        index_list.append((s, pattern, i, t))

    return index_list


model.TourWkendDof_idx = pyo.Set(dimen=4, initialize=TourWkendDof_idx_rule)
model.TourWkendDof = pyo.Var(model.TourWkendDof_idx, within=pyo.Boolean)


def Tour_mwdwWkendDof_idx_rule(M):
    index_list = []
    for s in M.TOURS:
        for i in [a for a in M.activeWIN if a == M.WIN_x[s]]:
            for t in [a for a in M.activeTT if a == M.TT_x[s]]:
                for p1 in pyo.sequence(M.max_mwdw_patterns):
                    for p2 in pyo.sequence(M.num_weekend_patterns[M.weekend_type[i, t], t]):
                        if (i, t) in M.okTourType and p1 <= M.num_mwdw_patterns[t]:
                            index_list.append((s, p1, p2, i, t))

    return index_list


model.Tour_mwdwWkendDof_idx = pyo.Set(dimen=5, initialize=Tour_mwdwWkendDof_idx_rule)
model.Tour_mwdwWkendDof = pyo.Var(model.Tour_mwdwWkendDof_idx, within=pyo.Boolean)


# .......................................................Obj. Function

# ## Objective function is just the sum of the tourshift variables - which makes it a search for a feasible solution.
# ##
# ##minimize total_cost :
# ##
# ## sum {s in 1..n_tours,i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS, t in TTYPES, p in PERIODS, d in DAYS :
# ##       t=TT_x[s] and p=WIN_x[s] and k in tt_length_x[t] and (i,j,k,t) in ok_shifts and
# ##       (i,j) in okWindowWepochs[p,d,k,t]} tourshift[s,i,j,w,k,t,p,d] ;

# Objective function
def objective_rule(M):
    obj1 = sum(M.TourShifts[s, i, j, w, k, t] for (s, i, j, w, k, t) in M.TourShifts_idx)
    return obj1


model.total_num_tours = pyo.Objective(rule=objective_rule, sense=pyo.minimize)


# One weekend days off pattern for each tour


def Tour_WkendDof_conservation_idx_rule(M):
    index_list = []
    for i in M.activeWIN:
        for t in M.activeTT:
            for pattern in pyo.sequence(M.max_weekend_patterns):
                if (i, t) in M.okTourType and pattern <= M.num_weekend_patterns[M.weekend_type[i, t], t]:
                    for s in M.TOURS:
                        if i == M.WIN_x[s] and t == M.TT_x[s]:
                            index_list.append((pattern, i, t))
                            break
    return index_list


def Tour_mwdw_WkendDof_conservation_idx_rule(M):
    index_list = []
    for i in M.activeWIN:
        for t in M.activeTT:
            for p1 in pyo.sequence(M.max_mwdw_patterns):
                for p2 in pyo.sequence(M.max_weekend_patterns):
                    if (i, t) in M.okTourType and p1 <= M.num_mwdw_patterns[t] and p2 <= M.num_weekend_patterns[
                        M.weekend_type[i, t], t]:
                        for s in M.TOURS:
                            if i == M.WIN_x[s] and t == M.TT_x[s]:
                                index_list.append((p1, p2, i, t))
                                break
    return index_list


model.Tour_WkendDof_conservation_idx = pyo.Set(dimen=3, initialize=Tour_WkendDof_conservation_idx_rule)
model.Tour_mwdw_WkendDof_conservation_idx = pyo.Set(dimen=4, initialize=Tour_mwdw_WkendDof_conservation_idx_rule)


def Tour_WkendDof_conservation_rule(M, pattern, i, t):
    """
    Sum over the tours within each (i,t) and make sure they add up to the WeekendDaysWorked variables
    """
    return (sum(M.TourWkendDof[s, pattern, i, t] for s in M.TOURS if i == M.WIN_x[s] and t == M.TT_x[s]) ==
            M.WeekendDaysWorked[i, t, pattern])


model.Tour_WkendDof_conservation = pyo.Constraint(model.Tour_WkendDof_conservation_idx,
                                                         rule=Tour_WkendDof_conservation_rule)


def Tour_mwdw_WkendDof_conservation_rule(M, p1, p2, i, t):
    """
    Sum over the tours within each (i,t) and make sure they add up to the MultiWeekDaysWorked variables
    """
    return (sum(M.Tour_mwdwWkendDof[s, p1, p2, i, t] for s in M.TOURS if i == M.WIN_x[s] and t == M.TT_x[s]) ==
            M.MultiWeekDaysWorked[i, t, p1, p2])


model.Tour_mwdw_WkendDof_conservation = pyo.Constraint(model.Tour_mwdw_WkendDof_conservation_idx,
                                                              rule=Tour_mwdw_WkendDof_conservation_rule)


def OnePatternPerTourShift_idx_rule(M):
    """
    Index based on activeTT so that we can create schedules for any subset
    of tour types used in the Phase1 problem.
    """
    index_list = []
    for s in M.TOURS:
        if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN:
            index_list.append(s)

    return index_list


model.OnePatternPerTourShift_idx = pyo.Set(dimen=1, initialize=OnePatternPerTourShift_idx_rule)


def OnePatternPerTourShift_rule(M, s):
    return sum(M.TourWkendDof[s, pattern, M.WIN_x[s], M.TT_x[s]]
               for pattern in
               pyo.sequence(M.num_weekend_patterns[M.weekend_type[M.WIN_x[s], M.TT_x[s]], M.TT_x[s]])) == 1


model.OnePatternPerTourShift = pyo.Constraint(model.OnePatternPerTourShift_idx,
                                                     rule=OnePatternPerTourShift_rule)


def One_mwdwPatternPerTourShift_idx_rule(M):
    """
    Index based on activeTT so that we can create schedules for any subset
    of tour types used in the Phase1 problem.
    """
    index_list = []
    for s in M.TOURS:
        if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN:
            index_list.append(s)

    return index_list


model.One_mwdwPatternPerTourShift_idx = pyo.Set(dimen=1, initialize=One_mwdwPatternPerTourShift_idx_rule)


def One_mwdwPatternPerTourShift_rule(M, s):
    return sum(M.Tour_mwdwWkendDof[s, p1, p2, M.WIN_x[s], M.TT_x[s]]
               for p1 in pyo.sequence(M.num_mwdw_patterns[M.TT_x[s]])
               for p2 in
               pyo.sequence(M.num_weekend_patterns[M.weekend_type[M.WIN_x[s], M.TT_x[s]], M.TT_x[s]])) == 1


model.One_mwdwPatternPerTourShift = pyo.Constraint(model.One_mwdwPatternPerTourShift_idx,
                                                          rule=One_mwdwPatternPerTourShift_rule)


# All days off patterns assigned

# subject to conserve_ID {d in 1..max_weekend_patterns,i in WINDOWS,t in okTTYPES: d <= num_weekend_patterns[weekend_type[i,t],t] and (i,t) in okTourType} :
# sum{s in 1..n_tours : t=TT_x[s] and i=WIN_x[s]} 
#        tourdof[s,d,i,t]= WeekendDaysWorked[d,i,t];

# 


# Weekend days worked for each tour determines number of weekend shifts assigned to each tour

###subject to tourshift_tourdof_integration1{s in 1..n_tours,i in WINDOWS,j in DAYS, w in WEEKS,t in TTYPES: 
###    t=TT_x[s] and i=WIN_x[s] and j in weekend[i,TT_x[s]]} : 

###     sum{k in tt_length_x[t],(p,d) in okWindowWepochs[i,j,k,t] } tourshift[s,p,d,w,k,t,i,j] =
###        (sum{g in 1..num_weekend_patterns[weekend_type[i,t],t]}tourdof[s,g,i,t]*A_wkend_days[g,j,w,t,weekend_type[i,t]]);

def Tour_ShiftWkendDof_integration1_idx_rule(M):
    index_list = []
    for s in M.TOURS:
        for i in [a for a in M.activeWIN if a == M.WIN_x[s]]:
            for t in [a for a in M.activeTT if a == M.TT_x[s]]:
                for d in M.weekend[i, t]:
                    for w in M.WEEKS:
                        index_list.append((s, i, d, w, t))
    return index_list


model.Tour_ShiftWkendDof_integration1_idx = \
    pyo.Set(dimen=5, initialize=Tour_ShiftWkendDof_integration1_idx_rule)


def weekend_integration_1_SS_idx_rule(M):
    return [(j, w, i, t) for j in M.DAYS
            for w in M.WEEKS
            for i in M.activeWIN
            for t in M.activeTT
            if j in M.weekend[i, t] and 1 in M.weekend[i, t] and 7 in M.weekend[i, t] \
            and (i, t, j) in M.okTourTypeDay]


model.weekend_integration_1_SS_idx = pyo.Set(dimen=4, initialize=weekend_integration_1_SS_idx_rule)


def action_check_WeekendDaysWorked_TourTypeDay_rule(M, j, w, i, t):
    val_w = pyo.value(
        sum(M.A_wkend_days[p, j, w, t, 1] * M.WeekendDaysWorked[i, t, p] for p in pyo.sequence(M.num_weekend_patterns[1, t])))
    val_d = M.TourTypeDay[i, t, j, w]

    if val_w == val_d:
        check = 'OK'
    else:
        check = 'BAD'

    msg = '({},{},{},{}) w={} d={} --> {}'.format(i, j, w, t, val_w, val_d, check)
    if check == 'BAD':
        print(msg)


model.action_check_WeekendDaysWorked_TourTypeDay = \
    pyo.BuildAction(model.weekend_integration_1_SS_idx,
                    rule=action_check_WeekendDaysWorked_TourTypeDay_rule)


def Tour_ShiftWkendDof_integration1_rule(M, s, i, j, w, t):
    return sum(M.TourShifts[s, p, d, q, k, t] for k in M.tt_length_x[t] for (p, d, q) in
               M.PotentialStartWindow[i, j, w, k, t]) == \
           sum(M.A_wkend_days[pattern, j, w, t, M.weekend_type[i, t]] * M.TourWkendDof[s, pattern, i, t] for pattern in
               pyo.sequence(M.num_weekend_patterns[M.weekend_type[i, t], t]))


model.Tour_ShiftWkendDof_integration1 = \
    pyo.Constraint(model.Tour_ShiftWkendDof_integration1_idx,
                   rule=Tour_ShiftWkendDof_integration1_rule)


# No more than one shift worked per day

# subject to tours_daily{s in 1..n_tours,j in DAYS,w in WEEKS} :
#  sum{k in tt_length_x[TT_x[s]],(p,d) in okWindowWepochs[WIN_x[s],j,k,TT_x[s]]}
#   tourshift[s,p,d,w,k,TT_x[s],WIN_x[s],j] <= 1;


def Tours_Daily_idx_rule(M):
    return [(s, j, w) for s in M.TOURS
            for j in M.DAYS
            for w in M.WEEKS
            if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN]


model.Tours_Daily_idx = pyo.Set(dimen=3, initialize=Tours_Daily_idx_rule)


def Tours_Daily_rule(M, s, j, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]]
               for k in M.tt_length_x[M.TT_x[s]] for (p, d, q) in
               M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]) <= 1


model.Tours_Daily = pyo.Constraint(model.Tours_Daily_idx,
                                          rule=Tours_Daily_rule)


## Needs to be generalized for start windows
# subject to tours_daily_conservation2{i in PERIODS, j in DAYS,w in WEEKS, k in LENGTHS, t in TTYPES: k in tt_length_x[t] and (i,j,k,t) in ok_shifts} :
#  sum{s in 1..n_tours: TT_x[s]=t and WIN_x[s]=i }
#   tourshift[s,i,j,w,k,t,WIN_x[s],j] = Shift[i,j,w,k,t];  
# 
#
#
## All the shifts need to get assigned to the tours
#
# /*
# subject to tours_daily_conservation2{p in PERIODS, d in DAYS,w in WEEKS, k in LENGTHS, t in TTYPES: k in tt_length_x[t]} :
#  sum{s in 1..n_tours, (i,j) in okWindowWepochs[p,d,k,t]: TT_x[s]=t and WIN_x[s]=p }
#   tourshift[s,i,j,w,k,t,p,d] = sum{(i,j) in okWindowWepochs[p,d,k,t] : (i,j,k,t) in ok_shifts}Shift[i,j,w,k,t]; 
# */

def Tours_Daily_conservation_idx_rule(M):
    index_list = []
    for p in M.PERIODS:
        for d in M.DAYS:
            for w in M.WEEKS:
                for t in M.activeTT:
                    for k in M.tt_length_x[t]:
                        if (p, d, w) in M.okStartWindowRoots[t, k]:
                            for s in M.TOURS:
                                if (p == M.WIN_x[s] and p in M.activeWIN and t == M.TT_x[s]):
                                    index_list.append((p, d, w, k, t))
                                    break
    return index_list


model.Tours_Daily_conservation_idx = \
    pyo.Set(dimen=5, initialize=Tours_Daily_conservation_idx_rule)


def Tours_Daily_conservation_rule(M, p, d, q, k, t):
    return sum(M.TourShifts[s, i, j, w, k, t]
               for s in M.TOURS for (i, j, w) in M.PotentialStartWindow[p, d, q, k, t] if
               p == M.WIN_x[s] and t == M.TT_x[s]) \
           == sum(M.Shift[i, j, w, k, t] for (i, j, w) in M.PotentialStartWindow[p, d, q, k, t])


model.Tours_Daily_conservation = \
    pyo.Constraint(model.Tours_Daily_conservation_idx,
                   rule=Tours_Daily_conservation_rule)


# =============== Tours_Weekly Lower and Upper Bounds on Days Worked ==================================================

# == For each (tour type) make sure the number of shifts assigned each week
# == satisfies associated lower and upper bounds.

## Tour type specific bounds on number of days worked over the weeks (both cumulative and non-cumulative)
# subject to tours_weekly_LB{t in 1..n_tours,w in WEEKS} :
#  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_min_dys_weeks[TT_x[t],w];
#   
# subject to tours_weekly_UB{t in 1..n_tours,w in WEEKS} :
#  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_max_dys_weeks[TT_x[t],w];

def Tours_Weekly_LB_idx_rule(M):
    return [(s, w) for s in M.TOURS for w in M.WEEKS if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN]


model.Tours_Weekly_LB_idx = pyo.Set(dimen=2,
                                           initialize=Tours_Weekly_LB_idx_rule)


def Tours_Weekly_LB_rule(M, s, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]]
               for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]]
               for (p, d, q) in M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]) >= M.tt_min_dys_weeks[
               M.TT_x[s], w]


# Instead of LB and UB, now that we have mwdw variables, we know exactly how many shifts worked each
# week for each tour.
def Tours_Weekly_rule(M, s, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]]
               for j in M.DAYS
               for k in M.tt_length_x[M.TT_x[s]]
               for (p, d, q) in M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]) == \
           sum(M.Tour_mwdwWkendDof[s, p1, p2, M.WIN_x[s], M.TT_x[s]] * M.A_mwdw[
               M.TT_x[s], p1, w] for p1 in
               pyo.sequence(M.num_mwdw_patterns[M.TT_x[s]]) for p2 in
               pyo.sequence(M.num_weekend_patterns[M.weekend_type[M.WIN_x[s], M.TT_x[s]], M.TT_x[s]]))


model.Tours_Weekly_LB = pyo.Constraint(model.Tours_Weekly_LB_idx,
                                              rule=Tours_Weekly_LB_rule)

model.Tours_Weekly = pyo.Constraint(model.Tours_Weekly_LB_idx,
                                           rule=Tours_Weekly_rule)


def Tours_Weekly_UB_idx_rule(M):
    return [(s, w) for s in M.TOURS for w in M.WEEKS if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN]


model.Tours_Weekly_UB_idx = pyo.Set(dimen=2,
                                           initialize=Tours_Weekly_UB_idx_rule)


def Tours_Weekly_UB_rule(M, s, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]]
               for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p, d, q) in
               M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]) <= M.tt_max_dys_weeks[M.TT_x[s], w]


model.Tours_Weekly_UB = pyo.Constraint(model.Tours_Weekly_UB_idx,
                                              rule=Tours_Weekly_UB_rule)


# =============== Tours_Total Lower and Upper Bounds on Days Worked ======================================================

# == For each (tour type) make sure the number of shifts assigned over cumulative
# == weeks satisfies associated lower and upper bounds.

# Seems like the following two constraints need to be generalized for > 2 weeks.
#
# subject to tours_tot_LB{t in 1..n_tours} :
#  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_min_cumul_dys_weeks[TT_x[t],n_weeks];
#   
# subject to tours_tot_UB{t in 1..n_tours} :
#  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_max_cumul_dys_weeks[TT_x[t],n_weeks];


def Tours_Total_LB_idx_rule(M):
    return [(s, w) for s in M.TOURS for w in M.WEEKS if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN]


model.Tours_Total_LB_idx = pyo.Set(dimen=2,
                                          initialize=Tours_Total_LB_idx_rule)


def Tours_Total_LB_rule(M, s, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]]
               for W in pyo.sequence(w) for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p, d, q) in
               M.PotentialStartWindow[M.WIN_x[s], j, W, k, M.TT_x[s]]) >= M.tt_min_cumul_dys_weeks[M.TT_x[s], w]


model.Tours_Total_LB = pyo.Constraint(model.Tours_Total_LB_idx,
                                             rule=Tours_Total_LB_rule)


def Tours_Total_UB_idx_rule(M):
    return [(s, w) for s in M.TOURS for w in M.WEEKS if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN]


model.Tours_Total_UB_idx = pyo.Set(dimen=2,
                                          initialize=Tours_Total_UB_idx_rule)


def Tours_Total_UB_rule(M, s, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]]
               for W in pyo.sequence(w)
               for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]]
               for (p, d, q) in M.PotentialStartWindow[M.WIN_x[s], j, W, k, M.TT_x[s]]) <= M.tt_max_cumul_dys_weeks[
               M.TT_x[s], w]


model.Tours_Total_UB = pyo.Constraint(model.Tours_Total_UB_idx,
                                             rule=Tours_Total_UB_rule)


# =============== Tours_Shiftlen_Weekly Lower and Upper Bounds on Days Worked ================================================

# == For each (tour type, shift length) make sure the number of shifts assigned each week
# == satisfies associated lower and upper bounds.

# Tour type and shift length specific bounds on number of days worked over the weeks (both cumulative and non-cumulative)
# subject to tours_shiftlen_weekly_LB{t in 1..n_tours,k in tt_length_x[TT_x[t]],w in WEEKS} :
#  sum{d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_shiftlen_min_dys_weeks[TT_x[t],k,w];
#   
# subject to tours_shiftlen_weekly_UB{t in 1..n_tours,k in tt_length_x[TT_x[t]],w in WEEKS} :
#  sum{d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_shiftlen_max_dys_weeks[TT_x[t],k,w];
#   

# and M.tt_shiftlen_min_dys_weeks[M.TT_x[s],k,w] > 0
def Tours_Shiftlen_Weekly_LB_idx_rule(M):
    index_list = []
    for s in M.TOURS:
        for k in M.tt_length_x[M.TT_x[s]]:
            for w in M.WEEKS:
                if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN:
                    if any((s, p, d, q, k, M.TT_x[s]) in model.TourShifts_idx for j in M.DAYS for (p, d, q)
                           in M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]):
                        index_list.append((s, k, w))

    return index_list


#    return [(s,k,w) for s in M.TOURS for k in M.tt_length_x[M.TT_x[s]] for w in M.WEEKS if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN]

model.Tours_Shiftlen_Weekly_LB_idx = pyo.Set(dimen=3,
                                                    initialize=Tours_Shiftlen_Weekly_LB_idx_rule)


def Tours_Shiftlen_Weekly_LB_rule(M, s, k, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]]
               for j in M.DAYS for (p, d, q) in M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]) >= \
           M.tt_shiftlen_min_dys_weeks[M.TT_x[s], k, w]


model.Tours_Shiftlen_Weekly_LB = \
    pyo.Constraint(model.Tours_Shiftlen_Weekly_LB_idx,
                   rule=Tours_Shiftlen_Weekly_LB_rule)


def Tours_Shiftlen_Weekly_UB_idx_rule(M):
    # return [(s,k,w) for s in M.TOURS for k in M.tt_length_x[M.TT_x[s]] for w in M.WEEKS if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN]
    index_list = []
    for s in M.TOURS:
        for k in M.tt_length_x[M.TT_x[s]]:
            for w in M.WEEKS:
                if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN:
                    if any((s, p, d, q, k, M.TT_x[s]) in model.TourShifts_idx for j in M.DAYS for (p, d, q) in
                           M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]):
                        index_list.append((s, k, w))

    return index_list


model.Tours_Shiftlen_Weekly_UB_idx = \
    pyo.Set(dimen=3, initialize=Tours_Shiftlen_Weekly_UB_idx_rule)


def Tours_Shiftlen_Weekly_UB_rule(M, s, k, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]]
               for j in M.DAYS for (p, d, q) in M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]) <= \
           M.tt_shiftlen_max_dys_weeks[M.TT_x[s], k, w]


model.Tours_Shiftlen_Weekly_UB = pyo.Constraint(model.Tours_Shiftlen_Weekly_UB_idx,
                                                       rule=Tours_Shiftlen_Weekly_UB_rule)


# =============== Tours_Shiftlen_Total Lower and Upper Bounds on Days Worked ==============================================

# == For each (tour type, shift length) make sure the number of shifts assigned over cumulative
# == weeks satisfies associated lower and upper bounds.

# subject to tours_shiftlen_tot_LB{t in 1..n_tours,k in tt_length_x[TT_x[t]]} :
#  sum{w in WEEKS,d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_shiftlen_min_cumul_dys_weeks[TT_x[t],k,n_weeks];
#   
# subject to tours_shiftlen_tot_UB{t in 1..n_tours,k in tt_length_x[TT_x[t]]} :
#  sum{w in WEEKS,d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_shiftlen_max_cumul_dys_weeks[TT_x[t],k,n_weeks]; 


def Tours_Shiftlen_Total_LB_idx_rule(M):
    index_list = []
    for s in M.TOURS:
        for k in M.tt_length_x[M.TT_x[s]]:
            for w in M.WEEKS:
                if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN:
                    if any((s, p, d, q, k, M.TT_x[s]) in model.TourShifts_idx for j in M.DAYS for (p, d, q)
                           in M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]):
                        index_list.append((s, k, w))

    return index_list


model.Tours_Shiftlen_Total_LB_idx = \
    pyo.Set(dimen=3,
            initialize=Tours_Shiftlen_Total_LB_idx_rule)


def Tours_Shiftlen_Total_LB_rule(M, s, k, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]]
               for W in range(1, w + 1) for j in M.DAYS for (p, d, q) in
               M.PotentialStartWindow[M.WIN_x[s], j, W, k, M.TT_x[s]]) >= M.tt_shiftlen_min_cumul_dys_weeks[
               M.TT_x[s], k, w]


model.Tours_Shiftlen_Total_LB = \
    pyo.Constraint(model.Tours_Shiftlen_Total_LB_idx,
                   rule=Tours_Shiftlen_Total_LB_rule)


def Tours_Shiftlen_Total_UB_idx_rule(M):
    index_list = []
    for s in M.TOURS:
        for k in M.tt_length_x[M.TT_x[s]]:
            for w in M.WEEKS:
                if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN:
                    if any((s, p, d, q, k, M.TT_x[s]) in model.TourShifts_idx for j in M.DAYS for (p, d, q)
                           in M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]):
                        index_list.append((s, k, w))

    return index_list


model.Tours_Shiftlen_Total_UB_idx = \
    pyo.Set(dimen=3,
            initialize=Tours_Shiftlen_Total_UB_idx_rule)


def Tours_Shiftlen_Total_UB_rule(M, s, k, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]]
               for W in range(1, w + 1) for j in M.DAYS for (p, d, q) in
               M.PotentialStartWindow[M.WIN_x[s], j, W, k, M.TT_x[s]]) <= M.tt_shiftlen_max_cumul_dys_weeks[
               M.TT_x[s], k, w]


model.Tours_Shiftlen_Total_UB = pyo.Constraint(model.Tours_Shiftlen_Total_UB_idx,
                                                      rule=Tours_Shiftlen_Total_UB_rule)


# =============== Tours_Weekly Lower and Upper Bounds on Periods Worked ===========================================

# == For each (tour type) make sure the number of periods assigned each week
# == satisfies associated lower and upper bounds.

# Tour type specific bounds on number of periods worked over the weeks (both cumulative and non-cumulative)
# subject to tours_weekly_prds_LB{t in 1..n_tours,w in WEEKS} :
#  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] >= tt_min_prds_weeks[TT_x[t],w];
#   
# subject to tours_weekly_prds_UB{t in 1..n_tours,w in WEEKS} :
#  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] <= tt_max_prds_weeks[TT_x[t],w];
#   

def Tours_Weekly_Prds_LB_idx_rule(M):
    return [(s, w) for s in M.TOURS for w in M.WEEKS if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN]


model.Tours_Weekly_Prds_LB_idx = \
    pyo.Set(dimen=2,
            initialize=Tours_Weekly_Prds_LB_idx_rule)


def Tours_Weekly_Prds_LB_rule(M, s, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]] * M.lengths[k]
               for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p, d, q) in
               M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]) >= M.tt_min_prds_weeks[M.TT_x[s], w]


model.Tours_Weekly_Prds_LB = pyo.Constraint(model.Tours_Weekly_Prds_LB_idx,
                                                   rule=Tours_Weekly_Prds_LB_rule)


def Tours_Weekly_Prds_UB_idx_rule(M):
    return [(s, w) for s in M.TOURS for w in M.WEEKS if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN]


model.Tours_Weekly_Prds_UB_idx = \
    pyo.Set(dimen=2,
            initialize=Tours_Weekly_Prds_UB_idx_rule)


def Tours_Weekly_Prds_UB_rule(M, s, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]] * M.lengths[k]
               for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p, d, q) in
               M.PotentialStartWindow[M.WIN_x[s], j, w, k, M.TT_x[s]]) <= M.tt_max_prds_weeks[M.TT_x[s], w]


model.Tours_Weekly_Prds_UB = \
    pyo.Constraint(model.Tours_Weekly_Prds_UB_idx,
                   rule=Tours_Weekly_Prds_UB_rule)


# =============== Tours_Total Lower and Upper Bounds on Periods Worked ===========================================

# == For each (tour type) make sure the number of periods assigned over cumulative weeks
# == satisfies associated lower and upper bounds.

# subject to tours_tot_prds_LB{t in 1..n_tours} :
#  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] >= tt_min_cumul_prds_weeks[TT_x[t],n_weeks];
#   
# subject to tours_tot_prds_UB{t in 1..n_tours} :
#  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] <= tt_max_cumul_prds_weeks[TT_x[t],n_weeks];  

def Tours_Total_Prds_LB_idx_rule(M):
    return [(s, w) for s in M.TOURS for w in M.WEEKS if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN]


model.Tours_Total_Prds_LB_idx = \
    pyo.Set(dimen=2,
            initialize=Tours_Total_Prds_LB_idx_rule)


def Tours_Total_Prds_LB_rule(M, s, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]] * M.lengths[k]
               for W in range(1, w + 1) for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p, d, q) in
               M.PotentialStartWindow[M.WIN_x[s], j, W, k, M.TT_x[s]]) >= M.tt_min_cumul_prds_weeks[M.TT_x[s], w]


model.Tours_Total_Prds_LB = pyo.Constraint(model.Tours_Total_Prds_LB_idx, rule=Tours_Total_Prds_LB_rule)


def Tours_Total_Prds_UB_idx_rule(M):
    return [(s, w) for s in M.TOURS for w in M.WEEKS if M.TT_x[s] in M.activeTT and M.WIN_x[s] in M.activeWIN]


model.Tours_Total_Prds_UB_idx = \
    pyo.Set(dimen=2,
            initialize=Tours_Total_Prds_UB_idx_rule)


def Tours_Total_Prds_UB_rule(M, s, w):
    return sum(M.TourShifts[s, p, d, q, k, M.TT_x[s]] * M.lengths[k]
               for W in range(1, w + 1) for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p, d, q) in
               M.PotentialStartWindow[M.WIN_x[s], j, W, k, M.TT_x[s]]) <= M.tt_max_cumul_prds_weeks[M.TT_x[s], w]


model.Tours_Total_Prds_UB = \
    pyo.Constraint(model.Tours_Total_Prds_UB_idx,
                   rule=Tours_Total_Prds_UB_rule)


# Ensures that each day worked in tour l is assigned exactly one shift
# of an appropriate length and from an appropriate start window.


# subject to tot_conserve_shifts
#   :
#  sum {s in 1..n_tours,i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS, t in TTYPES, p in PERIODS, d in DAYS :
#       t=TT_x[s] and p=WIN_x[s] and k in tt_length_x[t] and (i,j,k,t) in ok_shifts and 
#       (i,j) in okWindowWepochs[p,d,k,t]} tourshift[s,i,j,w,k,t,p,d] =
#   sum{m in PERIODS,n in DAYS,z in WEEKS, o in LENGTHS, t in TTYPES : 
#   (m,n,o,t) in ok_shifts and Shift[m,n,z,o,t]>0} Shift[m,n,z,o,t];

# def Tours_Total_conservation_rule(M):
#    
#    return sum(M.TourShifts[s,i,j,m,k,t,p,d] 
#        for s in M.TOURS 
#        for k in M.LENGTHS
#        for t in M.activeTT
#        for p in M.PERIODS
#        for d in M.DAYS
#        for w in M.WEEKS
#        for (i,j,m) in M.PotentialStartWindow[p,d,w,k,t] if p == M.WIN_x[s] and t == M.TT_x[s]) \
#                 <= sum(M.Shift[i,j,w,k,t] for (i,j,m) in M.PotentialStartWindow[p,d,w,k,t] if (i,j,w,k,t) in M.okShifts)
#                                    
# model.Tours_Total_conservation = pyo.Constraint(model.Tours_Total_conservation_idx,rule=Tours_Total_conservation_rule)


def main():
    pass


if __name__ == '__main__':
    main()
