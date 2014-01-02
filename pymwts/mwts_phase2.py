#-------------------------------------------------------------------------------
# Name:        mwts_phase2.py
# Purpose:
#
# Author:      isken
#
# Created:     27/11/2011
# Copyright:   (c) isken 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import sys

from coopr.pyomo import *
from pyutilib.misc import import_file

from mwts_utils import *


#model_phase2 = import_file('mwts_baseparams.py').model
model_phase2 = AbstractModel()
model_phase2.name = "mwts_phase2"

#### General parameters

infinity = float('inf')


model_phase2.n_prds_per_day = Param(within=PositiveIntegers)
model_phase2.n_days_per_week = Param(within=PositiveIntegers)
model_phase2.n_weeks = Param(within=PositiveIntegers)



def n_prds_per_week_init(M):
    return M.n_days_per_week()*M.n_prds_per_day()

model_phase2.n_prds_per_week = Param(within=PositiveIntegers, initialize=n_prds_per_week_init)

def n_prds_per_cycle_init(M):
    return M.n_weeks()*M.n_days_per_week()*M.n_prds_per_day()

model_phase2.n_prds_per_cycle = Param(within=PositiveIntegers, initialize=n_prds_per_cycle_init)

# For range sets, if start omitted, assumed range is 1..args[0]
model_phase2.PERIODS = RangeSet(1,model_phase2.n_prds_per_day)
model_phase2.WINDOWS = RangeSet(1,model_phase2.n_prds_per_day)
model_phase2.DAYS = RangeSet(1,model_phase2.n_days_per_week)
model_phase2.WEEKS = RangeSet(1,model_phase2.n_weeks)
model_phase2.WEEKENDS = RangeSet(1,2)

model_phase2.g_period = Param(model_phase2.PERIODS, model_phase2.DAYS, model_phase2.WEEKS, initialize=g_period_init)

# Chain related utility functions - I moved these to mwts_utils.
#def g_period_init(M, i, j, w):
#    return ((w-1)*M.n_days_per_week()*M.n_prds_per_day() + 
#        (j-1)*M.n_prds_per_day() + i)
#    
#model_phase2.g_period = Param(model_phase2.PERIODS, model_phase2.DAYS, model_phase2.WEEKS, initialize=g_period_init)
#
#def period_increment(M, i, j, w, incr):
#    p = M.g_period[i,j,w]
#    if (p + incr <= M.n_prds_per_cycle):
#        return p + incr
#    else:
#        return p + incr - M.n_prds_per_cycle
#
#def g_period_increment(M, p, incr):
#    if (p + incr <= M.n_prds_per_cycle):
#        return p + incr
#    else:
#        return p + incr - M.n_prds_per_cycle 
#    
#def g_period_difference(M, b_prd, e_prd):
#    if (e_prd >= b_prd):
#        return e_prd - b_prd + 1
#    else:
#        return M.n_prds_per_cycle + e_prd - b_prd + 1 
#    
#def g_prd_to_tuple(M, p):
##    param which_prd{p in 1..(n_days+1)*n_prds_per_day} :=
##   p-n_prds_per_day*(ceil(p/n_prds_per_day-1));
##
##param which_day{p in 1..(n_days+1)*n_prds_per_day} :=
##   (if p>n_prds_per_day*n_days then 1 else 1+ceil(p/n_prds_per_day-1));
#    
#    n_week = ((p-1) // M.n_prds_per_week) + 1
#    prds_remainder = p - (n_week-1) * M.n_prds_per_week
#    if (prds_remainder == 0):
#        n_day = 1
#    else:
#        n_day = ((prds_remainder-1) // M.n_prds_per_day) + 1
#                 
#    prds_remainder = prds_remainder - (n_day - 1) * M.n_prds_per_day 
#    if (prds_remainder == 0):
#        n_period = 1
#    else:
#        n_period = prds_remainder
    
      

model_phase2.bins = model_phase2.PERIODS * model_phase2.DAYS * model_phase2.WEEKS

def oneweek_bins_init(M):
    return [(i,j) for i in M.PERIODS
                  for j in M.DAYS]
    
model_phase2.oneweek_bins = Set(dimen=2, ordered=True, initialize=oneweek_bins_init)


#### Tour type related parameters

#-- Shift Lengths
model_phase2.n_lengths = Param(within=PositiveIntegers)    # Number of shift lengths
model_phase2.LENGTHS = RangeSet(1,model_phase2.n_lengths)
model_phase2.lengths = Param(model_phase2.LENGTHS)  # Vector of shift lengths

#-- Tour Types
model_phase2.n_tts = Param(within=PositiveIntegers)  # Number of different tour types
model_phase2.TTYPES = RangeSet(1,model_phase2.n_tts)
model_phase2.tt_length_x = Set(model_phase2.TTYPES,ordered=True,)  # Set of allowable length indices by tour type


#-- Weekend patterns

def maxwkend_init(M):
    maxpats = 2**(2*M.n_weeks)
    return maxpats

model_phase2.max_weekend_patterns = Param(initialize=maxwkend_init)

model_phase2.num_weekend_patterns = Param(model_phase2.WEEKENDS,model_phase2.TTYPES)       # Number of weekends worked patterns


## param A[i,j,w,t,e] = 1 if weekend pattern i calls for work on day j of week k for tour type t having weekend type e and 0 otherwise

def A_idx_rule(M):
    return [(i,j,w,t,e) for i in range(1,M.max_weekend_patterns + 1)
                 for j in M.DAYS
                 for w in range(1,M.n_weeks + 1)
                 for t in range(1,M.n_tts + 1)
                 for e in range(1,3) if i <= M.num_weekend_patterns[e,t]]

model_phase2.A_idx = Set(dimen=5,ordered=True,initialize=A_idx_rule)
model_phase2.A = Param(model_phase2.A_idx, default=0.0)
#model_phase2.A.deactivate()


### Bounds on days and shifts worked over the week

model_phase2.tt_min_dys_weeks = Param(model_phase2.TTYPES, model_phase2.WEEKS, default=0.0)       # Minimum number of days worked by week by tour type
model_phase2.tt_max_dys_weeks = Param(model_phase2.TTYPES, model_phase2.WEEKS, default=1e+6)         # Maximum number of days worked by week by tour type

model_phase2.tt_min_cumul_dys_weeks = Param(model_phase2.TTYPES, model_phase2.WEEKS, default=0.0) # Minimum number of days worked by cumulative weeks by tour type
model_phase2.tt_max_cumul_dys_weeks = Param(model_phase2.TTYPES, model_phase2.WEEKS, default=1e+6) # Maximum number of days worked by cumulative weeks by tour type


# Weekends consisting of a Fri and Sat imply updated lower bounds on some of the daily tour type variables. 
# These were the key to modeling weekends worked patterns.

def FriSat_idx_rule(M):
#    return [(t,e,w) for t in M.activeTT for e in range(1,M.num_weekend_patterns[2,t] + 1) for w in model_phase2.WEEKS]
    index_list =[]
    for t in M.activeTT:
        numpats = M.num_weekend_patterns[2,t]
        for p in range(1,numpats+1):
            for w in M.WEEKS:
                index_list.append((t,p,w))
                    
    return index_list

model_phase2.FriSat_idx = Set(dimen=3,ordered=True,initialize=FriSat_idx_rule)

model_phase2.FriSat_min_dys_weeks = Param(model_phase2.FriSat_idx, default=0)
model_phase2.FriSat_min_cumul_dys_weeks = Param(model_phase2.FriSat_idx, default=0)



model_phase2.tt_shiftlen_min_dys_weeks = Param(model_phase2.TTYPES, model_phase2.LENGTHS, model_phase2.WEEKS, default=0)        # Minimum number of days worked by week by shiftlen by tour type
model_phase2.tt_shiftlen_max_dys_weeks = Param(model_phase2.TTYPES, model_phase2.LENGTHS, model_phase2.WEEKS, default=1e+6)         # Maximum number of days worked by week by shiftlen by tour type

model_phase2.tt_shiftlen_min_cumul_dys_weeks = Param(model_phase2.TTYPES, model_phase2.LENGTHS, model_phase2.WEEKS, default=0)         # Minimum number of days worked by cumulative weeks by shiftlen by tour type
model_phase2.tt_shiftlen_max_cumul_dys_weeks = Param(model_phase2.TTYPES, model_phase2.LENGTHS, model_phase2.WEEKS, default=1e+6)         # Maximum number of days worked by cumulative weeks by shiftlen by tour type

model_phase2.tt_min_prds_weeks = Param(model_phase2.TTYPES, model_phase2.WEEKS, default=0)        # Minimum number of periods worked by week by tour type
model_phase2.tt_max_prds_weeks = Param(model_phase2.TTYPES, model_phase2.WEEKS, default=1e+6)         # Maximum number of periods worked by week by tour type

model_phase2.tt_min_cumul_prds_weeks = Param(model_phase2.TTYPES, model_phase2.WEEKS, default=0)         # Minimum number of periods worked by cumulative weeks by tour type
model_phase2.tt_max_cumul_prds_weeks = Param(model_phase2.TTYPES, model_phase2.WEEKS, default=1e+6)         # Minimum number of periods worked by cumulative weeks by tour type

model_phase2.tt_shiftlen_min_prds_weeks = Param(model_phase2.TTYPES, model_phase2.LENGTHS, model_phase2.WEEKS, default=0)        # Minimum number of periods worked by week by tour type by shift length
model_phase2.tt_shiftlen_max_prds_weeks = Param(model_phase2.TTYPES, model_phase2.LENGTHS, model_phase2.WEEKS, default=1e+6)       # Minimum number of periods worked by week by tour type by shift length

model_phase2.tt_shiftlen_min_cumul_prds_weeks = Param(model_phase2.TTYPES, model_phase2.LENGTHS, model_phase2.WEEKS, default=0)         # Minimum number of periods worked by cumulative weeks by tour type by shift length
model_phase2.tt_shiftlen_max_cumul_prds_weeks = Param(model_phase2.TTYPES, model_phase2.LENGTHS, model_phase2.WEEKS, default=1e+6)        # Minimum number of periods worked by cumulative weeks by tour type by shift length


### Bounds on tour type variables
model_phase2.tt_lb =  Param(model_phase2.TTYPES)       # RHS from .MIX
model_phase2.tt_ub =  Param(model_phase2.TTYPES,  default=infinity)

def activeTT_init(M):
    return [t for t in M.TTYPES if M.tt_ub[t] > 0]
    
model_phase2.activeTT = Set(dimen=1,ordered=True,initialize=activeTT_init)

### To the above, we'll add params and sets to allow direct modeling of side constraints
### of the form sum{subset of tour types} =, >=, <= some bound

model_phase2.tt_parttime = Param(model_phase2.TTYPES)     # 1 for part-time, 0 for full-time


# Allowable shift start times  - note that these are tour type specific
model_phase2.allow_start = Param(model_phase2.PERIODS,model_phase2.DAYS,model_phase2.LENGTHS,model_phase2.TTYPES,default=0.0)

def okShifts_rule(M):
    return [(i,j,w,k,t) for i in M.PERIODS
                        for j in M.DAYS
                        for w in M.WEEKS
                        for k in M.LENGTHS
                        for t in M.activeTT
                        if M.allow_start[i,j,k,t] > 0 and k in M.tt_length_x[t]]


model_phase2.okShifts = Set(dimen=5,ordered=True,initialize=okShifts_rule)


def okShiftTypes_rule(M):
    return [(i,j,k,t) for i in M.PERIODS
                        for j in M.DAYS
                        for k in M.LENGTHS
                        for t in M.activeTT                        
                        if M.allow_start[i,j,k,t] > 0 and k in M.tt_length_x[t]]


model_phase2.okShiftTypes = Set(dimen=4,ordered=True,initialize=okShiftTypes_rule)


# Limits on part time labor and limits on total labor

model_phase2.max_parttime_frac = Param()                   # Maximum fraction of labor hours covered by part-time employees
model_phase2.labor_budget = Param()                   # Maximum labor expenditure

# -----------------------------------------------------------------------
# COST RELATED PARAMETERS
# -----------------------------------------------------------------------

model_phase2.tt_cost_multiplier = Param(model_phase2.TTYPES)           # Tour type differential

model_phase2.cu1 = Param()
model_phase2.cu2 = Param()        
model_phase2.usb = Param()

# -----------------------------------------------------------------------
# Weekend RELATED PARAMETERS
# -----------------------------------------------------------------------

model_phase2.midnight_thresh = Param(model_phase2.TTYPES, default=1e+6) # Need to generalize for any number of periods per day

def weekend_init(M,i,t):

    result = []
    lens = [M.lengths[k] for k in M.tt_length_x[t]]
    maxlen = max(lens)
    if i+maxlen-1>=M.midnight_thresh[t]:
        result.append(6)
    else:
        result.append(1)
    result.append(7)
    return result

model_phase2.weekend = Set(model_phase2.WINDOWS, model_phase2.TTYPES, ordered=True, initialize=weekend_init)
 
def weekend_type_init(M,i,t):

    result = 1
    lens = [M.lengths[k] for k in M.tt_length_x[t]]
    maxlen = max(lens)
    if i+maxlen-1>=M.midnight_thresh[t]:
        result = 2
    
    return result
        
          
model_phase2.weekend_type = Param(model_phase2.WINDOWS, model_phase2.TTYPES, initialize=weekend_type_init)
  
# -----------------------------------------------------------------------
# Coverage related PARAMETERS
# -----------------------------------------------------------------------

# Target and minimum staffing levels - this is week specific. We can always allow user to input
# a single week and then repeat it for the other weeks.

model_phase2.dmd_staff = Param(model_phase2.PERIODS,model_phase2.DAYS,model_phase2.WEEKS)
model_phase2.min_staff = Param(model_phase2.PERIODS,model_phase2.DAYS,model_phase2.WEEKS)

# -----------------------------------------------------------------------
# START WINDOWS - should these be tour type specific since allow start is tour type specific?
#
# To start with,I made them week specific but not sure if they should be. Seems they should be 
# period, day, tour type.
#
# Maybe best strategy is to do the variables and constraints first and then work backwards to
# define windows appropriately.
# -----------------------------------------------------------------------

model_phase2.g_start_window_width = Param()                            # Width of start-time windows



##/**** Beginning of each start window (in total periods from Sunday @ midnight)****/
##param b_window_wepoch{i in PERIODS,j in DAYS} := n_prds_per_day*(j-1)+i;
##
##/**** End of each start window (in total periods from Sunday @ midnight) ****/
##param e_window_wepoch{i in PERIODS,j in DAYS} :=
## ( if n_prds_per_day*(j-1)+i+width <= n_prds_per_day*n_days then
##    n_prds_per_day*(j-1)+i+width
##  else
##    (n_prds_per_day*(j-1)+i+width )-n_prds_per_day*n_days);



# This version spans multiple weeks
def b_window_epoch_init(M, i, j, w):
    return M.n_days_per_week*M.n_prds_per_day*(w-1) + M.n_prds_per_day*(j-1)+i

model_phase2.b_window_epoch = Param(model_phase2.bins,initialize=b_window_epoch_init)


def e_window_epoch_init(M,i,j,w):
    epoch = M.n_days_per_week * M.n_prds_per_day*(w-1) + M.n_prds_per_day * (j-1)+i+M.g_start_window_width
    if epoch <= M.n_prds_per_cycle:
        e_window = epoch
    else:
        e_window = epoch-M.n_prds_per_cycle
    
    return e_window

model_phase2.e_window_epoch = Param(model_phase2.bins,initialize=e_window_epoch_init)

###/**** The set WindowWepochs{(i,j) in {PERIODS,DAYS}} contains all pairs 
###   (l=period,m=day) within the start window which begins in period (i,j). ****/
###
###set WindowWepochs{i in PERIODS,j in DAYS}
###       within {PERIODS,DAYS} := {(l,m) in {PERIODS,DAYS}:
###(  (  (n_prds_per_day*(m-1)+l>=b_window_wepoch[i,j]) and
###      (n_prds_per_day*(m-1)+l<=
###    (if b_window_wepoch[i,j]<=e_window_wepoch[i,j] then e_window_wepoch[i,j]
###     else n_prds_per_day*n_days)
###      )
###   ) or
###   ( (n_prds_per_day*(m-1)+l>=
###       (if b_window_wepoch[i,j]<=e_window_wepoch[i,j] then b_window_wepoch[i,j]
###    else 1)
###     ) and (n_prds_per_day*(m-1)+l<=e_window_wepoch[i,j])
###   )
###)       };

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
   
model_phase2.PotentialGlobalStartWindow = Set(model_phase2.PERIODS,model_phase2.DAYS,model_phase2.WEEKS,dimen=3,ordered=True,initialize=PotentialGlobalStartWindow_init)




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
# Defn:  PotentialStartWindow[i,j,w,k,t] contains all the bin triplets within g_start_window_width periods of
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


model_phase2.PotentialStartWindow_idx = Set(dimen=5,ordered=True,initialize=PotentialStartWindow_idx_rule)

def PotentialStartWindow_init(M,i,j,w,k,t):
    return [(l,m,n) for (l,m,n) in M.PotentialGlobalStartWindow[i,j,w]
                  if M.allow_start[i,j,k,t] > 0 and M.allow_start[l,m,k,t] > 0]
    
model_phase2.PotentialStartWindow = Set(model_phase2.PotentialStartWindow_idx,ordered=True,initialize=PotentialStartWindow_init,dimen=3)
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


model_phase2.okStartWindowRoots_idx = Set(dimen=2,ordered=True,initialize=okStartWindowRoots_idx_rule)
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
    

model_phase2.okStartWindowRoots = Set(model_phase2.okStartWindowRoots_idx,dimen=3,ordered=True,initialize=okStartWindowRoots_init)

##
##
### CHAINS - NON-OVERLAPPING SEQUENCES OF (i,j) period,day PAIRS THAT
###          CAN BE ISOLATED FOR COORDINATING ix AND DWT VARIABLES.
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
    for (i,j,w) in M.okStartWindowRoots[t,k]:
        for (p,q,r) in (M.okStartWindowRoots[t,k] - Set(initialize=[(i,j,w)])):
            if (i,j,w) not in M.PotentialStartWindow[p,q,r,k,t]:
                window_list.append((i,j,w))
         
    return window_list
#    
#
model_phase2.bchain = Set(model_phase2.okStartWindowRoots_idx,dimen=3,ordered=True,initialize=bchain_init)
#
##set echain {t in okTTYPES,k in LENGTHS,i in PERIODS,j in DAYS:
## (i,j) in bchain[t,k]} 
## := setof{(w,x) in {PERIODS,DAYS}: (w,x) not in (bchain[t,k] diff {(i,j)}) and
## (w,x) in okWindowBeginnings[t,k] and
## forall{(p,q) in (okWindowBeginnings[t,k] diff {(w,x)})} (p,q) not in (okWindowWepochs[w,x,k,t]) and
##     (period[w,x]>=period[i,j] and 
##     forall{(n,o) in bchain[t,k]: period[n,o]>period[i,j]} period[w,x]<
##  period[n,o] or 
##  (period[w,x]<period[i,j] and 
##  forall{(n,o) in bchain[t,k]} (period[w,x]< 
##   period[n,o] and period[n,o]<=period[i,j]) ) )
##  } (w,x) ;
#
#
def chain_idx_rule(M):
    return [(t,k,i,j,w) for t in M.activeTT
                      for k in M.LENGTHS
                      for i in M.PERIODS
                      for j in M.DAYS
                      for w in M.WEEKS
                      if (t,k) in M.okStartWindowRoots_idx and (i,j,w) in M.bchain[t,k]]


model_phase2.chain_idx = Set(dimen=5, ordered=True, initialize=chain_idx_rule)
#
def echain_init(M,t,k,i,j,w):
    window_list =[]
    # Compute global period of (i,j,w)
    g_prd = M.g_period[i, j, w]
    prd = g_prd_to_tuple(M, g_prd)
    done = False
    
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


model_phase2.echain = Set(model_phase2.chain_idx,ordered=True,dimen=3,initialize=echain_init)


def n_links_init(M,t,k,i,j,w):
    
    # Compute global period of (i,j,w)
    b_g_prd = M.g_period[i, j, w]
    e_prd = M.echain[t, k, i, j, w][1]
    e_g_prd = M.g_period[e_prd[0], e_prd[1], e_prd[2]]
        
    return g_period_difference(M,b_g_prd,e_g_prd)


model_phase2.n_links = Param(model_phase2.chain_idx,initialize=n_links_init)


def chain_init(M,t,k,i,j,w):
    window_list =[(i,j,w)]
    # Compute global period of (i,j,w)
    g_prd = M.g_period[i, j, w]
    done = False
    
    steps = 1
    while steps <= M.n_links[t,k,i,j,w]:
        g_prd_next = g_period_increment(M,g_prd,steps)
        prd_next = g_prd_to_tuple(M, g_prd_next)
        if prd_next in M.okStartWindowRoots[t,k]:
            window_list.append(prd_next)
        steps = steps + 1
                
        
    return window_list


model_phase2.chain = Set(model_phase2.chain_idx,ordered=True,dimen=3,initialize=chain_init)


def link_idx_rule(M):
      
    return [(t,k,i,j,w,m) for t in M.activeTT
                          for k in M.LENGTHS
                          for i in M.PERIODS
                          for j in M.DAYS
                          for w in M.WEEKS
                          for m in range(1,M.n_prds_per_cycle)
                          if (t,k) in M.okStartWindowRoots_idx and (i,j,w) in M.bchain[t,k]
                              and m <= M.n_links[t,k,i,j,w]]
                          


model_phase2.link_idx = Set(dimen=6, ordered=True, initialize=link_idx_rule)

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


model_phase2.link = Set(model_phase2.link_idx,ordered=True,dimen=3,initialize=link_init)

#set linkspan{t in okTTYPES, k in tt_length_x[t], (i,j) in bchain[t,k], m in 1..numlinks[t,k,i,j] } dimen 2 := 
#    {(n,o) in {PERIODS,DAYS}: forall{(p,d) in link[t,k,i,j,m]} (n,o) in WindowWepochs[p,d] union
#        if m=1 then WindowWepochs[p,d] else linkspan[t,k,i,j,m-1]};

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


model_phase2.linkspan = Set(model_phase2.link_idx,dimen=3,ordered=True,initialize=linkspan_init)


##### Phase 1 Shift variables --> Phase 2 parameters



model_phase2.Shift = Param(model_phase2.okShifts, default=0)

# Shift[i,j,w,k,t] = Number of shifts of length k starting in period i
# of day j in week w for a tour of type t


model_phase2.okTourType = model_phase2.WINDOWS * model_phase2.activeTT

def TourType_idx_rule(M):
    return [(i,t) for i in M.WINDOWS                        
                      for t in M.activeTT
                      if (i,t) in M.okTourType]
                        


model_phase2.TourType_idx = Set(dimen=2,initialize=TourType_idx_rule)



model_phase2.TourType = Param(model_phase2.TourType_idx, default=0)


    
##### Daily tour type variables
model_phase2.okDailyTourType = model_phase2.WINDOWS * model_phase2.TTYPES * model_phase2.DAYS


#    /* DailyTourType[i,t,d] Number of employees working tour type t
#       starting in window i and working day d in week w*/

def DailyTourType_idx_rule(M):
    return [(i,t,j,w) for i in M.WINDOWS                        
                      for t in M.activeTT
                      for j in M.DAYS
                      for w in M.WEEKS
                      if (i,t,j) in M.okDailyTourType]
                        


model_phase2.DailyTourType_idx = Set(dimen=4,initialize=DailyTourType_idx_rule)



model_phase2.DailyTourType = Param(model_phase2.DailyTourType_idx)



##### Daily shift worked variables

#var DailyShiftWorked{i in 1..n_windows,t in okTTYPES,k in tt_length_x[t],d in DAYS,w in WEEKS : (i,t,d) in okDailyTourType}
#    >= 0, integer;
    
# TODO - modify this index when get width>0 working. These are windows, not period start times.
# Just using okShifts for now for w=0 case.  
def DailyShiftWorked_idx_rule(M):
    return [(i,t,k,j,w) for i in M.WINDOWS                        
                        for t in M.activeTT
                        for k in M.tt_length_x[t]
                        for j in M.DAYS
                        for w in M.WEEKS
                        if (i,t,j) in M.okDailyTourType]
                        


model_phase2.DailyShiftWorked_idx = Set(dimen=5,initialize=DailyShiftWorked_idx_rule)


model_phase2.DailyShiftWorked = Param(model_phase2.DailyShiftWorked_idx, default=0)

##### Weekend Days off variables   

#var WeekendDaysWorked{p in 1..max_weekend_patterns,i in 1..n_windows,t in okTTYPES : 
#                      (i,t) in okTourType and p <= num_weekend_patterns[weekend_type[i,t],t]} 
#       >= 0, integer ;
       
def ok_weekenddaysworked_idx_rule(M):
    index_list =[]
    for (i,t) in M.okTourType:
        for pattern in sequence(M.max_weekend_patterns):
            weekendtype = M.weekend_type[i,t]
            if pattern <= M.num_weekend_patterns[weekendtype,t]:
                index_list.append((i,t,pattern))
    
    return index_list
                        
model_phase2.ok_weekenddaysworked_idx = Set(dimen=3,initialize=ok_weekenddaysworked_idx_rule)



model_phase2.WeekendDaysWorked = Param(model_phase2.ok_weekenddaysworked_idx, default=0)  
    
    # WeekendDaysWorked[d,i,t] = Number of employees working days-off patterns d
    # in start window i and of tour type t
    
    
#param n_tours:=sum{i in 1..n_windows,t in okTTYPES :(i,t) in okTourType } TourType[i,t];

model_phase2.n_tours = Param(within=PositiveIntegers)
model_phase2.TOURS = RangeSet(1,model_phase2.n_tours)

###set index_tours{i in 1..n_windows,t in okTTYPES:(i,t) in okTourType }  :=
### { (sum{k in 1..i, j in 1..(if k<i then n_tts else t):
###     j in okTTYPES and (k,j) in okTourType}(TourType[k,j])-TourType[i,t]+1)..
###   (sum{k in 1..i, j in 1..(if k<i then n_tts else t):
###     j in okTTYPES and (k,j) in okTourType}TourType[k,j]) };
###set WIN_sx{l in 1..n_tours}:={i in 1..n_windows:
###    exists{t in okTTYPES} (i,t) in okTourType
###      and l in index_tours[i,t]};
###
###
###set TT_sx{l in 1..n_tours}:={t in okTTYPES:
###    exists{(i,t) in okTourType}
###      l in index_tours[i,t]};
###
###param WIN_x{l in 1..n_tours}:=min{i in WIN_sx[l]} i;
###param TT_x{l in 1..n_tours}:=min{i in TT_sx[l]} i;



model_phase2.WIN_x = Param(model_phase2.TOURS)
model_phase2.TT_x = Param(model_phase2.TOURS)




#......................................................Dec. Variables


###var tourshift{s in 1..n_tours, i in PERIODS, j in DAYS, w in WEEKS, k in LENGTHS, t in TTYPES, p in PERIODS, d in DAYS :
###       t=TT_x[s] and p=WIN_x[s] and k in tt_length_x[t] and (i,j,k,t) in ok_shifts and 
###       (i,j) in okWindowWepochs[p,d,k,t]} binary ;
###    
###    # tourshift[s,i,j,w,k,t,p,d,q] = 1 if an x[i,j,w,k,t] shift is assigned to tour s in window (p,d,q)
###    #            = 0 otherwise
### 
###

def TourShifts_idx_rule(M):
    index_list = []
    for s in M.TOURS:
        for i in M.PERIODS: 
            for j in M.DAYS:                
                for w in M.WEEKS:
                    for t in [a for a in M.activeTT if a == M.TT_x[s]]:
                        for k in M.tt_length_x[t]:
                            for (p,d,q) in [(x,y,z) for (x,y,z) in M.okStartWindowRoots[t,k] if x == M.WIN_x[s]]:
                            #for p in M.WINDOWS:
                            #    for d in M.DAYS:
                                if (i,j,w) in M.PotentialStartWindow[p,d,q,k,t]:
                                    index_list.append((s,i,j,w,k,t))
                        
    return index_list

model_phase2.TourShifts_idx = Set(dimen=6,initialize=TourShifts_idx_rule)

model_phase2.TourShifts = Var(model_phase2.TourShifts_idx, within=Boolean)




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
        for i in [a for a in M.WINDOWS if a == M.WIN_x[s]]:
            for t in [a for a in M.activeTT if a == M.TT_x[s]]: 
                for pattern in sequence(M.num_weekend_patterns[M.weekend_type[i,t],t]):
                    if (i,t) in M.okTourType:
                        index_list.append((s,pattern,i,t))
                                       
    return index_list
                
model_phase2.TourWkendDof_idx = Set(dimen=4,initialize=TourWkendDof_idx_rule)

model_phase2.TourWkendDof = Var(model_phase2.TourWkendDof_idx, within=Boolean)

#.......................................................Obj. Function

### Objective function is just the sum of the tourshift variables - which makes it a search for a feasible solution.
###
###minimize total_cost :
###        
### sum {s in 1..n_tours,i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS, t in TTYPES, p in PERIODS, d in DAYS :
###       t=TT_x[s] and p=WIN_x[s] and k in tt_length_x[t] and (i,j,k,t) in ok_shifts and 
###       (i,j) in okWindowWepochs[p,d,k,t]} tourshift[s,i,j,w,k,t,p,d] ; 
 
# Objective function
def objective_rule(M):
    obj1 = sum(M.TourShifts[s,i,j,w,k,t] for (s,i,j,w,k,t) in M.TourShifts_idx)
    return obj1

model_phase2.total_num_tours = Objective(rule=objective_rule, sense=minimize)

# One days off pattern for each tour 

#subject to one_tourdof{s in 1..n_tours} :
#    sum{d in 1..num_weekend_patterns[weekend_type[WIN_x[s],TT_x[s]],TT_x[s]]} tourdof[s,d,WIN_x[s],TT_x[s]]=1;


def OnePatternPerTourShift_rule(M,s):
    
    return sum(M.TourWkendDof[s,pattern,M.WIN_x[s],M.TT_x[s]] 
        for pattern in sequence(M.num_weekend_patterns[M.weekend_type[M.WIN_x[s],M.TT_x[s]],M.TT_x[s]])) == 1

model_phase2.OnePatternPerTourShift = Constraint(model_phase2.TOURS,rule=OnePatternPerTourShift_rule)

# All days off patterns assigned 

#subject to conserve_ID {d in 1..max_weekend_patterns,i in WINDOWS,t in okTTYPES: d <= num_weekend_patterns[weekend_type[i,t],t] and (i,t) in okTourType} :
# sum{s in 1..n_tours : t=TT_x[s] and i=WIN_x[s]} 
#        tourdof[s,d,i,t]= WeekendDaysWorked[d,i,t];

# Sum over the tours within each (i,t) and make sure they add up to the WeekendDaysWorked variables
def Tour_WkendDof_conservation_idx_rule(M):
    index_list = []
    for i in M.WINDOWS:
        for t in M.activeTT: 
            for pattern in sequence(M.max_weekend_patterns):
                if (i,t) in M.okTourType and pattern <= M.num_weekend_patterns[M.weekend_type[i,t],t]:
                    index_list.append((pattern,i,t))
    return index_list   
    
model_phase2.Tour_WkendDof_conservation_idx = Set(dimen=3,initialize=Tour_WkendDof_conservation_idx_rule)

def Tour_WkendDof_conservation_rule(M,pattern,i,t):
    return sum(M.TourWkendDof[s,pattern,i,t] for s in M.TOURS if i == M.WIN_x[s] and t == M.TT_x[s]) == \
        M.WeekendDaysWorked[i,t,pattern]
                              
model_phase2.Tour_WkendDof_conservation = Constraint(model_phase2.Tour_WkendDof_conservation_idx,rule=Tour_WkendDof_conservation_rule)





# Weekend days worked for each tour determines number of weekend shifts assigned to each tour

###subject to tourshift_tourdof_integration1{s in 1..n_tours,i in WINDOWS,j in DAYS, w in WEEKS,t in TTYPES: 
###    t=TT_x[s] and i=WIN_x[s] and j in weekend[i,TT_x[s]]} : 
       
###     sum{k in tt_length_x[t],(p,d) in okWindowWepochs[i,j,k,t] } tourshift[s,p,d,w,k,t,i,j] =
###        (sum{g in 1..num_weekend_patterns[weekend_type[i,t],t]}tourdof[s,g,i,t]*A[g,j,w,t,weekend_type[i,t]]);

def Tour_ShiftWkendDof_integration1_idx_rule(M):
    index_list = []
    for s in M.TOURS:
        for i in [a for a in M.WINDOWS if a == M.WIN_x[s]]:
            for t in [a for a in M.activeTT if a == M.TT_x[s]]: 
                for d in M.weekend[i,t]:
                    for w in M.WEEKS: 
                        index_list.append((s,i,d,w,t))
    return index_list
    
    
model_phase2.Tour_ShiftWkendDof_integration1_idx = Set(dimen=5,initialize=Tour_ShiftWkendDof_integration1_idx_rule)

def weekend_integration_1_SS_idx_rule(M):
    return [(j,w,i,t) for j in M.DAYS
                    for w in M.WEEKS
                    for i in M.WINDOWS                        
                    for t in M.activeTT
                    if j in M.weekend[i,t] and 1 in M.weekend[i,t] and 7 in M.weekend[i,t] \
                      and (i,t,j) in M.okDailyTourType]
                      
model_phase2.weekend_integration_1_SS_idx = Set(dimen=4,initialize=weekend_integration_1_SS_idx_rule) 
                     
def action_check_WeekendDaysWorked_DailyTourType_rule(M,j,w,i,t):
    val_w = value(sum(M.A[p,j,w,t,1] * M.WeekendDaysWorked[i,t,p] for p in sequence(M.num_weekend_patterns[1,t])))
    val_d = M.DailyTourType[i,t,j,w]
    
    if val_w == val_d:
        check = 'OK'
    else:
        check = 'BAD'

    msg = '({},{},{},{}) w={} d={} --> {}'.format(i,j,w,t,val_w,val_d,check)
    if check == 'BAD':
        print msg
    
model_phase2.action_check_WeekendDaysWorked_DailyTourType = \
    BuildAction(model_phase2.weekend_integration_1_SS_idx,rule=action_check_WeekendDaysWorked_DailyTourType_rule)

def Tour_ShiftWkendDof_integration1_rule(M,s,i,j,w,t):
    return sum(M.TourShifts[s,p,d,q,k,t] for k in M.tt_length_x[t] for (p,d,q) in M.PotentialStartWindow[i,j,w,k,t]) == \
                              sum(M.A[pattern,j,w,t,M.weekend_type[i,t]] * M.TourWkendDof[s,pattern,i,t] for pattern in sequence(M.num_weekend_patterns[M.weekend_type[i,t],t]))
                 
model_phase2.Tour_ShiftWkendDof_integration1 = Constraint(model_phase2.Tour_ShiftWkendDof_integration1_idx,rule=Tour_ShiftWkendDof_integration1_rule)






# No more than one shift worked per day

#subject to tours_daily{s in 1..n_tours,j in DAYS,w in WEEKS} :
#  sum{k in tt_length_x[TT_x[s]],(p,d) in okWindowWepochs[WIN_x[s],j,k,TT_x[s]]}
#   tourshift[s,p,d,w,k,TT_x[s],WIN_x[s],j] <= 1;


def Tours_Daily_idx_rule(M):
    return [(s,j,w) for s in M.TOURS                        
                    for j in M.DAYS
                    for w in M.WEEKS]

    
    
model_phase2.Tours_Daily_idx = Set(dimen=3,initialize=Tours_Daily_idx_rule)

def Tours_Daily_rule(M,s,j,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]] 
                 for k in M.tt_length_x[M.TT_x[s]] for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,w,k,M.TT_x[s]]) <= 1
                                    
model_phase2.Tours_Daily = Constraint(model_phase2.Tours_Daily_idx,rule=Tours_Daily_rule)


## Needs to be generalized for start windows
#subject to tours_daily_conservation2{i in PERIODS, j in DAYS,w in WEEKS, k in LENGTHS, t in TTYPES: k in tt_length_x[t] and (i,j,k,t) in ok_shifts} :
#  sum{s in 1..n_tours: TT_x[s]=t and WIN_x[s]=i }
#   tourshift[s,i,j,w,k,t,WIN_x[s],j] = Shift[i,j,w,k,t];  
# 
#
#
## All the shifts need to get assigned to the tours
#
#/*
#subject to tours_daily_conservation2{p in PERIODS, d in DAYS,w in WEEKS, k in LENGTHS, t in TTYPES: k in tt_length_x[t]} :
#  sum{s in 1..n_tours, (i,j) in okWindowWepochs[p,d,k,t]: TT_x[s]=t and WIN_x[s]=p }
#   tourshift[s,i,j,w,k,t,p,d] = sum{(i,j) in okWindowWepochs[p,d,k,t] : (i,j,k,t) in ok_shifts}Shift[i,j,w,k,t]; 
#*/

def Tours_Daily_conservation_idx_rule(M):
    index_list = []
    for p in M.PERIODS: 
        for d in M.DAYS:                
            for w in M.WEEKS:
                for t in M.activeTT:
                    for k in M.tt_length_x[t]:
                        if (p,d,w) in M.okStartWindowRoots[t,k]:
                            index_list.append((p,d,w,k,t))
    return index_list

model_phase2.Tours_Daily_conservation_idx = Set(dimen=5,initialize=Tours_Daily_conservation_idx_rule)

def Tours_Daily_conservation_rule(M,p,d,q,k,t):
    
    return sum(M.TourShifts[s,i,j,w,k,t] 
                 for s in M.TOURS for (i,j,w) in M.PotentialStartWindow[p,d,q,k,t] if p == M.WIN_x[s] and t == M.TT_x[s]) \
                 == sum(M.Shift[i,j,w,k,t] for (i,j,w) in M.PotentialStartWindow[p,d,q,k,t])
                                    
model_phase2.Tours_Daily_conservation = Constraint(model_phase2.Tours_Daily_conservation_idx,rule=Tours_Daily_conservation_rule)

#=============== Tours_Weekly Lower and Upper Bounds on Days Worked ==================================================

#== For each (tour type) make sure the number of shifts assigned each week
#== satisfies associated lower and upper bounds.

## Tour type specific bounds on number of days worked over the weeks (both cumulative and non-cumulative)
#subject to tours_weekly_LB{t in 1..n_tours,w in WEEKS} :
#  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_min_dys_weeks[TT_x[t],w];
#   
#subject to tours_weekly_UB{t in 1..n_tours,w in WEEKS} :
#  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_max_dys_weeks[TT_x[t],w];

def Tours_Weekly_LB_idx_rule(M):
    return [(s,w) for s in M.TOURS for w in M.WEEKS]

model_phase2.Tours_Weekly_LB_idx = Set(dimen=2,initialize=Tours_Weekly_LB_idx_rule)

def Tours_Weekly_LB_rule(M,s,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]] 
                 for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,w,k,M.TT_x[s]]) >= M.tt_min_dys_weeks[M.TT_x[s],w]

                                    
model_phase2.Tours_Weekly_LB = Constraint(model_phase2.Tours_Weekly_LB_idx,rule=Tours_Weekly_LB_rule)

def Tours_Weekly_UB_idx_rule(M):
    return [(s,w) for s in M.TOURS for w in M.WEEKS]

model_phase2.Tours_Weekly_UB_idx = Set(dimen=2,initialize=Tours_Weekly_UB_idx_rule)

def Tours_Weekly_UB_rule(M,s,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]] 
                 for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,w,k,M.TT_x[s]]) <= M.tt_max_dys_weeks[M.TT_x[s],w]

                                    
model_phase2.Tours_Weekly_UB = Constraint(model_phase2.Tours_Weekly_UB_idx,rule=Tours_Weekly_UB_rule)

#=============== Tours_Total Lower and Upper Bounds on Days Worked ======================================================

#== For each (tour type) make sure the number of shifts assigned over cumulative
#== weeks satisfies associated lower and upper bounds.

 #Seems like the following two constraints need to be generalized for > 2 weeks.
#
#subject to tours_tot_LB{t in 1..n_tours} :
#  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_min_cumul_dys_weeks[TT_x[t],n_weeks];
#   
#subject to tours_tot_UB{t in 1..n_tours} :
#  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_max_cumul_dys_weeks[TT_x[t],n_weeks];


def Tours_Total_LB_idx_rule(M):
    return [(s,w) for s in M.TOURS for w in M.WEEKS]

model_phase2.Tours_Total_LB_idx = Set(dimen=2,initialize=Tours_Total_LB_idx_rule)

def Tours_Total_LB_rule(M,s,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]] 
                 for W in sequence(w) for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,W,k,M.TT_x[s]]) >= M.tt_min_cumul_dys_weeks[M.TT_x[s],w]

                                    
model_phase2.Tours_Total_LB = Constraint(model_phase2.Tours_Total_LB_idx,rule=Tours_Total_LB_rule)


def Tours_Total_UB_idx_rule(M):
    return [(s,w) for s in M.TOURS for w in M.WEEKS]

model_phase2.Tours_Total_UB_idx = Set(dimen=2,initialize=Tours_Total_UB_idx_rule)

def Tours_Total_UB_rule(M,s,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]] 
                 for W in sequence(w) for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,W,k,M.TT_x[s]]) <= M.tt_max_cumul_dys_weeks[M.TT_x[s],w]

                                    
model_phase2.Tours_Total_UB = Constraint(model_phase2.Tours_Total_UB_idx,rule=Tours_Total_UB_rule)


#=============== Tours_Shiftlen_Weekly Lower and Upper Bounds on Days Worked ================================================

#== For each (tour type, shift length) make sure the number of shifts assigned each week
#== satisfies associated lower and upper bounds.

# Tour type and shift length specific bounds on number of days worked over the weeks (both cumulative and non-cumulative)
#subject to tours_shiftlen_weekly_LB{t in 1..n_tours,k in tt_length_x[TT_x[t]],w in WEEKS} :
#  sum{d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_shiftlen_min_dys_weeks[TT_x[t],k,w];
#   
#subject to tours_shiftlen_weekly_UB{t in 1..n_tours,k in tt_length_x[TT_x[t]],w in WEEKS} :
#  sum{d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_shiftlen_max_dys_weeks[TT_x[t],k,w];
#   


def Tours_Shiftlen_Weekly_LB_idx_rule(M):
    return [(s,k,w) for s in M.TOURS for k in M.tt_length_x[M.TT_x[s]] for w in M.WEEKS]

model_phase2.Tours_Shiftlen_Weekly_LB_idx = Set(dimen=3,initialize=Tours_Shiftlen_Weekly_LB_idx_rule)

def Tours_Shiftlen_Weekly_LB_rule(M,s,k,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]] 
                 for j in M.DAYS for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,w,k,M.TT_x[s]]) >= M.tt_shiftlen_min_dys_weeks[M.TT_x[s],k,w]

                                    
model_phase2.Tours_Shiftlen_Weekly_LB = Constraint(model_phase2.Tours_Shiftlen_Weekly_LB_idx,rule=Tours_Shiftlen_Weekly_LB_rule)

def Tours_Shiftlen_Weekly_UB_idx_rule(M):
    return [(s,k,w) for s in M.TOURS for k in M.tt_length_x[M.TT_x[s]] for w in M.WEEKS]

model_phase2.Tours_Shiftlen_Weekly_UB_idx = Set(dimen=3,initialize=Tours_Shiftlen_Weekly_UB_idx_rule)

def Tours_Shiftlen_Weekly_UB_rule(M,s,k,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]] 
                 for j in M.DAYS for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,w,k,M.TT_x[s]]) <= M.tt_shiftlen_max_dys_weeks[M.TT_x[s],k,w]

                                    
model_phase2.Tours_Shiftlen_Weekly_UB = Constraint(model_phase2.Tours_Shiftlen_Weekly_UB_idx,rule=Tours_Shiftlen_Weekly_UB_rule)


#=============== Tours_Shiftlen_Total Lower and Upper Bounds on Days Worked ==============================================

#== For each (tour type, shift length) make sure the number of shifts assigned over cumulative
#== weeks satisfies associated lower and upper bounds.

#subject to tours_shiftlen_tot_LB{t in 1..n_tours,k in tt_length_x[TT_x[t]]} :
#  sum{w in WEEKS,d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_shiftlen_min_cumul_dys_weeks[TT_x[t],k,n_weeks];
#   
#subject to tours_shiftlen_tot_UB{t in 1..n_tours,k in tt_length_x[TT_x[t]]} :
#  sum{w in WEEKS,d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_shiftlen_max_cumul_dys_weeks[TT_x[t],k,n_weeks]; 


def Tours_Shiftlen_Total_LB_idx_rule(M):
    return [(s,k,w) for s in M.TOURS for k in M.tt_length_x[M.TT_x[s]] for w in M.WEEKS]

model_phase2.Tours_Shiftlen_Total_LB_idx = Set(dimen=3,initialize=Tours_Shiftlen_Total_LB_idx_rule)

def Tours_Shiftlen_Total_LB_rule(M,s,k,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]] 
                 for W in range(1,w+1) for j in M.DAYS for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,W,k,M.TT_x[s]]) >= M.tt_shiftlen_min_cumul_dys_weeks[M.TT_x[s],k,w]

model_phase2.Tours_Shiftlen_Total_LB = Constraint(model_phase2.Tours_Shiftlen_Total_LB_idx,rule=Tours_Shiftlen_Total_LB_rule)

def Tours_Shiftlen_Total_UB_idx_rule(M):
    return [(s,k,w) for s in M.TOURS for k in M.tt_length_x[M.TT_x[s]] for w in M.WEEKS]

model_phase2.Tours_Shiftlen_Total_UB_idx = Set(dimen=3,initialize=Tours_Shiftlen_Total_UB_idx_rule)

def Tours_Shiftlen_Total_UB_rule(M,s,k,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]] 
                 for W in range(1,w+1) for j in M.DAYS for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,W,k,M.TT_x[s]]) <= M.tt_shiftlen_max_cumul_dys_weeks[M.TT_x[s],k,w]

                                    
model_phase2.Tours_Shiftlen_Total_UB = Constraint(model_phase2.Tours_Shiftlen_Total_UB_idx,rule=Tours_Shiftlen_Total_UB_rule)


#=============== Tours_Weekly Lower and Upper Bounds on Periods Worked ===========================================

#== For each (tour type) make sure the number of periods assigned each week
#== satisfies associated lower and upper bounds.

# Tour type specific bounds on number of periods worked over the weeks (both cumulative and non-cumulative)
#subject to tours_weekly_prds_LB{t in 1..n_tours,w in WEEKS} :
#  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] >= tt_min_prds_weeks[TT_x[t],w];
#   
#subject to tours_weekly_prds_UB{t in 1..n_tours,w in WEEKS} :
#  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] <= tt_max_prds_weeks[TT_x[t],w];
#   

def Tours_Weekly_Prds_LB_idx_rule(M):
    return [(s,w) for s in M.TOURS for w in M.WEEKS]

model_phase2.Tours_Weekly_Prds_LB_idx = Set(dimen=2,initialize=Tours_Weekly_Prds_LB_idx_rule)

def Tours_Weekly_Prds_LB_rule(M,s,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]]*M.lengths[k] 
                 for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,w,k,M.TT_x[s]]) >= M.tt_min_prds_weeks[M.TT_x[s],w]

                                    
model_phase2.Tours_Weekly_Prds_LB = Constraint(model_phase2.Tours_Weekly_Prds_LB_idx,rule=Tours_Weekly_Prds_LB_rule)

def Tours_Weekly_Prds_UB_idx_rule(M):
    return [(s,w) for s in M.TOURS for w in M.WEEKS]

model_phase2.Tours_Weekly_Prds_UB_idx = Set(dimen=2,initialize=Tours_Weekly_Prds_UB_idx_rule)

def Tours_Weekly_Prds_UB_rule(M,s,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]]*M.lengths[k] 
                 for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,w,k,M.TT_x[s]]) <= M.tt_max_prds_weeks[M.TT_x[s],w]

model_phase2.Tours_Weekly_Prds_UB = Constraint(model_phase2.Tours_Weekly_Prds_UB_idx,rule=Tours_Weekly_Prds_UB_rule)


#=============== Tours_Total Lower and Upper Bounds on Periods Worked ===========================================

#== For each (tour type) make sure the number of periods assigned over cumulative weeks
#== satisfies associated lower and upper bounds.

#subject to tours_tot_prds_LB{t in 1..n_tours} :
#  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] >= tt_min_cumul_prds_weeks[TT_x[t],n_weeks];
#   
#subject to tours_tot_prds_UB{t in 1..n_tours} :
#  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
#   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] <= tt_max_cumul_prds_weeks[TT_x[t],n_weeks];  

def Tours_Total_Prds_LB_idx_rule(M):
    return [(s,w) for s in M.TOURS for w in M.WEEKS]

model_phase2.Tours_Total_Prds_LB_idx = Set(dimen=2,initialize=Tours_Total_Prds_LB_idx_rule)

def Tours_Total_Prds_LB_rule(M,s,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]]*M.lengths[k] 
                 for W in range(1,w+1) for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,W,k,M.TT_x[s]]) >= M.tt_min_cumul_prds_weeks[M.TT_x[s],w]

                                    
model_phase2.Tours_Total_Prds_LB = Constraint(model_phase2.Tours_Total_Prds_LB_idx,rule=Tours_Total_Prds_LB_rule)


def Tours_Total_Prds_UB_idx_rule(M):
    return [(s,w) for s in M.TOURS for w in M.WEEKS]

model_phase2.Tours_Total_Prds_UB_idx = Set(dimen=2,initialize=Tours_Total_Prds_UB_idx_rule)

def Tours_Total_Prds_UB_rule(M,s,w):
    
    return sum(M.TourShifts[s,p,d,q,k,M.TT_x[s]]*M.lengths[k] 
                 for W in range(1,w+1) for j in M.DAYS for k in M.tt_length_x[M.TT_x[s]] for (p,d,q) in M.PotentialStartWindow[M.WIN_x[s],j,W,k,M.TT_x[s]]) <= M.tt_max_cumul_prds_weeks[M.TT_x[s],w]

                                    
model_phase2.Tours_Total_Prds_UB = Constraint(model_phase2.Tours_Total_Prds_UB_idx,rule=Tours_Total_Prds_UB_rule)




 # Ensures that each day worked in tour l is assigned exactly one shift
# of an appropriate length and from an appropriate start window.


#subject to tot_conserve_shifts
#   :
#  sum {s in 1..n_tours,i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS, t in TTYPES, p in PERIODS, d in DAYS :
#       t=TT_x[s] and p=WIN_x[s] and k in tt_length_x[t] and (i,j,k,t) in ok_shifts and 
#       (i,j) in okWindowWepochs[p,d,k,t]} tourshift[s,i,j,w,k,t,p,d] =
#   sum{m in PERIODS,n in DAYS,z in WEEKS, o in LENGTHS, t in TTYPES : 
#   (m,n,o,t) in ok_shifts and Shift[m,n,z,o,t]>0} Shift[m,n,z,o,t];



#def Tours_Total_conservation_rule(M):
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
#model_phase2.Tours_Total_conservation = Constraint(model_phase2.Tours_Total_conservation_idx,rule=Tours_Total_conservation_rule)










#
#
def main():
    pass

if __name__ == '__main__':
    main()
