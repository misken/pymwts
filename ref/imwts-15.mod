#/***********************************************************************
#                                                                      
#  Program:  imwts.mod
#  Language: AMPL/GMPL                                                       
#  Author:   Mark Isken                                                 
#  Create Date:     06/01/2011
#                                                                      
#  imwts-13.mod generates implicit multi-week tour scheduling models.
#  This version is the second main approach attempted. imwts-01-04.mod was a family
#  of models that accomodated a variety of scheduling characteristics but was
#  not fully capable of handling weekends via implicit methods only.
#  In this new model, we attempt to use explicit weekends patterns.
#                                                                      
#  This version accomodates:
#
#  (1) Tour types consisting of one or more shift lengths and limits on number of days worked and
#       number of periods worked per week and over cumulative weeks.
#  (2) Allowable shift start times that are tour type specific.
#  (3) Two-tiered understaffing costs and variables
#  (4) Maximum fraction of part-time personnel
#  (5) Budget restrictions
#  (6) Multi week planning horizon
#  (7) Number of periods per day can be 1,2,4,6,8,12,24,48, or 96
#  (8) Coverage constraints using standard set covering constraints.
#  (9) Explicit days off patterns are modeled for weekends.  
# 
#  The objective is to minimize the sum of tour type costs
#  subject to constraints on minimum staff
#  size, desired coverage levels, maximum level of part-time personnel
#  and budget restrictions.
#
#  MathProg specific comments
#  --------------------------
#  This version of the model is modified to be compatible with  
#  MathProg (an open source AMPL parser) which is included 
#  in the GNU GLPK package. http://www.gnu.org/software/glpk/glpk.html
# 
#  (1) No ordered sets - old model used 'circular' keyword (unnecessarily)
#  (2) No iterated union operator in MathProg. Required rewriting the 
#      chains_tot and chain_sweep_l constraints.
#
# v04 differs from v03 in that weekends constraints are only applied to full
# time tour types.
#v05 uses weekends off patterns
#v06 uses weekends off pattern constraint specific lower bounds on DWT
#
#v11 (assuming that v10 is truly working) generalizes some things that I
# shortcutted in an attempt to get the basic MWITS model working in terms
# of handling weekends off constraints. 

# (1) Start windows have been modelled but not thoroughly tested
# (2) Weekends off patterns have been generalized for multiple weeks
# (3) Since weekends of type (Fri, Sat) might require some modification of 
#     FriSat_min_dys_weeks and FriSat_max_dys_weeks, so far I've just
#     tested with (Sat,Sun) weekend types.
# (4) Need to set up a whole bunch of test scenarios with many tour types.
# (5) No machinery to create input files
# (6) No machinery for visualizing/checking solution schedules
# (7) No post processing to deal with max workstretch type constraints
#**********************************************************************/

# -----------------------------------------------------------------------
# GENERAL AND SCHEDULING RELATED SETS AND PARAMETERS
# -----------------------------------------------------------------------

#### General parameters

param BigNumber := 10000000;
param n_prds_per_day;
param n_days;
param n_weeks;
set PERIODS:={1..n_prds_per_day};
set DAYS:={1..n_days};
set WEEKS:={1..n_weeks};
set WEEKENDS:={1..2};

# The following gets used in the chain related constructions and may need further modification
param period_multiweek {(i,j,w) in {PERIODS,DAYS,WEEKS}}:=n_prds_per_day*n_days*(w-1)+n_prds_per_day*(j-1)+i;
param period {(i,j) in {PERIODS,DAYS}}:=n_prds_per_day*(j-1)+i;


#### Tour type related parameters

param n_lengths;                # Number of shift lengths
set LENGTHS:=1..n_lengths;
param lengths{LENGTHS}  ;  # Vector of shift lengths


param n_tts;                    # Number of different tour types
set TTYPES:=1..n_tts;
set tt_length_x{TTYPES};          # Set of allowable length indices by tour type

param num_weekend_patterns{TTYPES};       # Number of weekends worked patterns
param max_weekend_patterns := 2**(2*n_weeks);

# param A[i,j,w,t,e] = 1 if weekend pattern i calls for work on day j of week k for tour type t having weekend type e and 0 otherwise
param A{i in 1..max_weekend_patterns,j in DAYS,w in WEEKS,t in TTYPES,e in WEEKENDS: i <= num_weekend_patterns[t]} default 0;


param tt_min_dys_weeks{TTYPES,WEEKS} default 0;         # Minimum number of days worked by week by tour type
param tt_max_dys_weeks{TTYPES,WEEKS} default BigNumber;         # Maximum number of days worked by week by tour type

param tt_min_cumul_dys_weeks{TTYPES,WEEKS} default 0;         # Minimum number of days worked by cumulative weeks by tour type
param tt_max_cumul_dys_weeks{TTYPES,WEEKS} BigNumber 0;         # Maximum number of days worked by cumulative weeks by tour type

# Weekends consisting of a Fri and Sat imply updated lower bounds on some of the daily tour type variables. These were the key to modeling weekends
# worked patterns.
param FriSat_min_dys_weeks{t in TTYPES,1..num_weekend_patterns[t],WEEKS} default 0;
param FriSat_min_cumul_dys_weeks{t in TTYPES,1..num_weekend_patterns[t],WEEKS} default 0;


param tt_shiftlen_min_dys_weeks{TTYPES,LENGTHS,WEEKS} default 0;         # Minimum number of days worked by week by shiftlen by tour type
param tt_shiftlen_max_dys_weeks{TTYPES,LENGTHS,WEEKS} default 0;         # Maximum number of days worked by week by shiftlen by tour type

param tt_shiftlen_min_cumul_dys_weeks{TTYPES,LENGTHS,WEEKS} default 0;         # Minimum number of days worked by cumulative weeks by shiftlen by tour type
param tt_shiftlen_max_cumul_dys_weeks{TTYPES,LENGTHS,WEEKS} default 0;         # Maximum number of days worked by cumulative weeks by shiftlen by tour type

param tt_min_prds_weeks{TTYPES,WEEKS} default 0;        # Minimum number of periods worked by week by tour type
param tt_max_prds_weeks{TTYPES,WEEKS} default 0;         # Minimum number of periods worked by week by tour type

param tt_min_cumul_prds_weeks{TTYPES,WEEKS} default 0;         # Minimum number of periods worked by cumulative weeks by tour type
param tt_max_cumul_prds_weeks{TTYPES,WEEKS} default 0;         # Minimum number of periods worked by cumulative weeks by tour type



param tt_rel{TTYPES};           # Relationships from .MIX 
				#  (-1 --> "<=", 0 --> "=", 1 --> ">="
param tt_rhs{TTYPES};           # RHS from .MIX

### To the above, we'll add params and sets to allow direct modeling of side constraints
### of the form sum{subset of tour types} =, >=, <= some bound

param tt_parttime{TTYPES};      # 1 for part-time, 0 for full-time

set okTTYPES :=
 {t in TTYPES:tt_rel[t]<>0 or tt_rhs[t]<>0};
 
# Allowable shift start times  - note that these are tour type specific
param allow_start{PERIODS,DAYS,LENGTHS,TTYPES} default 0;
set ok_shifts:={i in PERIODS,j in DAYS,k in LENGTHS, t in okTTYPES: allow_start[i,j,k,t]>0};



				

# -----------------------------------------------------------------------
# START WINDOWS - should these be tour type specific since allow start is tour type specific?
#
# To start with,I made them week specific but not sure if they should be. Seems they should be 
# period, day, tour type.
#
# Maybe best strategy is to do the variables and constraints first and then work backwards to
# define windows appropriately.
# -----------------------------------------------------------------------

param width;                            # Width of start-time windows
param n_windows := n_prds_per_day;      # Number of start windows
set WINDOWS := {1..n_windows};

param midnight_thresh{TTYPES} default BigNumber;  # Need to generalize for any number of periods per day

set weekend{i in WINDOWS, t in TTYPES}:=
 { (if i+max{k in tt_length_x[t]}(lengths[k])-1>=midnight_thresh[t] then 6 else 1)}
  union {7};
  
  
/**** Beginning of each start window (in total periods from Sunday @ midnight)****/
###param b_window{i in PERIODS,j in DAYS,w in WEEKS} := n_prds_per_day*n_days*(w-1)+n_prds_per_day*(j-1)+i;

/**** Beginning of each start window (in total periods from Sunday @ midnight)****/
param b_window_wepoch{i in PERIODS,j in DAYS} := n_prds_per_day*(j-1)+i;

/**** End of each start window (in total periods from Sunday @ midnight) ****/
param e_window_wepoch{i in PERIODS,j in DAYS} :=
 ( if n_prds_per_day*(j-1)+i+width <= n_prds_per_day*n_days then
    n_prds_per_day*(j-1)+i+width
  else
    (n_prds_per_day*(j-1)+i+width )-n_prds_per_day*n_days);

#display b_window_wepoch;
#display e_window_wepoch;
#display {i in PERIODS,j in DAYS} e_window_wepoch[i,j]-b_window_wepoch[i,j];

/**** The set WindowWepochs{(i,j) in {PERIODS,DAYS}} contains all pairs 
   (l=period,m=day) within the start window which begins in period (i,j). ****/

set WindowWepochs{i in PERIODS,j in DAYS}
	   within {PERIODS,DAYS} := {(l,m) in {PERIODS,DAYS}:
(  (  (n_prds_per_day*(m-1)+l>=b_window_wepoch[i,j]) and
      (n_prds_per_day*(m-1)+l<=
	(if b_window_wepoch[i,j]<=e_window_wepoch[i,j] then e_window_wepoch[i,j]
	 else n_prds_per_day*n_days)
      )
   ) or
   ( (n_prds_per_day*(m-1)+l>=
       (if b_window_wepoch[i,j]<=e_window_wepoch[i,j] then b_window_wepoch[i,j]
	else 1)
     ) and (n_prds_per_day*(m-1)+l<=e_window_wepoch[i,j])
   )
)       };

#display WindowWepochs;

/**** The set okWindowWepochs{i in PERIODS,j in DAYS,k in LENGTHS,t in TTYPES } creates shift
length and tour type specific sets of windows that start in (i,j) pairs which have
(i,j,k) as an allowable shift start time. ****/

set okWindowWepochs{i in PERIODS,j in DAYS,k in LENGTHS,t in TTYPES}
 :={(l,m) in WindowWepochs[i,j]: allow_start[i,j,k,t]>0 and allow_start[l,m,k,t]>0};

#display okWindowWepochs;

/**** ok_window_beginnings is the set of start windows in which there is
at least one period in which a shift of length k can start
and which are not subsets of some other window. 

*/

set okWindowBeginnings{t in okTTYPES, k in tt_length_x[t]} :=
  setof{(p,q) in {PERIODS,DAYS}: (p,q,k,t) in ok_shifts and
   forall{(i,j) in {PERIODS,DAYS} diff {(p,q)}:allow_start[i,j,k,t]>0 }
    (not 
     ({(l,m) in okWindowWepochs[p,q,k,t]} 
      within {(n,o) in okWindowWepochs[i,j,k,t]}))  } (p,q);


#display okWindowBeginnings;

# -----------------------------------------------------------------------
# USEFUL RELATIONSHIPS BETWEEN PERIODS AND DAYS 
# USED IN START WINDOW RELATED CONSTRAINTS
# -----------------------------------------------------------------------


param which_prd{p in 1..(n_days+1)*n_prds_per_day} :=
   p-n_prds_per_day*(ceil(p/n_prds_per_day-1));

param which_day{p in 1..(n_days+1)*n_prds_per_day} :=
   (if p>n_prds_per_day*n_days then 1 else 1+ceil(p/n_prds_per_day-1));

# -----------------------------------------------------------------------
# CHAINS - NON-OVERLAPPING SEQUENCES OF (i,j) period,day PAIRS THAT
#          CAN BE ISOLATED FOR COORDINATING ix AND DWT VARIABLES.
#
# Let's wait on generalizing these for multiple weeks. Get model working
# for start window width of 0.
# -----------------------------------------------------------------------


set bchain {t in okTTYPES, k in tt_length_x[t]} := setof{(w,j) in (okWindowBeginnings[t,k]):
	forall{(p,q) in (okWindowBeginnings[t,k] diff {(w,j)})} (w,j) not in (okWindowWepochs[p,q,k,t])} (w,j) ;
	
#display bchain;


set echain {t in okTTYPES,k in LENGTHS,i in PERIODS,j in DAYS:
 (i,j) in bchain[t,k]} 
 := setof{(w,x) in {PERIODS,DAYS}: (w,x) not in (bchain[t,k] diff {(i,j)}) and
 (w,x) in okWindowBeginnings[t,k] and
 forall{(p,q) in (okWindowBeginnings[t,k] diff {(w,x)})} (p,q) not in (okWindowWepochs[w,x,k,t]) and
 	(period[w,x]>=period[i,j] and 
 	forall{(n,o) in bchain[t,k]: period[n,o]>period[i,j]} period[w,x]<
  period[n,o] or 
  (period[w,x]<period[i,j] and 
  forall{(n,o) in bchain[t,k]} (period[w,x]< 
   period[n,o] and period[n,o]<=period[i,j]) ) )
  } (w,x) ;
  
#display echain;


param numlinks{t in okTTYPES, k in tt_length_x[t],i in PERIODS,j in DAYS:(i,j) in bchain[t,k] }
 := sum{(l,m) in echain[t,k,i,j]}
  (if period[l,m]<period[i,j] then 
    n_prds_per_day*n_days+period[l,m] 
   else period[l,m]) - period[i,j]+1;

#display numlinks;

set chain{i in PERIODS, j in DAYS, t in okTTYPES, k in tt_length_x[t]: (i,j) in bchain[t,k]}
 := setof{(l,m) in {PERIODS,DAYS}: 
     ((l,m) in okWindowBeginnings[t,k] and 
     period[i,j]<=sum{(n,o) in echain[t,k,i,j]}period[n,o] and
     period[l,m]>=period[i,j] and 
     period[l,m]<=sum{(n,o) in echain[t,k,i,j]}period[n,o]) 
     or
     ((l,m) in okWindowBeginnings[t,k] and
     period[i,j]>sum{(n,o) in echain[t,k,i,j]}period[n,o] and
     ((period[l,m]>=period[i,j] and period[l,m]<=n_prds_per_day*n_days) or
     (period[l,m]<=sum{(n,o) in echain[t,k,i,j]}period[n,o] )) )
     } (l,m) ; 

#display chain;

/* NOTE: The 'dimen 2' in the following two sets is needed due to recursive nature of 
         their definitions. Both AMPL and GMPL have some problems if dimen 2 is omitted.
*/


set link{t in okTTYPES, k in tt_length_x[t], (i,j) in bchain[t,k], m in 1..numlinks[t,k,i,j] } dimen 2 := 
	setof{(n,o) in (if m=1 then chain[i,j,t,k] else (chain[i,j,t,k] diff link[t,k,i,j,m-1])): 
	forall{(p,d) in chain[i,j,t,k] : 
		period[p,d] > (if m=1 then 0 else max{(r,s) in link[t,k,i,j,m-1]} period[r,s]) } 
		period[n,o] <= period[p,d] and period[n,o] >= 
			(if m=1 then 0 else max{(r,s) in link[t,k,i,j,m-1]} period[r,s])} (n,o);
			
#display link;


set linkspan{t in okTTYPES, k in tt_length_x[t], (i,j) in bchain[t,k], m in 1..numlinks[t,k,i,j] } dimen 2 := 
	{(n,o) in {PERIODS,DAYS}: forall{(p,d) in link[t,k,i,j,m]} (n,o) in WindowWepochs[p,d] union
		if m=1 then WindowWepochs[p,d] else linkspan[t,k,i,j,m-1]};

#display linkspan;

# -----------------------------------------------------------------------
# VARIABLE DOMAINS        
# -----------------------------------------------------------------------



/*---------------------------------------------------------------------------------------------------- /
  OLD - The domain of the WT variables is defined to be those (window,ttype) pairs having at 
  least one days-off pattern with all window,day pairs (for days worked) in ok_window_beginnings.

  NEW - The domain of the WT variables is defined to be those (window,ttype) pairs having at 
  least tt_min_dys_week{ttype} in ok_window_beginnings[ttype].  
 ----------------------------------------------------------------------------------------------------*/

	  
/* set okTourType := setof{i in WINDOWS,t in okTTYPES :
	sum{j in DAYS} (if (i,j) in ok_window_beginnings[t] then 1 else 0) >= tt_min_dys_week[t]} (i,t);*/
	
set okTourType := setof{i in WINDOWS,t in okTTYPES} (i,t);

/*-----------------------------------------------------------------------------------------------------/
  OLD - The domain of the DWT variables is defined to be those (window,ttype,day) trips having at 
  least one days-off pattern calling for work on day p and that have with all 
  window,day pairs (for days worked) in ok_window_beginnings.
  
    NEW - The domain of the DWT variables is defined to be those (window,ttype,day) trips having 
          (window,day) in ok_window_beginnings[ttype]. 
 ----------------------------------------------------------------------------------------------------*/

/*set okDailyTourType := setof{i in WINDOWS,t in okTTYPES, d in DAYS :
	(i,d) in ok_window_beginnings[t]} (i,t,d);*/


set okDailyTourType := setof{i in WINDOWS,t in okTTYPES, d in DAYS} (i,t,d);

	  
# -----------------------------------------------------------------------
# SCHEDULING ENVIRONMENT PARAMETERS        
# -----------------------------------------------------------------------

# Target and minimum staffing levels - this is week specific. We can always allow user to input
# a single week and then repeat it for the other weeks.

param dmd_staff{PERIODS,DAYS,WEEKS};  
param tot_dmd := sum{(i,j,w) in {PERIODS,DAYS,WEEKS}} dmd_staff[i,j,w];
param min_staff{PERIODS,DAYS,WEEKS};  # Minimum staffing levels

# Limits on part time labor and limits on total labor

param ptfrac_tog;               # This will get set to one if there are part-time tour types   
param ptfrac;                   # Maximum fraction of labor hours covered by part-time employees
param budget;                   # Maximum labor expenditure


# -----------------------------------------------------------------------
# COST RELATED PARAMETERS
# -----------------------------------------------------------------------

param tt_diff{TTYPES};           # Tour type differential

param cu1;
param cu2;        
param usb;


# -----------------------------------------------------------------------
# DECISION VARIABLES
# -----------------------------------------------------------------------



##### Shift variables
	
var Shift{i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS,t in okTTYPES: (i,j,k,t) in ok_shifts} >= 0, integer;

	# Shift[i,j,w,k,t] = Number of shifts of length k starting in period i
	# of day j in week w for a tour of type t

##### Windowed tour type variables

var TourType{i in 1..n_windows,t in okTTYPES : (i,t) in okTourType}
	>= 0, integer;
	
	/* TourType[i,j] Number of employees working tour type j
	   starting in window i  */

##### Daily tour type variables

var DailyTourType{i in 1..n_windows,t in okTTYPES,d in DAYS,w in WEEKS : (i,t,d) in okDailyTourType}
	>= 0, integer;

	/* DailyTourType[i,t,d] Number of employees working tour type t
	   starting in window i and working day d in week w*/
	   
##### Days off variables  

/*  These do not exist in this model. It's only weekend days that get modelled explicitly. */

##### Daily shift worked variables

var DailyShiftWorked{i in 1..n_windows,t in okTTYPES,k in tt_length_x[t],d in DAYS,w in WEEKS : (i,t,d) in okDailyTourType}
	>= 0, integer;


##### Weekend Days off variables   

var WeekendDaysWorked{p in 1..max_weekend_patterns,i in 1..n_windows,t in okTTYPES : (i,t) in okTourType and p <= num_weekend_patterns[t]} 
       >= 0, integer ;
 
    # WeekendDaysWorked[d,i,t] = Number of employees working days-off patterns d
    # in start window i and of tour type t
	
	
	
##### Coverage related variables

	# For computational convenience, we broke up the calculation of coverage in each period into four separate cases related
	# to various types of overlap (or lack of). See coverage constraints.
	
var cov1{i in PERIODS,j in DAYS,w in WEEKS} >=0;
var cov2{i in PERIODS,j in DAYS,w in WEEKS} >=0;
var cov3{i in PERIODS,j in DAYS,w in WEEKS} >=0;
var cov4{i in PERIODS,j in DAYS,w in WEEKS} >=0;

var under1{i in PERIODS,j in DAYS,w in WEEKS} >=0,
				   <=min(dmd_staff[i,j,w]-
				   min(dmd_staff[i,j,w],min_staff[i,j,w]),usb);

var under2{i in PERIODS,j in DAYS,w in WEEKS} >=0
				   <=
      max(0,(if usb>0 then 1 else 0)*(dmd_staff[i,j,w]-
		     min(dmd_staff[i,j,w],min_staff[i,j,w])-usb));

# -----------------------------------------------------------------------
# OBJECTIVE FUNCTION
# -----------------------------------------------------------------------

minimize total_cost :

	sum{i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS,t in okTTYPES: (i,j,k,t) in ok_shifts} Shift[i,j,w,k,t]*lengths[k]*tt_diff[t] + 
         sum{i in PERIODS,j in DAYS,w in WEEKS}(under1[i,j,w]*cu1 + under2[i,j,w]*cu2);

# -----------------------------------------------------------------------
# CONSTRAINTS
# -----------------------------------------------------------------------


##### Budget constraints

subject to max_budget : 

	sum{i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS,t in TTYPES: (i,j,k,t) in ok_shifts} Shift[i,j,w,k,t]*lengths[k]*tt_diff[t] <= budget;

##### Coverage constraints

# Breaking them up into four different constraints, one for each case in terms of handling end of day horizon wrapping
subject to coverage1{i in PERIODS,j in DAYS,w in WEEKS} :
  cov1[i,j,w] = sum{l in LENGTHS,t in okTTYPES,p in 0..lengths[l]-1: (i-p)>0} if allow_start[(i-p),j,l,t]>0 then Shift[(i-p),j,w,l,t] else 0;
  
subject to coverage2{i in PERIODS,j in DAYS,w in WEEKS} :
  cov2[i,j,w] = sum{l in LENGTHS,t in okTTYPES,p in 0..lengths[l]-1: (i-p)<=0 and j>1} if allow_start[n_prds_per_day+(i-p),j-1,l,t]>0 then 
      Shift[n_prds_per_day+(i-p),j-1,w,l,t] else 0;
	  
subject to coverage3{i in PERIODS,j in DAYS,w in WEEKS} :
  cov3[i,j,w] = sum{l in LENGTHS,t in okTTYPES,p in 0..lengths[l]-1: (i-p)<=0 and j=1 and w>1} if allow_start[n_prds_per_day+(i-p),n_days,l,t]>0 then 
      Shift[n_prds_per_day+(i-p),n_days,w-1,l,t] else 0;
  
subject to coverage4{i in PERIODS,j in DAYS,w in WEEKS} :
  cov4[i,j,w] = sum{l in LENGTHS,t in okTTYPES,p in 0..lengths[l]-1: (i-p)<=0 and j=1 and w=1} if allow_start[n_prds_per_day+(i-p),n_days,l,t]>0 then 
    Shift[n_prds_per_day+(i-p),n_days,n_weeks,l,t] else 0;
	
# One final coverage constraint that combines the pieces	
subject to coverage{i in PERIODS,j in DAYS,w in WEEKS} :
  cov1[i,j,w] + cov2[i,j,w] + cov3[i,j,w] + cov4[i,j,w] + under1[i,j,w] + under2[i,j,w] >= dmd_staff[i,j,w];

# WT_bounds - bounds from .MIX file

subject to TourType_LB{t in TTYPES :tt_rhs[t]>=0}:
	  sum{(i,x) in okTourType : x=t}
	   TourType[i,t] >= ( if tt_rel[t] = -1 then 0 else tt_rhs[t] );

subject to TourType_UB{t in TTYPES :tt_rhs[t]>=0}:
	  sum{(i,x) in okTourType : x=t}
	  TourType[i,t] <= ( if tt_rel[t] = 1 then BigNumber else tt_rhs[t] ); 


# Chains - coordinates DWT and ix within each chain

# How to deal with multiple weeks?
# How to deal with fact that single tour type can have multiple shift lengths? It's a problem since we CANNOT equate sums
# of shift variables for a specific shift length with DailyTourType variables.

subject to chains_sweep_l{e in WEEKS, t in okTTYPES, k in tt_length_x[t], (b,j) in bchain[t,k],
			 p in period[b,j]..period[b,j],
			 w in 0..(numlinks[t,k,b,j]-1-(p-period[b,j])):
   width>0 } :
    sum{ (l,m) in linkspan[t,k,b,j,w+1]:
     (l,m,k,t) in ok_shifts} Shift[l,m,e,k,t] >=
     sum{u in p..p+w: (which_prd[u],which_day[u]) in okWindowBeginnings[t,k]
     and sum{(l,m) in WindowWepochs[which_prd[u],which_day[u]]} allow_start[l,m,k,t]>0} 
      DailyShiftWorked[which_prd[u],t,k,which_day[u],e] ; 


subject to chains_sweep_u{e in WEEKS,t in okTTYPES, k in tt_length_x[t],(b,j) in bchain[t,k],
			 w in 0..(numlinks[t,k,b,j]-1):
   width>0 } :
    sum{ i in period[b,j]..period[b,j]+w :
     (which_prd[i],which_day[i],k,t) in ok_shifts} 
      Shift[which_prd[i],which_day[i],e,k,t] <=
     sum{u in period[b,j]..period[b,j]+w : (which_prd[u],which_day[u]) in okWindowBeginnings[t,k] and 
      sum{(l,m) in WindowWepochs[which_prd[u],which_day[u]]} allow_start[l,m,k,t]>0} 
      DailyShiftWorked[which_prd[u],t,k,which_day[u],e] ; 




/**** ENABLED - AVOIDS USING ITERATED UNION ****/
/* NOTE: chains_tot constraints must appear after chains_sweep_l to avoid some
         sort of exception handling error in ampl.exe. It seems related to fact
         that linkspan is defined recursively and chain_sweep_l forces calculation
         of the terms in the natural order whereas chains_tot immediately asks
         for the 'last' element in the recursively defined set. */


subject to chains_tot{e in WEEKS, t in okTTYPES, k in tt_length_x[t],(i,j) in bchain[t,k]} :
 sum{(l,m) in linkspan[t,k,i,j,numlinks[t,k,i,j]]
   :(l,m,k,t) in ok_shifts}
    Shift[l,m,e,k,t] = sum{(n,o) in chain[i,j,t,k]:
      sum{(p,d) in WindowWepochs[n,o]} allow_start[p,d,k,t]>0
      } DailyShiftWorked[n,t,k,o,e];


# Integrate days-off and shift scheduling sub-problems

subject to weekend_integration_1{j in DAYS,w in WEEKS,i in 1..n_windows,t in okTTYPES : 
    j in weekend[i,t] and 1 in weekend[i,t] and 7 in weekend[i,t] and (i,t,j) in okDailyTourType} :
      sum{p in 1..num_weekend_patterns[t]} A[p,j,w,t,1]*WeekendDaysWorked[p,i,t] = DailyTourType[i,t,j,w];
	  
subject to weekend_integration_2{j in DAYS,w in WEEKS,i in 1..n_windows,
	t in okTTYPES : j in weekend[i,t] and 6 in weekend[i,t] and 7 in weekend[i,t] and (i,t,j) in okDailyTourType} :
      sum{p in 1..num_weekend_patterns[t]} A[p,j,w,t,2]*WeekendDaysWorked[p,i,t] = DailyTourType[i,t,j,w];

# Each tour variable must get a WeekendDaysWorked pattern assigned to it.
subject to weekend_total{(i,t) in okTourType} :
    (sum{p in 1..num_weekend_patterns[t]} WeekendDaysWorked[p,i,t]) - TourType[i,t]  = 0;
   
# .............................................................................................
# Integrate shift, days worked, and tour type variables 
# .............................................................................................

# Make sure number of shifts scheduled each day for each tour type doesn't exceed the number of 
# tour types scheduled to work that day

/*subject to shift_DWT_UB{j in DAYS, i in 1..n_windows, t in TTYPES, w in WEEKS : (i,t,j) in okDailyTourType} :
    sum{k in LENGTHS: k in tt_length_x[t] and (i,j,k,t) in ok_shifts} Shift[i,j,w,k,t] <= DailyTourType[i,t,j,w];
*/
/*
subject to shift_DWT{j in DAYS, i in 1..n_windows, t in TTYPES, w in WEEKS : (i,t,j) in okDailyTourType} :
    sum{k in LENGTHS: k in tt_length_x[t] and (i,j,k,t) in ok_shifts} Shift[i,j,w,k,t] = DailyTourType[i,t,j,w];
*/
# ---- Make sure number of shifts scheduled each week for each tour is equal to the number of days worked each week

subject to shift_DWT_dailyconservation{i in 1..n_windows, j in DAYS, w in WEEKS, t in TTYPES : (i,t,j) in okDailyTourType } :
    sum{k in LENGTHS: k in tt_length_x[t] and (i,j,k,t) in ok_shifts} Shift[i,j,w,k,t] = DailyTourType[i,t,j,w];

# ------------------------------------------------------------------------------------------------

# ---- Make sure number of people working on any given day <= number of employees
subject to DWT_WT_UB{j in DAYS,i in 1..n_windows, t in okTTYPES, w in WEEKS : (i,t,j) in okDailyTourType} :
    DailyTourType[i,t,j,w] <= TourType[i,t];
	
	
subject to DWT_DSW{j in DAYS,i in 1..n_windows, t in okTTYPES, w in WEEKS : (i,t,j) in okDailyTourType} :
    sum{k in tt_length_x[t]}DailyShiftWorked[i,t,k,j,w] = DailyTourType[i,t,j,w];
# ------------------------------------------------------------------------------------------------	
	
	
# ---- Min and max bounds on periods worked over cumulative number of weeks	

# ---- Min and max bounds on periods worked each week	
subject to prds_worked_weekly_LB{i in 1..n_windows, t in okTTYPES, w in WEEKS } :
    sum{j in DAYS,k in LENGTHS: (i,j,k,t) in ok_shifts} Shift[i,j,w,k,t]*lengths[k] >= TourType[i,t]*tt_min_prds_weeks[t,w];
	
subject to prds_worked_weekly_UB{i in 1..n_windows, t in okTTYPES, w in WEEKS } :
    sum{j in DAYS,k in LENGTHS: (i,j,k,t) in ok_shifts} Shift[i,j,w,k,t]*lengths[k] <= TourType[i,t]*tt_max_prds_weeks[t,w];
	
subject to prds_worked_cumul_weekly_LB{i in 1..n_windows, t in okTTYPES, w in WEEKS } :
    sum{j in DAYS,k in LENGTHS,z in 1..w: (i,j,k,t) in ok_shifts} Shift[i,j,z,k,t]*lengths[k] >= TourType[i,t]*tt_min_cumul_prds_weeks[t,w];
	
subject to prds_worked_cumul_weekly_UB{i in 1..n_windows, t in okTTYPES, w in WEEKS } :
    sum{j in DAYS,k in LENGTHS,z in 1..w: (i,j,k,t) in ok_shifts} Shift[i,j,z,k,t]*lengths[k] <= TourType[i,t]*tt_max_cumul_prds_weeks[t,w];
	


subject to prds_ID_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, w in WEEKS:
 6 in weekend[i,t] and 7 in weekend[i,t]} :
    sum{j in DAYS,k in LENGTHS: (i,j,k,t) in ok_shifts} Shift[i,j,w,k,t]*lengths[k] >= 
	sum{p in 1..7}WeekendDaysWorked[p,i,t]*FriSat_min_dys_weeks[t,p,w]*min{k in tt_length_x[t]}lengths[k];
	
subject to prds_ID_cumul_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, w in WEEKS:
 6 in weekend[i,t] and 7 in weekend[i,t]} :
    sum{j in DAYS,k in LENGTHS,z in 1..w: (i,j,k,t) in ok_shifts} Shift[i,j,z,k,t]*lengths[k] >= 
	sum{p in 1..7}WeekendDaysWorked[p,i,t]*FriSat_min_cumul_dys_weeks[t,p,w]*min{k in tt_length_x[t]}lengths[k];

# ------------------------------------------------------------------------------------------------

# ---- OK - Min and max bounds on shifts worked each week
subject to shift_WT_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, k in LENGTHS, w in WEEKS } :
    sum{d in DAYS: (i,d,k,t) in ok_shifts} Shift[i,d,w,k,t] >= TourType[i,t]*tt_shiftlen_min_dys_weeks[t,k,w];
	
subject to shift_WT_weeklyconservation_UB{i in 1..n_windows, t in okTTYPES, k in LENGTHS, w in WEEKS } :
    sum{d in DAYS: (i,d,k,t) in ok_shifts} Shift[i,d,w,k,t] <= TourType[i,t]*tt_shiftlen_max_dys_weeks[t,k,w];


# ---- OK - Min and max bounds on shifts worked over cumulative number of weeks
subject to shift_WT_cumul_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, k in LENGTHS, w in WEEKS } :
    sum{d in DAYS,z in 1..w: (i,d,k,t) in ok_shifts} Shift[i,d,z,k,t] >= TourType[i,t]*tt_shiftlen_min_cumul_dys_weeks[t,k,w];
	
subject to shift_WT_cumul_weeklyconservation_UB{i in 1..n_windows, t in okTTYPES, k in LENGTHS, w in WEEKS } :
    sum{d in DAYS,z in 1..w: (i,d,k,t) in ok_shifts} Shift[i,d,z,k,t] <= TourType[i,t]*tt_shiftlen_max_cumul_dys_weeks[t,k,w];




# ---- Min and max bounds on days worked each week

subject to DWT_WT_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, w in WEEKS } :
    sum{d in DAYS: (i,t,d) in okDailyTourType} DailyTourType[i,t,d,w] >= TourType[i,t]*tt_min_dys_weeks[t,w];
	
subject to DWT_WT_weeklyconservation_UB{i in 1..n_windows, t in okTTYPES, w in WEEKS } :
    sum{d in DAYS: (i,t,d) in okDailyTourType} DailyTourType[i,t,d,w] <= TourType[i,t]*tt_max_dys_weeks[t,w];

subject to DWT_ID_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, w in WEEKS:
 6 in weekend[i,t] and 7 in weekend[i,t]} :
    sum{d in DAYS: (i,t,d) in okDailyTourType} DailyTourType[i,t,d,w] >= 
	sum{p in 1..num_weekend_patterns[t]}WeekendDaysWorked[p,i,t]*FriSat_min_dys_weeks[t,p,w];
	
subject to DWT_ID_cumul_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, w in WEEKS:
 6 in weekend[i,t] and 7 in weekend[i,t]} :
    sum{d in DAYS, z in 1..w: (i,t,d) in okDailyTourType} DailyTourType[i,t,d,z] >= 
	sum{p in 1..num_weekend_patterns[t]}WeekendDaysWorked[p,i,t]*FriSat_min_cumul_dys_weeks[t,p,w];
	

# ---- Min and max bounds on days worked over cumulative number of weeks

subject to DWT_WT_cumul_weeklyconservation_LB{i in 1..n_windows, t in okTTYPES, w in WEEKS } :
    sum{d in DAYS,z in 1..w: (i,t,d) in okDailyTourType} DailyTourType[i,t,d,z] >= TourType[i,t]*tt_min_cumul_dys_weeks[t,w];
	
subject to DWT_WT_cumul_weeklyconservation_UB{i in 1..n_windows, t in okTTYPES, w in WEEKS } :
    sum{d in DAYS,z in 1..w: (i,t,d) in okDailyTourType} DailyTourType[i,t,d,z] <= TourType[i,t]*tt_max_cumul_dys_weeks[t,w];


# Part-time fraction constraints


subject to max_ptfrac :
   (  sum{i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS,t in TTYPES: ptfrac_tog > 0 and (i,j,k,t) in ok_shifts} 
	   ( if tt_parttime[t] > 0 then 
	    Shift[i,j,w,k,t]*lengths[k]*(1.0-1.0/ptfrac) else
	        Shift[i,j,w,k,t]*lengths[k] ) ) >= 0.0 ;
			


/*******************************************************************/

solve;

printf "# Start TourType\n";
printf "param\nTourType";
for {t in okTTYPES}
{
	
	printf " [*,%i]\n",t;
	printf "   :=\n";

	for {w in WINDOWS: (w,t) in okTourType and TourType[w,t]>0}
	{
	    printf "%5i%5i",w,TourType[w,t];	
	    printf "\n";
	}

}
printf ";\n\n";


printf "\n";
printf "# Start WeekendDaysWorked\n";
printf "param\nWeekendDaysWorked";
for {t in TTYPES}
{
	
	printf " [*,*,%i] (tr)\n",t;
	printf "   :";
	for {p in (1..num_weekend_patterns[t])}
	{
		printf "%4i",p;
	}
	printf "  := \n";



	for {i in PERIODS: (i,t) in okTourType}
	{
		printf "%3i",i;
		for {p in (1..num_weekend_patterns[t])}
		{
			printf "%5i",WeekendDaysWorked[p,i,t];	
		}
		printf "\n";
	}

}

printf ";\n\n";

printf "# Start DailyTourType\n";
printf "param\nDailyTourType";
for {w in WEEKS}
{
for {t in okTTYPES}
{
	
	printf " [*,%i,*,%i]\n",t,w;
	printf "   :";
	for {j in DAYS}
	{
		printf "%4i",j;
	}
	printf "  := \n";



	for {i in WINDOWS: exists{d in DAYS}(i,t,d) in okDailyTourType}
	{
		printf "%5i",i;
		for {j in DAYS}
		{
			printf "%5i",DailyTourType[i,t,j,w];	
		}
		printf "\n";
	}

}
}
printf ";\n\n";




printf "# Start shift\n";
printf "param\nShift";
for{w in WEEKS}
{
for{t in TTYPES}
{
for {k in LENGTHS}
{
	
	printf " [*,*,%i,%i,%i]\n",w,k,t;
	printf "   :";
	for {j in DAYS}
	{
		printf "%4i",j;
	}
	printf "  := \n";



	for {i in PERIODS: exists{d in DAYS}(i,d,k,t) in ok_shifts}
	{
		printf "%3i",i;
		for {j in DAYS}
		{
			printf "%5i",Shift[i,j,w,k,t];	
		}
		printf "\n";
	}

}
}
}
printf ";\n\n";

printf "/*";
#display DailyShiftWorked;
printf "*/";



printf "/*";
printf "# Start shift\n";

for{w in WEEKS}
{
for{t in TTYPES}
{
for {k in LENGTHS}
{
    
	for {i in PERIODS: exists{d in DAYS}(i,d,k,t) in ok_shifts}
	{
		
		for {j in DAYS}
		{
		    
			
			for {s in 1..Shift[i,j,w,k,t]}
			{
			printf "%i,%i,%i,%i,%i\n",i,j,w,k,t;
            }
            			
		}
		
	}

}
}
}
printf "*/";
printf "\n\n";

end;

