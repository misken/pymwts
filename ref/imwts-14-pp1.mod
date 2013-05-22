#/***********************************************************************
#                                                                      
#  Program:  imwts-14-pp1.mod
#  Language: AMPL/GMPL                                                       
#  Author:   Mark Isken                                                 
#  Create Date:     11/10/2011
#
#  btspp1-mp for post-processing btswidth-mp.mod                       *
#  solutions to create tour schedules.                                 *
#                                                                      *
#  Days off patterns NOT explicitly considered.                            *
#                                                                      *
#***********************************************************************
# -----------------------------------------------------------------------
# GENERAL AND SCHEDULING RELATED SETS AND PARAMETERS
# -----------------------------------------------------------------------

#### General parameters

param BigNumber := 10000000;
param n_prds_per_day;
param n_days_per_week;
param n_weeks;
set PERIODS:={1..n_prds_per_day};
set DAYS:={1..n_days_per_week};
set WEEKS:={1..n_weeks};
set WEEKENDS:={1..2};

# The following gets used in the chain related constructions and may need further modification
param period_multiweek {(i,j,w) in {PERIODS,DAYS,WEEKS}}:=n_prds_per_day*n_days_per_week*(w-1)+n_prds_per_day*(j-1)+i;
param period {(i,j) in {PERIODS,DAYS}}:=n_prds_per_day*(j-1)+i;

#### Tour type related parameters

param n_lengths;                # Number of shift lengths
set LENGTHS:=1..n_lengths;
param lengths{LENGTHS}  ;  # Vector of shift lengths


param n_tts;                    # Number of different tour types
set TTYPES:=1..n_tts;
set tt_length_x{TTYPES};          # Set of allowable length indices by tour type

param max_weekend_patterns := 2**(2*n_weeks);
param num_weekend_patterns{WEEKENDS,TTYPES} default 0;       # Number of weekends worked patterns


# param A[i,j,w,t,e] = 1 if weekend pattern i calls for work on day j of week k for tour type t having weekend type e and 0 otherwise
param A{i in 1..max_weekend_patterns,j in DAYS,w in WEEKS,t in TTYPES,e in WEEKENDS: i <= num_weekend_patterns[e, t]} default 0;

param tt_min_dys_weeks{TTYPES,WEEKS} default 0;         # Minimum number of days worked by week by tour type
param tt_max_dys_weeks{TTYPES,WEEKS} default BigNumber;         # Maximum number of days worked by week by tour type

param tt_min_cumul_dys_weeks{TTYPES,WEEKS} default 0;         # Minimum number of days worked by cumulative weeks by tour type
param tt_max_cumul_dys_weeks{TTYPES,WEEKS} default BigNumber;         # Maximum number of days worked by cumulative weeks by tour type

# Weekends consisting of a Fri and Sat imply updated lower bounds on some of the daily tour type variables. These were the key to modeling weekends
# worked patterns.
param FriSat_min_dys_weeks{t in TTYPES,1..num_weekend_patterns[2,t],WEEKS} default 0;
param FriSat_min_cumul_dys_weeks{t in TTYPES,1..num_weekend_patterns[2,t],WEEKS} default 0;


param tt_shiftlen_min_dys_weeks{TTYPES,LENGTHS,WEEKS} default 0;         # Minimum number of days worked by week by shiftlen by tour type
param tt_shiftlen_max_dys_weeks{TTYPES,LENGTHS,WEEKS} default 0;         # Maximum number of days worked by week by shiftlen by tour type

param tt_shiftlen_min_cumul_dys_weeks{TTYPES,LENGTHS,WEEKS} default 0;         # Minimum number of days worked by cumulative weeks by shiftlen by tour type
param tt_shiftlen_max_cumul_dys_weeks{TTYPES,LENGTHS,WEEKS} default 0;         # Maximum number of days worked by cumulative weeks by shiftlen by tour type

param tt_min_prds_weeks{TTYPES,WEEKS} default 0;        # Minimum number of periods worked by week by tour type
param tt_max_prds_weeks{TTYPES,WEEKS} default 0;         # Minimum number of periods worked by week by tour type

param tt_min_cumul_prds_weeks{TTYPES,WEEKS} default 0;         # Minimum number of periods worked by cumulative weeks by tour type
param tt_max_cumul_prds_weeks{TTYPES,WEEKS} default 0;         # Minimum number of periods worked by cumulative weeks by tour type

param tt_shiftlen_min_prds_weeks{TTYPES,LENGTHS,WEEKS} default 0;        # Minimum number of periods worked by week by tour type by shift length
param tt_shiftlen_max_prds_weeks{TTYPES,LENGTHS,WEEKS} default 0;         # Minimum number of periods worked by week by tour type by shift length

param tt_shiftlen_min_cumul_prds_weeks{TTYPES,LENGTHS,WEEKS} default 0;         # Minimum number of periods worked by cumulative weeks by tour type by shift length
param tt_shiftlen_max_cumul_prds_weeks{TTYPES,LENGTHS,WEEKS} default 0;         # Minimum number of periods worked by cumulative weeks by tour type by shift length


param tt_lb{TTYPES} default 0;           # RHS from .MIX
param tt_ub{TTYPES} default BigNumber;
### To the above, we'll add params and sets to allow direct modeling of side constraints
### of the form sum{subset of tour types} =, >=, <= some bound


param tt_parttime{TTYPES};      # 1 for part-time, 0 for full-time

set okTTYPES := {t in TTYPES};  # A placeholder for now
 
# Allowable shift start times  - note that these are tour type specific
param allow_start{PERIODS,DAYS,LENGTHS,TTYPES} default 0;
set ok_shifts:={i in PERIODS,j in DAYS,k in LENGTHS, t in okTTYPES: allow_start[i,j,k,t]>0};


param g_start_window_width;                            # Width of start-time windows
param n_windows := n_prds_per_day;      # Number of start windows
set WINDOWS := {1..n_windows};

param midnight_thresh{TTYPES} default BigNumber;  # Need to generalize for any number of periods per day

set weekend{i in WINDOWS, t in TTYPES}:=
 { (if i+max{k in tt_length_x[t]}(lengths[k])-1>=midnight_thresh[t] then 6 else 1)}
  union {7};
  
param weekend_type{i in WINDOWS, t in TTYPES}:=
  (if i+max{k in tt_length_x[t]}(lengths[k])-1>=midnight_thresh[t] then 2 else 1);



/**** Beginning of each start window (in total periods from Sunday @ midnight)****/
###param b_window{i in PERIODS,j in DAYS,w in WEEKS} := n_prds_per_day*n_days_per_week*(w-1)+n_prds_per_day*(j-1)+i;

/**** Beginning of each start window (in total periods from Sunday @ midnight)****/
param b_window_wepoch{i in PERIODS,j in DAYS} := n_prds_per_day*(j-1)+i;

/**** End of each start window (in total periods from Sunday @ midnight) ****/
param e_window_wepoch{i in PERIODS,j in DAYS} :=
 ( if n_prds_per_day*(j-1)+i+g_start_window_width <= n_prds_per_day*n_days_per_week then
    n_prds_per_day*(j-1)+i+g_start_window_width
  else
    (n_prds_per_day*(j-1)+i+g_start_window_width )-n_prds_per_day*n_days_per_week);

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
	 else n_prds_per_day*n_days_per_week)
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


param which_prd{p in 1..(n_days_per_week+1)*n_prds_per_day} :=
   p-n_prds_per_day*(ceil(p/n_prds_per_day-1));

param which_day{p in 1..(n_days_per_week+1)*n_prds_per_day} :=
   (if p>n_prds_per_day*n_days_per_week then 1 else 1+ceil(p/n_prds_per_day-1));
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
    n_prds_per_day*n_days_per_week+period[l,m] 
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
     ((period[l,m]>=period[i,j] and period[l,m]<=n_prds_per_day*n_days_per_week) or
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

# Do I need this???
set window_shifts{i in WINDOWS, k in LENGTHS, t in TTYPES: k in tt_length_x[t]} :=
 {(p,d) in {PERIODS,DAYS}: allow_start[p,d,k,t]>0 and 
     ((p,d) in WindowWepochs[i,d] or 
      (d>1 and (p,d) in WindowWepochs[i,d-1]) or
      (d=1 and (p,d) in WindowWepochs[i,n_days_per_week]))};
     


# -----------------------------------------------------------------------
# VARIABLE DOMAINS        
# -----------------------------------------------------------------------



/*---------------------------------------------------------------------------------------------------- /
  OLD - The domain of the TourType variables is defined to be those (window,ttype) pairs having at 
  least one days-off pattern with all window,day pairs (for days worked) in ok_window_beginnings.

  NEW - The domain of the TourType variables is defined to be those (window,ttype) pairs having at 
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

param max_parttime_frac_tog;               # This will get set to one if there are part-time tour types   
param max_parttime_frac;                   # Maximum fraction of labor hours covered by part-time employees
param labor_budget;                   # Maximum labor expenditure


# -----------------------------------------------------------------------
# COST RELATED PARAMETERS
# -----------------------------------------------------------------------

param tt_cost_multiplier{TTYPES};           # Tour type differential

param cu1;
param cu2;        
param usb;



# -----------------------------------------------------------------------
# DECISION VARIABLES
# -----------------------------------------------------------------------



##### Shift variables
	
param Shift{i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS,t in okTTYPES: (i,j,k,t) in ok_shifts} default 0;

	# Shift[i,j,w,k,t] = Number of shifts of length k starting in period i
	# of day j in week w for a tour of type t

##### Windowed tour type variables

param TourType{i in 1..n_windows,t in okTTYPES : (i,t) in okTourType} default 0;
	
	/* TourType[i,j] Number of employees working tour type j
	   starting in window i  */

##### Daily tour type variables

param DailyTourType{i in 1..n_windows,t in okTTYPES,d in DAYS,w in WEEKS : (i,t,d) in okDailyTourType} default 0;

	/* DailyTourType[i,t,d] Number of employees working tour type t
	   starting in window i and working day d in week w*/
	   
##### Days off variables  

/*  These do not exist in this model. It's only weekend days that get modelled explicitly. */

##### Daily shift worked variables

param DailyShiftWorked{i in 1..n_windows,t in okTTYPES,k in tt_length_x[t],d in DAYS,w in WEEKS : (i,t,d) in okDailyTourType} default 0;


##### Weekend Days off variables   

param WeekendDaysWorked{p in 1..max_weekend_patterns,i in 1..n_windows,t in okTTYPES : (i,t) in okTourType and p <= num_weekend_patterns[weekend_type[i,t],t]} 
       default 0;
       
       
       
       
       
param n_tours:=sum{i in 1..n_windows,t in okTTYPES :(i,t) in okTourType } TourType[i,t];

set index_tours{i in 1..n_windows,t in okTTYPES:(i,t) in okTourType }  :=
 { (sum{k in 1..i, j in 1..(if k<i then n_tts else t):
	 j in okTTYPES and (k,j) in okTourType}(TourType[k,j])-TourType[i,t]+1)..
   (sum{k in 1..i, j in 1..(if k<i then n_tts else t):
	 j in okTTYPES and (k,j) in okTourType}TourType[k,j]) };


set WIN_sx{l in 1..n_tours}:={i in 1..n_windows:
    exists{t in okTTYPES} (i,t) in okTourType
      and l in index_tours[i,t]};


set TT_sx{l in 1..n_tours}:={t in okTTYPES:
    exists{(i,t) in okTourType}
      l in index_tours[i,t]};

 param WIN_x{l in 1..n_tours}:=min{i in WIN_sx[l]} i;
 param TT_x{l in 1..n_tours}:=min{i in TT_sx[l]} i;


#......................................................Dec. Variables

#TODO - generalize for start windows
var tourshift{s in 1..n_tours, i in PERIODS, j in DAYS, w in WEEKS, k in LENGTHS, t in TTYPES, p in PERIODS, d in DAYS :
       t=TT_x[s] and p=WIN_x[s] and k in tt_length_x[t] and (i,j,k,t) in ok_shifts and 
       (i,j) in okWindowWepochs[p,d,k,t]} binary ;
    
    # tourshift[s,i,j,w,k,t,p,d] = 1 if an x[i,j,w,k,t] shift is assigned to tour s in window (p,d)
    #            = 0 otherwise
 


var tourdof{l in 1..n_tours,d in 1..max_weekend_patterns,i in WINDOWS,t in okTTYPES :
       t=TT_x[l] and i=WIN_x[l] and (i,t) in okTourType and d <= num_weekend_patterns[weekend_type[i,t],t] } binary ;
    
    # tourdof[l,d,i,t] = 1 if an ID[d,i,t] days off variable is assigned to tour l
    #            = 0 otherwise
    # tourdof ==> "binary days-off"  


#.......................................................Obj. Function

# Objective function is just the sum of the tourshift variables - which makes it a search for a feasible solution.

minimize total_cost :
 	   
 sum {s in 1..n_tours,i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS, t in TTYPES, p in PERIODS, d in DAYS :
       t=TT_x[s] and p=WIN_x[s] and k in tt_length_x[t] and (i,j,k,t) in ok_shifts and 
       (i,j) in okWindowWepochs[p,d,k,t]} tourshift[s,i,j,w,k,t,p,d] ;
		
		
		

# Weekend days worked for each tour determines number of weekend shifts assigned to each tour

subject to tourshift_tourdof_integration1{s in 1..n_tours,i in WINDOWS,j in DAYS, w in WEEKS,t in TTYPES: 
    t=TT_x[s] and i=WIN_x[s] and j in weekend[i,TT_x[s]]} :        
     sum{k in tt_length_x[t],(p,d) in okWindowWepochs[i,j,k,t] } tourshift[s,p,d,w,k,t,i,j] =
        (sum{g in 1..num_weekend_patterns[weekend_type[i,t],t]}tourdof[s,g,i,t]*A[g,j,w,t,weekend_type[i,t]]);


# One days off pattern for each tour 

subject to one_tourdof{s in 1..n_tours} :
	sum{d in 1..num_weekend_patterns[weekend_type[WIN_x[s],TT_x[s]],TT_x[s]]} tourdof[s,d,WIN_x[s],TT_x[s]]=1;


# All days off patterns assigned 

subject to conserve_ID {d in 1..max_weekend_patterns,i in WINDOWS,t in okTTYPES: d <= num_weekend_patterns[weekend_type[i,t],t] and (i,t) in okTourType} :
 sum{s in 1..n_tours : t=TT_x[s] and i=WIN_x[s]} 
        tourdof[s,d,i,t]= WeekendDaysWorked[d,i,t];


# No more than one shift worked per day

subject to tours_daily{s in 1..n_tours,j in DAYS,w in WEEKS} :
  sum{k in tt_length_x[TT_x[s]],(p,d) in okWindowWepochs[WIN_x[s],j,k,TT_x[s]]}
   tourshift[s,p,d,w,k,TT_x[s],WIN_x[s],j] <= 1;
   
 
# Needs to be generalized for start windows
subject to tours_daily_conservation2{i in PERIODS, j in DAYS,w in WEEKS, k in LENGTHS, t in TTYPES: k in tt_length_x[t] and (i,j,k,t) in ok_shifts} :
  sum{s in 1..n_tours: TT_x[s]=t and WIN_x[s]=i }
   tourshift[s,i,j,w,k,t,WIN_x[s],j] = Shift[i,j,w,k,t];  
 


# All the shifts need to get assigned to the tours

/*
subject to tours_daily_conservation2{p in PERIODS, d in DAYS,w in WEEKS, k in LENGTHS, t in TTYPES: k in tt_length_x[t]} :
  sum{s in 1..n_tours, (i,j) in okWindowWepochs[p,d,k,t]: TT_x[s]=t and WIN_x[s]=p }
   tourshift[s,i,j,w,k,t,p,d] = sum{(i,j) in okWindowWepochs[p,d,k,t] : (i,j,k,t) in ok_shifts}Shift[i,j,w,k,t]; 
*/



# Tour type specific bounds on number of days worked over the weeks (both cumulative and non-cumulative)
subject to tours_weekly_LB{t in 1..n_tours,w in WEEKS} :
  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_min_dys_weeks[TT_x[t],w];
   
subject to tours_weekly_UB{t in 1..n_tours,w in WEEKS} :
  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_max_dys_weeks[TT_x[t],w];
   
   

 
 
 #Seems like the following two constraints need to be generalized for > 2 weeks.

subject to tours_tot_LB{t in 1..n_tours} :
  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_min_cumul_dys_weeks[TT_x[t],n_weeks];
   
subject to tours_tot_UB{t in 1..n_tours} :
  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_max_cumul_dys_weeks[TT_x[t],n_weeks];
   
   

 
 
 
# Tour type and shift length specific bounds on number of days worked over the weeks (both cumulative and non-cumulative)
subject to tours_shiftlen_weekly_LB{t in 1..n_tours,k in tt_length_x[TT_x[t]],w in WEEKS} :
  sum{d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_shiftlen_min_dys_weeks[TT_x[t],k,w];
   
subject to tours_shiftlen_weekly_UB{t in 1..n_tours,k in tt_length_x[TT_x[t]],w in WEEKS} :
  sum{d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_shiftlen_max_dys_weeks[TT_x[t],k,w];
   
subject to tours_shiftlen_tot_LB{t in 1..n_tours,k in tt_length_x[TT_x[t]]} :
  sum{w in WEEKS,d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] >= tt_shiftlen_min_cumul_dys_weeks[TT_x[t],k,n_weeks];
   
subject to tours_shiftlen_tot_UB{t in 1..n_tours,k in tt_length_x[TT_x[t]]} :
  sum{w in WEEKS,d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d] <= tt_shiftlen_max_cumul_dys_weeks[TT_x[t],k,n_weeks]; 
 
 

  
  
# Tour type specific bounds on number of periods worked over the weeks (both cumulative and non-cumulative)
subject to tours_weekly_prds_LB{t in 1..n_tours,w in WEEKS} :
  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] >= tt_min_prds_weeks[TT_x[t],w];
   
subject to tours_weekly_prds_UB{t in 1..n_tours,w in WEEKS} :
  sum{k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] <= tt_max_prds_weeks[TT_x[t],w];
   
subject to tours_tot_prds_LB{t in 1..n_tours} :
  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] >= tt_min_cumul_prds_weeks[TT_x[t],n_weeks];
   
subject to tours_tot_prds_UB{t in 1..n_tours} :
  sum{w in WEEKS,k in tt_length_x[TT_x[t]],d in DAYS,(i,j) in okWindowWepochs[WIN_x[t],d,k,TT_x[t]] }
   tourshift[t,i,j,w,k,TT_x[t],WIN_x[t],d]*lengths[k] <= tt_max_cumul_prds_weeks[TT_x[t],n_weeks];  
  
  #####################################################################  
  
 # Ensures that each day worked in tour l is assigned exactly one shift
# of an appropriate length and from an appropriate start window.


subject to tot_conserve_shifts
   :
  sum {s in 1..n_tours,i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS, t in TTYPES, p in PERIODS, d in DAYS :
       t=TT_x[s] and p=WIN_x[s] and k in tt_length_x[t] and (i,j,k,t) in ok_shifts and 
       (i,j) in okWindowWepochs[p,d,k,t]} tourshift[s,i,j,w,k,t,p,d] =
   sum{m in PERIODS,n in DAYS,z in WEEKS, o in LENGTHS, t in TTYPES : 
   (m,n,o,t) in ok_shifts and Shift[m,n,z,o,t]>0} Shift[m,n,z,o,t];


#display tourshift;

solve;

#display tourshift;
/*
printf{p in PERIODS, d in DAYS,w in WEEKS, k in LENGTHS, t in TTYPES: k in tt_length_x[t]} :
  'w=%2i k=%2i t=%2i p=%2i d=%2i sumtours=%3i\n',w,k,t,p,d,sum{s in 1..n_tours, (i,j) in okWindowWepochs[p,d,k,t]: TT_x[s]=t and WIN_x[s]=p }
   tourshift[s,i,j,w,k,t,p,d] ;
   
printf{p in PERIODS, d in DAYS,w in WEEKS, k in LENGTHS, t in TTYPES: k in tt_length_x[t]} :
  'w=%2i k=%2i t=%2i p=%2i d=%2i sumshifts=%3i\n',w,k,t,p,d,sum{(i,j) in okWindowWepochs[p,d,k,t] : (i,j,k,t) in ok_shifts}Shift[i,j,w,k,t];   
    
*/
printf 'n_weeks %2i\n',n_weeks;
printf 'n_tours %4i\n',n_tours;

printf 'PP4\n';
printf {l in 1..n_tours, w in WEEKS, i in PERIODS, j in DAYS, d in DAYS, p in PERIODS, k in LENGTHS, t in TTYPES  : 
	k in tt_length_x[TT_x[l]] and t=TT_x[l] and 
     (i,j) in okWindowWepochs[WIN_x[l],d,k,TT_x[l]] and p=WIN_x[l] and tourshift[l,i,j,w,k,t,p,d] = 1}: 
	'%3i%3i%3i%3i%3i%3i%3i%3i\n',l,t,w,i,j,lengths[k],p,d;

/*
display {l in 1..n_tours,i in PERIODS,j in DAYS,w in WEEKS,k in LENGTHS, t in TTYPES : 
       t=TT_x[l] and (i,j,k,TT_x[l]) in ok_shifts and shift[i,j,w,k,TT_x[l]]>0 and k in tt_length_x[TT_x[l]] and 
        (i,j) in window_shifts[WIN_x[l],k,TT_x[l]] and tourshift[l,i,j,w,k,t]>0} tourshift[l,i,j,w,k,t] ;
*/
printf '\nTTS\n';
printf {t in 1..n_tours}: '%3i\n', TT_x[t];

/*
printf 'PP2\n';
printf {t in 1..n_tours,w in WEEKS}: '%3i%3i%3i%3i%3i%3i%3i%3i%4i%3i%3i%3i%3i%3i%3i%3i\n', 
  w,start_time[t,1,w], start_time[t,2,w], start_time[t,3,w], start_time[t,4,w], start_time[t,5,w], start_time[t,6,w], start_time[t,7,w],
  1, 
  shift_length[t,1,w], shift_length[t,2,w], shift_length[t,3,w], shift_length[t,4,w], shift_length[t,5,w], shift_length[t,6,w], shift_length[t,7,w];
*/

end;