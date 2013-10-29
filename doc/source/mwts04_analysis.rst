Analysis of mwts04 experiment
=============================

File locations
--------------

Main experiment directory is 

/home/mark/Documents/research/MultiWeek/pymwts-exps/exps/mwts04/

All of the final output files are in:

/home/mark/Documents/research/MultiWeek/pymwts-exps/exps/mwts04/outputs

The original analysis files are in:

/home/mark/Documents/research/MultiWeek/pymwts-exps/exps/mwts04/mwts04out_original_nonlive

The working analysis files are in Dropbox in:

/home/mark/Dropbox/MWTS/mwts04out

The paper is in /home/mark/Dropbox/MWTS/ and the entire folder is under
version control with git and github.





The key output database is:

-rw-rw-r-- 1 mark mark 4613120 Aug  8  2012 mwts04_d2456_cbc.db


Overall feasibility tests
-------------------------

Many, but not all problems planned to be run, got run.

Of those that were run, most but not all solved. 

According to the output database:

optimal	optimal	3420
optimal	unknown	69

According Excel pivot table supposedly based on same database:

optimal:optimal	3420
optimal:unknown	69
Total Result	3489

There was a problem with a computed sol_status field in the Excel file leading to only
68 optimal:unknown records. After rebuilding the string formula for that column, all is good
and Excel file shows 69 unsolved problems. Now the challenge is to figure out why these
69 problems are phase2 infeasible.

mwts04_d2456_unsolved.csv

Unsolved problems
-----------------

Common features of the unsolved problems - each unsolved problem uses either
mix 6,7 or 8 (or several of these) along with other mixes.

Here are the relevant pieces from mwts04_ash01.mix. Notice that I truncated the
ash part::

      - ttnum: 6
        description: part time 10hr
        
        is_parttime: 1
        tt_lb: 0
        tt_ub: 1000
        tt_cost_multiplier: 1.0
        
       
        shiftlengths: 
          - numbins: 20
            allowable_starttimes: 
              - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...          - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]



            min_shiftlen_days_week: [2,2,2,2]
            max_shiftlen_days_week: [3,3,3,3]
            
            min_shiftlen_cumul_days_week: [2,4,6,8]
            max_shiftlen_cumul_days_week: [3,6,9,12]
        
            min_shiftlen_prds_week: [40,40,40,40]
            max_shiftlen_prds_week: [60,60,60,60]
        
            min_shiftlen_cumul_prds_week: [40,80,120,160]
            max_shiftlen_cumul_prds_week: [60,120,180,240]

      - ttnum: 7
        description: half-time 8hr shifts
        
        is_parttime: 1
        tt_lb: 0
        tt_ub: 1000
        tt_cost_multiplier: 1.0
        
       
        shiftlengths: 
          - numbins: 16
            allowable_starttimes: 
              - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...          - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]



            min_shiftlen_days_week: [2,2,2,2]
            max_shiftlen_days_week: [3,3,3,3]
            
            min_shiftlen_cumul_days_week: [2,5,7,10]
            max_shiftlen_cumul_days_week: [3,5,8,10]
        
            min_shiftlen_prds_week: [32,32,32,32]
            max_shiftlen_prds_week: [48,48,48,48]
        
            min_shiftlen_cumul_prds_week: [32,80,112,160]
            max_shiftlen_cumul_prds_week: [48,80,128,160]

      - ttnum: 8
        description: 3334 12-hr FT
        
        is_parttime: 0
        tt_lb: 0
        tt_ub: 1000
        tt_cost_multiplier: 1.0
        
       
        shiftlengths: 
          - numbins: 24
            allowable_starttimes: 
              - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...          - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]


            min_shiftlen_days_week: [3,3,3,3]
            max_shiftlen_days_week: [4,4,4,4]
            
            min_shiftlen_cumul_days_week: [3,6,9,13]
            max_shiftlen_cumul_days_week: [4,7,10,13]

            min_shiftlen_prds_week: [72,72,72,72]
            max_shiftlen_prds_week: [96,96,96,96]
        
            min_shiftlen_cumul_prds_week: [72,144,216,312]
            max_shiftlen_cumul_prds_week: [96,168,240,312]




Mix 6 debugging
---------------

Mix 6: part-time 10hr shifts (either 2 or 3 days worked per week)

Let's find problems that used mix 6 but not 7 or 8 and were unsolved. What
did the phase1 solution produce that caused no feasible phase2 solution to be
found? Is it a phase1 problem or a phase2 problem?


It's only 

* mwts04_d02_t46_a01_ptub_moderate
* mwts04_d02_t46_a01_ptub_loose
* mwts04_d02_t46_a01_noptub_moderate
* mwts04_d02_t46_a01_noptub_loose



Mix 7 debugging
---------------

Mix 7: 1/2 time 8hr shifts (either a 2-3 or a 3-2)




Mix 8 debugging
---------------

Mix 8: 3334 12-hr FT








































