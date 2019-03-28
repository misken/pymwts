# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:08:59 2012

@author: mark
"""
import yaml
import StringIO
import sqlite3
from os import chdir
from os import getcwd
import itertools
import time
import sys

from pymwtsio.mwts_makedat import mwts_createdat
from pymwtsio.mwts_makedat import csvrow_to_yaml

def mwts_create_yni(fn_yni,
               scenario_name,
               znotes,
               n_prds_per_day,
               n_days_per_week,
               n_weeks,
               filename_dmd,
               filename_min,
               filename_mix,
               filename_wkd,
               filename_ttbounds,
               labor_budget,
               understaff_cost_1,
               understaff_cost_2,
               understaff_1_ub):

    """
    Generate a single YAML yni scenario file for mwts02 experiments.
    """

    D_scenario = {'scenario_name': scenario_name,'znotes':[n for n in znotes]}
    D_time = {'n_prds_per_day': n_prds_per_day, 'n_days_per_week': n_days_per_week, 'n_weeks': n_weeks}
    D_reqd_files = {'filename_dmd': filename_dmd,
                    'filename_min': filename_min,
                    'filename_mix': filename_mix,
                    'filename_wkd': filename_wkd}
                    
    D_opt_files = {'filename_ttbounds': filename_ttbounds}
    
    D_cost = {'labor_budget': labor_budget,
              'understaff_cost_1' : understaff_cost_1,
              'understaff_cost_2': understaff_cost_2,
              'understaff_1_ub': understaff_1_ub}
    
    D_yni = {'scenario': D_scenario, 'time': D_time, 
             'reqd_files': D_reqd_files,
             'opt_files': D_opt_files,
             'cost': D_cost}
    
    fout = open(fn_yni,"w") 
    print yaml.dump(D_yni,fout,default_flow_style=False)
    fout.close()



def mwts_create_yni_files(expt, suffix, db_problemlist, tbl_problemlist):
    """
    Generates a bunch of yni files for the mwts02 experiments. Problem info is read from a sqlite3 db.
    """

    # Connect to the problem list database.
    conn = sqlite3.connect(db_problemlist)
    # Plug in Row so we can access values by column name
    conn.row_factory = sqlite3.Row
    
    # Step through the problem list database and create a yni file for each record.
    cur = conn.cursor()
    cur.execute('select * from ' + tbl_problemlist)
    rows = cur.fetchall()
#    print r.keys()
#    for member in r:
#        print member

    for r in rows:
        scenario_name = r['problem']
        fn_yni = 'exps/' + expt + '/inputs/yni/' + scenario_name + '.yni'
        znotes = ['timestamp','other info']
        n_prds_per_day = 48
        n_days_per_week = 7
        n_weeks = 4
        filename_dmd = r['dmd_file']
        filename_min = r['min_file']
        filename_mix = r['basemix_file']
        filename_wkd = r['wkd_file']
        filename_ttbounds = r['ttbnd_file']
        labor_budget = r['budget']
        understaff_cost_1 = r['us1']
        understaff_cost_2 = r['us2']
        understaff_1_ub = r['usb']
        
        mwts_create_yni(fn_yni,
               scenario_name,
               znotes,
               n_prds_per_day,
               n_days_per_week,
               n_weeks,
               filename_dmd,
               filename_min,
               filename_mix,
               filename_wkd,
               filename_ttbounds,
               labor_budget,
               understaff_cost_1,
               understaff_cost_2,
               understaff_1_ub)

    conn.close()


def mwts_create_dat_files(expt, db_problemlist, tbl_problemlist, pnumlower=1, pnumupper=1000000, maxtocreate=1000000):
    """
    Generates a bunch of dat files for the mwts02 experiments. Problem info is read from a sqlite3 db.
    """

    # Connect to the problem list database.
    conn = sqlite3.connect(db_problemlist)
    # Plug in Row so we can access values by column name
    conn.row_factory = sqlite3.Row
    
    # Step through the problem list database and create a dat file for each record.
    cur = conn.cursor()
    cur.execute('select * from ' + tbl_problemlist + ' where seq>=' + str(pnumlower) + ' and seq<=' +str(pnumupper))
    rows = cur.fetchall()
#    print r.keys()
#    for member in r:
#        print member
    n = 0
    for r in rows:
        scenario_name = r['problem']
        fn_yni = 'exps/' + expt + '/inputs/yni/' + scenario_name + '.yni'
        fn_dat = 'exps/' + expt + '/inputs/dat/' + scenario_name + '.dat'
        
        result = mwts_createdat(fn_yni, fn_dat)
        print scenario_name, result, time.clock() 
        n += 1
        if n == maxtocreate:
            conn.close()
            sys.exit()

    conn.close()

def mwts_create_runpy_files(expt, suffix, db_problemlist, tbl_problemlist, pnumlower=1, pnumupper=1000000, maxtocreate=1000000):
    """
    Generates a bunch of dat files for the mwts04 experiments. Problem info is read from a sqlite3 db.
    """

    # Connect to the problem list database.
    conn = sqlite3.connect(db_problemlist)
    # Plug in Row so we can access values by column name
    conn.row_factory = sqlite3.Row
    
    
    cur = conn.cursor()
    cur.execute('select * from ' + tbl_problemlist + ' where seq>=' + str(pnumlower) + ' and seq<=' +str(pnumupper))
    rows = cur.fetchall()
#    print r.keys()
#    for member in r:
#        print member
    n = 0
    fn_bat = '../' + expt + suffix + '.sh'
    print fn_bat
    f_bat = open(fn_bat,"w")
    for r in rows:
        scenario_name = r['problem']
        print scenario_name
        fn_out = 'exps/' + expt + '/outputs/' + scenario_name + '.out'
        fn_run = 'run_' + scenario_name + '.py'
        f_run = open(fn_run,"w")
        f_run.write('import solvemwts\n')
        solve_cmd = 'solvemwts.solvemwts("' + scenario_name + '",'
        solve_cmd = solve_cmd + '"exps/' + expt + '/inputs/dat/' + scenario_name + '.dat",'
        solve_cmd = solve_cmd + '"exps/' + expt + '/outputs/"' + ','
        solve_cmd = solve_cmd + '"' + r['solver'] + '"' + ','
        solve_cmd = solve_cmd + str(r['timelimit']) + ','
        solve_cmd = solve_cmd + str(r['mipgap']) + ','
        solve_cmd = solve_cmd + 'results_db="' + db_problemlist + '")'
        f_run.write(solve_cmd)
        f_run.close()
        f_bat.write('python ' + fn_run + ' > ' + fn_out + '\n')
        n += 1
        if n == maxtocreate:
            conn.close()
            f_bat.close()
            sys.exit()

    conn.close()
    f_bat.close()


def mwts_create_ttbnd_combos(stub,
               ttnum_ub):

    """
    Generate a YAML tt_bnd for mwts04 experiments.
    """
#range(1,ttnum_ub+1)
    for r in range(1,ttnum_ub+1):
        L=list(itertools.combinations([1,2,3,4,5,6,7,8],r))
        
        for subset in L:
            suffix = ''.join(map(str,subset))
            L_tt = []
            for c in range(1,ttnum_ub+1):
                if c in subset:
                    tt_ub = 1000
                else:
                    tt_ub = 0
                    
                L_tt.append({'ttnum':c,'tt_lb':0,'tt_ub':tt_ub})
               
            D_bnd = {'tourtypes':L_tt}        
            fn_ttbnd = 'exps/mwts04/inputs/mix/'+stub+suffix+'.bnd'        
            
            fout = open(fn_ttbnd,"w") 
            print yaml.dump(D_bnd,fout,default_flow_style=False)
            fout.close()

def mwts_create_main_combos(
               dmd_ub,
               ttnum_ub,
               wkends,
               ash_ub):

    """
    Generate a YAML tt_bnd for mwts04 experiments.
    """
    fn_mainlist = 'exps/mwts04/inputs/mainlist.txt'
    fout = open(fn_mainlist,"w")
    
    tt_list = []
    for r in range(1,ttnum_ub+1):
        L=list(itertools.combinations([1,2,3,4,5,6,7,8],r))
            
        for subset in L:
            suffix = ''.join(map(str,subset))
            tt_list.append(suffix)
    
    for d in range(1,dmd_ub+1):
        for t in tt_list:
            for w in wkends:
                for a in range(1,ash_ub+1):
                    lineout = '{} {} {} {}'.format(d,''.join(t),w,a) + '\n'
                 
                    fout.write(lineout)
    fout.close()
#========================================================================================






def main():
    
#    mwtsdir = getcwd()
#    print mwtsdir
#    chdir('/home/mark/Documents/research/MultiWeek/pymwts/pymwts/exps/mwts02/')
#    

   
#    ashdir = '/home/mark/Documents/research/MultiWeek/pymwts/pymwts/exps/mwts04/inputs/ash/'
#    csv = ashdir+'ash_all.csv'
#    yml = ashdir+'ash_all.yaml'
#        
#    ash = StringIO.StringIO()  
#    print >>ash, csvrow_to_yaml(csv)
#    fout = open(yml,"w")
#    print >>fout, ash.getvalue()
#    fout.close()
    
    #mwts_create_ttbnd_combos("mwts04_",8)  
    #mwts_create_main_combos(6,8,['tight','moderate','loose'],2)
    
    maxtocreate = 500
    for num in range(1500,12500,maxtocreate):
        suffix = '_' + str(num) + '_' + str(num + maxtocreate)
        mwts_create_runpy_files('mwts04',suffix,'exps/mwts04/mwts04_d2456_cbc.db','problem_list',num,num + maxtocreate,1000000)
     
    
    
    
    
    
    
#    mwts_create_dat_files('mwts04','exps/mwts04/mwts04_d2456.db','problem_list',1500,2000,maxtocreate)

#    mwts_createdat('exps/mwts04/inputs/yni/mwts04_d02_t12345678_a01_ptub_moderate.yni', 
#                   'exps/mwts04/inputs/dat/mwts04_d02_t12345678_a01_ptub_moderate.dat')
#    
#    chdir(mwtsdir) 
#    print getcwd()
#    mwts_create_runpy_files('mwts04', 
#                             'exps/mwts04/mwts04_d2456.db', 
#
#                                   'problem_list',maxtocreate)

if __name__ == '__main__':
    main()
