import sys
sys.path.insert(0,'/home/mark/Documents/research/MultiWeek/pymwts/pymwts')

import solvemwts
solvemwts.solvemwts("mwts04_d02_t168_a02_noptub_loose",
                    "../tests/inputs/mwts04_d02_t168_a02_noptub_loose.dat",
                    "../tests/outputs/",
                    "gurobi",
                    1200.0,
                    0.02,
                    results_db="../tests/mwts04_d2456_cbc_testing.db")


