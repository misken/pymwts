import sys
sys.path.insert(0,'/home/mark/Documents/research/MultiWeek/pymwts/pymwts')

import solvemwts
solvemwts.solvemwts("test_mwts04_d02_t27_a01_ptub_moderate",
                    "../tests/inputs/mwts04_d02_t27_a01_ptub_moderate.dat",
                    "../tests/outputs/",
                    "cbc",
                    1200.0,
                    0.02,
                    results_db="../tests/mwts04_d2456_cbc_testing.db")


