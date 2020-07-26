import sys
sys.path.insert(0, "/home/mark/Documents/research/MultiWeek/pymwts")
import pymwts.solvemwts
pymwts.solvemwts.solvemwts("mwts07_d02_t1_a01_noptub_moderate",
                    "./width_inputs/mwts07_d02_t1_a01_noptub_moderate.dat",
                    "./width_outputs/",
                    "cbc",
                    100.0,
                    0.02,
                    results_db="../tests/mwts05_d2456_testing.db",
                    force_solve=True,
                    write_phase1_instance=True)


