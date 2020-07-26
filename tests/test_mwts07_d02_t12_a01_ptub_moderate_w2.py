import sys
sys.path.insert(0, "/home/mark/Documents/research/MultiWeek/pymwts")
import pymwts.solvemwts
pymwts.solvemwts.solvemwts("mwts07_d02_t12_a01_ptub_moderate_w2",
                    "./width_inputs/mwts07_d02_t12_a01_ptub_moderate_w2.dat",
                    "./width_outputs/",
                    "gurobi",
                    1200.0,
                    0.02,
                    results_db="../tests/mwts05_d2456_testing.db", force_solve=True,
                    debug_start_windows=True,
                    write_phase1_instance=False,
                    write_phase2_instance=False)


