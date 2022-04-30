import sys
sys.path.insert(0, "/home/mark/Documents/research/MultiWeek/pymwts")
import pymwts.solvemwts
pymwts.solvemwts.solvemwts("scenario3a_tt123",
                    "./width_inputs/scenario3a_tt123.dat",
                    "./width_outputs/",
                    "gurobi",
                    100.0,
                    0.02,
                    results_db="../tests/mwts05_d2456_testing.db", force_solve=True,
                    debug_start_windows=True,
                    write_phase1_instance=True,
                    write_phase2_instance=False)


