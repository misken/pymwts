
import solvemwts
solvemwts.solvemwts("mwts05_d02_t3_a01_ptub_tight",
                    "../tests/inputs/mwts05_d02_t3_a01_ptub_tight.dat",
                    "../tests/outputs/"
                    ,"gurobi",
                    1200.0,
                    0.02,
                    results_db="../tests/mwts05_d2456_testing.db")