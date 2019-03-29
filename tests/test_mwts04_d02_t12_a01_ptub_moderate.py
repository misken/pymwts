import solvemwts
solvemwts.solvemwts("mwts04_d02_t12_a01_ptub_moderate",
                    "../tests/inputs/mwts04_d02_t12_a01_ptub_moderate.dat",
                    "../tests/outputs/",
                    "cbc",
                    360.0,
                    0.05,
                    results_db="../tests/mwts04_d2456_cbc_testing.db")
