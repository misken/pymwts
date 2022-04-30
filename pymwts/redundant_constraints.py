# -----------------------------------------------------------------------------
# Almost certainly redundant constraints - binary activators available --------

# Bounds on days worked each week (TTD and TT integration and conservation)

# TODO - actually need to think through possible redundancy of the weekly and cumulative weekly
# constraints on TTD and TT now that the mwdw vars explicitly determine the number of people
# working each day - i.e. the TTD values. Still need these constraints on TTDS_TT since those
# are shift length specific variables.

# Index for both lower and upper bound versions of these constraints are the same

# def TTD_TT_weeklyconservation_idx_rule(M):
#     """
#     Index is (window, tour type, week).
#
#     This index is used in several sets of TTD_TT related constraints.
#
#     :param M: Model
#     :return: Constraint index rule
#     """
#
#     return [(i, t, w) for (i, t) in M.okTourType
#             for w in M.WEEKS]
#
#
# model.TTD_TT_weeklyconservation_idx = pyo.Set(
#     dimen=3, initialize=TTD_TT_weeklyconservation_idx_rule)


# Lower bound on TTD vars based on minimum number of days worked per week
# TODO: These bounds are redundant given TTD_MWDW constraints. Binary deactivation code already added.
# def TTD_TT_weeklyconservation_LB_rule(M, i, t, w):
#     return sum(M.TourTypeDay[i, t, d, w] for d in M.DAYS) >= M.TourType[i, t] * M.tt_min_dys_weeks[t, w]
#
#
# model.TTD_TT_weeklyconservation_LB = \
#     pyo.Constraint(model.TTD_TT_weeklyconservation_idx, rule=TTD_TT_weeklyconservation_LB_rule)


# Upper bound on TTD vars based on maximum number of days worked per week
# TODO: These bounds are redundant given TTD_MWDW constraints Binary deactivation code already added.
# def TTD_TT_weeklyconservation_UB_rule(M, i, t, w):
#     return sum(M.TourTypeDay[i, t, d, w] for d in M.DAYS) <= M.TourType[i, t] * M.tt_max_dys_weeks[t, w]
#
#
# model.TTD_TT_weeklyconservation_UB = \
#     pyo.Constraint(model.TTD_TT_weeklyconservation_idx, rule=TTD_TT_weeklyconservation_UB_rule)


# Cumulative (over weeks) versions of the lower and upper bound constraints immediately above
# TODO: These bounds are redundant given TTD_MWDW constraints
# def TTD_TT_cumul_weeklyconservation_LB_rule(M, i, t, w):
#     return sum(M.TourTypeDay[i, t, d, z] for d in M.DAYS for z in pyo.sequence(w)) >= \
#            M.TourType[i, t] * M.tt_min_cumul_dys_weeks[t, w]
#
#
# model.TTD_TT_cumul_weeklyconservation_LB = \
#     pyo.Constraint(model.TTD_TT_weeklyconservation_idx, rule=TTD_TT_cumul_weeklyconservation_LB_rule)
#
#
# def TTD_TT_cumul_weeklyconservation_UB_rule(M, i, t, w):
#     return sum(M.TourTypeDay[i, t, d, z] for d in M.DAYS for z in pyo.sequence(w)) <= \
#            M.TourType[i, t] * M.tt_max_cumul_dys_weeks[t, w]
#
#
# model.TTD_TT_cumul_weeklyconservation_UB = \
#     pyo.Constraint(model.TTD_TT_weeklyconservation_idx,
#                    rule=TTD_TT_cumul_weeklyconservation_UB_rule)


# The TTDS_TT_weeklyconservation constraints feel redundant given that TTDS summed over lengths is
# equal to TTD and we've already got TTD_TT constraints above.
# TODO - test for redundancy. Binary deactivation code already added.


# def TTDS_TT_weeklyconservation_idx_rule(M):
#     return [(i, t, w) for (i, t) in M.okTourType for w in M.WEEKS]
#
#
# model.TTDS_TT_weeklyconservation_idx = pyo.Set(
#     dimen=3, initialize=TTDS_TT_weeklyconservation_idx_rule)
#
#
# def TTDS_TT_weeklyconservation_LB_rule(M, i, t, w):
#     return sum(M.TourTypeDayShift[i, t, k, d, w]
#                for d in M.DAYS for k in M.tt_length_x[t] if
#                (i, t, k, d, w) in M.TourTypeDayShift_idx) >= \
#            M.TourType[i, t] * M.tt_min_dys_weeks[t, w]
#
#
# model.TTDS_TT_weeklyconservation_LB = \
#     pyo.Constraint(
#         model.TTDS_TT_weeklyconservation_idx,
#         rule=TTDS_TT_weeklyconservation_LB_rule)
#
#
# def TTDS_TT_weeklyconservation_UB_rule(M, i, t, w):
#     return sum(M.TourTypeDayShift[i, t, k, d, w]
#                for d in M.DAYS for k in M.tt_length_x[t] if
#                (i, t, k, d, w) in M.TourTypeDayShift_idx) <= \
#            M.TourType[i, t] * M.tt_max_dys_weeks[t, w]
#
#
# model.TTDS_TT_weeklyconservation_UB = \
#     pyo.Constraint(
#         model.TTDS_TT_weeklyconservation_idx,
#         rule=TTDS_TT_weeklyconservation_UB_rule)
#
#
# def TTDS_TT_cumul_weeklyconservation_LB_rule(M, i, t, w):
#     return sum(M.TourTypeDayShift[i, t, k, d, z]
#                for d in M.DAYS for k in M.tt_length_x[t]
#                for z in pyo.sequence(w)
#                if (i, t, k, d, z) in M.TourTypeDayShift_idx) >= \
#            M.TourType[i, t] * M.tt_min_cumul_dys_weeks[t, w]
#
#
# model.TTDS_TT_cumul_weeklyconservation_LB = \
#     pyo.Constraint(
#         model.TTDS_TT_weeklyconservation_idx,
#         rule=TTDS_TT_cumul_weeklyconservation_LB_rule)
#
#
# def TTDS_TT_cumul_weeklyconservation_UB_rule(M, i, t, w):
#     return sum(M.TourTypeDayShift[i, t, k, d, z]
#                for d in M.DAYS for k in M.tt_length_x[t]
#                for z in pyo.sequence(w)
#                if (i, t, k, d, z) in M.TourTypeDayShift_idx) <= \
#            M.TourType[i, t] * M.tt_max_cumul_dys_weeks[t, w]
#
#
# model.TTDS_TT_cumul_weeklyconservation_UB = \
#     pyo.Constraint(
#         model.TTDS_TT_weeklyconservation_idx,
#         rule=TTDS_TT_cumul_weeklyconservation_UB_rule)

# These should be the shift length specific versions of the above 4 constraints
# TODO: Aren't these redundant given shiftlen versions of num days worked constraints
# like this? period level constraints only seem relevant in an overall sense.
# def prds_worked_shiflen_weekly_idx_rule(M):
#     """
#     Index is (window, tour type, shift length, week).
#
#     :param M: Model
#     :return: Constraint index rule
#     """
#     return [(i, t, k, w) for (i, t) in M.okTourType
#             for k in M.tt_length_x[t]
#             for w in M.WEEKS]
#
#
# model.prds_worked_shiflen_weekly_idx = pyo.Set(
#     dimen=4, initialize=prds_worked_shiflen_weekly_idx_rule)
#
#
# def prds_worked_shiflen_weekly_LB_rule(M, i, t, k, w):
#     """
#     Coordinate TTDS and TT vars for num periods worked each week as
#     part of a given shift length - lower bounds.
#
#     :param M: Model
#     :param i: start window
#     :param t: tour type
#     :param k: shift length
#     :param w: week
#     :return: Constraint rule
#     """
#
#     return sum(M.TourTypeDayShift[i, t, k, d, w] * M.lengths[k]
#                for d in M.DAYS if
#                (i, t, k, d, w) in M.TourTypeDayShift_idx) >= \
#            M.TourType[i, t] * M.tt_shiftlen_min_prds_weeks[t, k, w]
#
#
# model.prds_worked_shiflen_weekly_LB = \
#     pyo.Constraint(
#         model.prds_worked_shiflen_weekly_idx,
#         rule=prds_worked_shiflen_weekly_LB_rule)
#
#
# def prds_worked_shiflen_weekly_UB_rule(M, i, t, k, w):
#     """
#     Coordinate TTDS and TT vars for num periods worked each week as
#     part of a given shift length - upper bounds.
#
#     :param M: Model
#     :param i: start window
#     :param t: tour type
#     :param k: shift length
#     :param w: week
#     :return: Constraint rule
#     """
#
#     return sum(M.TourTypeDayShift[i, t, k, d, w] * M.lengths[k]
#                for d in M.DAYS if
#                (i, t, k, d, w) in M.TourTypeDayShift_idx) <= \
#            M.TourType[i, t] * M.tt_shiftlen_max_prds_weeks[t, k, w]
#
#
# model.prds_worked_shiflen_weekly_UB = \
#     pyo.Constraint(model.prds_worked_shiflen_weekly_idx,
#                    rule=prds_worked_shiflen_weekly_UB_rule)
#
#
# # Cumulative versions of the above 2 constraints
# def prds_worked_cumul_shiflen_weekly_LB_rule(M, i, t, k, w):
#     """
#     Coordinate TTDS and TT vars for num periods worked each week as
#     part of a given shift length - cumulative lower bounds.
#
#     :param M: Model
#     :param i: start window
#     :param t: tour type
#     :param k: shift length
#     :param w: week
#     :return: Constraint rule
#     """
#
#     return sum(M.TourTypeDayShift[i, t, k, d, z] * M.lengths[k]
#                for d in M.DAYS for z in pyo.sequence(w)
#                if (i, t, k, d, z) in M.TourTypeDayShift_idx) >= \
#            M.TourType[i, t] * M.tt_shiftlen_min_cumul_prds_weeks[t, k, w]
#
#
# model.prds_worked_cumul_shiflen_weekly_LB = \
#     pyo.Constraint(model.prds_worked_shiflen_weekly_idx,
#                    rule=prds_worked_cumul_shiflen_weekly_LB_rule)
#
#
# def prds_worked_cumul_shiflen_weekly_UB_rule(M, i, t, k, w):
#     """
#     Coordinate TTDS and TT vars for num periods worked each week as
#     part of a given shift length - cumulative upper bounds.
#
#     :param M: Model
#     :param i: start window
#     :param t: tour type
#     :param k: shift length
#     :param w: week
#     :return: Constraint rule
#     """
#
#     return sum(M.TourTypeDayShift[i, t, k, d, z] * M.lengths[k]
#                for d in M.DAYS for z in pyo.sequence(w)
#                if (i, t, k, d, z) in M.TourTypeDayShift_idx) <= \
#            M.TourType[i, t] * M.tt_shiftlen_max_cumul_prds_weeks[t, k, w]
#
#
# model.prds_worked_cumul_shiflen_weekly_UB = \
#     pyo.Constraint(model.prds_worked_shiflen_weekly_idx,
#                    rule=prds_worked_cumul_shiflen_weekly_UB_rule)
