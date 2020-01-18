"""
Shared model components and functions
"""

# Author: misken
# License: TBD


def n_prds_per_week_init(M):
    """
    Initialize convenience parameter n_prds_per_week
    """
    return M.n_days_per_week() * M.n_prds_per_day()


def n_prds_per_cycle_init(M):
    """
    Initialize convenience parameter n_prds_per_cycle where cycle may include
    one or more weeks.
    """
    return M.n_weeks() * M.n_days_per_week() * M.n_prds_per_day()


def epoch_init(M, i, j, w):
    """
    Initialize epoch index from daily period, day of week and week of cycle.

    :param M: Model
    :param i: period of day
    :param j: day of week
    :param w: week of cycle
    :return: period of cycle in 1..n_prds_per_cycle
    """
    return (w - 1) * M.n_days_per_week() * M.n_prds_per_day() + (j - 1) * M.n_prds_per_day() + i


def epoch_increment(M, p, incr):
    """
    Increment period of cycle, circling as necessary

    :param M: Model
    :param p: period
    :param incr: Number of periods to increment by
    :return: incremented period
    """
    if p + incr <= M.n_prds_per_cycle:
        return p + incr
    else:
        return p + incr - M.n_prds_per_cycle


def epoch_difference(M, b_prd, e_prd):
    """
    Compute difference between two epochs accounting for end of horizon circling

    :param M: Model
    :param b_prd: global start period
    :param e_prd: global end period
    :return: end - start + 1 (+ n_prds_per_cycle if spans end of horizon)
    """

    if e_prd >= b_prd:
        return e_prd - b_prd + 1
    else:
        return M.n_prds_per_cycle + e_prd - b_prd + 1


def epoch_to_tuple(M, p):
    """
    Convert epoch to (daily period, day, week) tuple

    :param M: Model
    :param p: epoch
    :return: (daily period, day, week) tuple
    """

    week = ((p - 1) // M.n_prds_per_week.value) + 1
    prds_remainder = p - (week - 1) * M.n_prds_per_week
    if prds_remainder == 0:
        day = 1
    else:
        day = ((prds_remainder - 1) // M.n_prds_per_day.value) + 1

    prds_remainder = prds_remainder - (day - 1) * M.n_prds_per_day
    if prds_remainder == 0:
        period = 1
    else:
        period = prds_remainder

    return period, day, week
