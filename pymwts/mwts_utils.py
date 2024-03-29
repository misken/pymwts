"""
Shared model components and functions
"""

# Author: misken
# License: TBD

import io

# import pyomo
# import numpy as np
#
# import mwts_shared
from pymwts.pymwtsio.mwts_makedat import scalar_to_param, list_to_param


def shift_to_param(param_name, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Phase 1 Shift variable to a GMPL representation of a parameter.

    Used by `solvemwts` in Phase 2 dat file creation.

    :param param_name: name (str) of paramter in GMPL file
    :param inst: Model instance
    :param reverseidx: True to reverse the order of the indexes (essentially transposing the matrix)
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: GMPL dat code for list parameter either as a StringIO
        object or a string.
    """

    param = 'param ' + param_name + ' default 0 :=\n'
    for (i, j, w, k, t) in inst.okShifts:
        try:
            val = int(round(inst.Shift[i, j, w, k, t]()))
        except TypeError:
            val = inst.Shift[i, j, w, k, t]()

        if val > 0:
            pos_list = [str(p) for p in (i, j, w, k, t)]

            if reverseidx:
                pos_list.reverse()

            data_row = ' '.join(pos_list) + ' ' + str(val) + '\n'
            param += data_row

    param += ";\n"
    if isStringIO:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def tourtype_to_param(param_name, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Phase 1 TourType variable to a GMPL representation of a parameter.

    Used by `solvemwts` in Phase 2 dat file creation.

    :param param_name: name (str) of parameter in GMPL file
    :param inst: Model instance
    :param reverseidx: True to reverse the order of the indexes (essentially transposing the matrix)
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: GMPL dat code for list parameter either as a StringIO
        object or a string.
    """

    param = 'param ' + param_name + ' default 0 :=\n'
    for (i, t) in inst.TourType_idx:
        try:
            val = int(round(inst.TourType[i, t]()))
        except TypeError:
            val = inst.TourType[i, t]()

        if val > 0:
            pos_list = [str(p) for p in (i, t)]

            if reverseidx:
                pos_list.reverse()

            data_row = ' '.join(pos_list) + ' ' + str(val) + '\n'
            param += data_row

    param += ";\n"
    if isStringIO:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def tourtypeday_to_param(param_name, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Phase 1 TourTypeDay variable to a GMPL representation of a parameter.

    Used by `solvemwts` in Phase 2 dat file creation.

    :param param_name: name (str) of parameter in GMPL file
    :param inst: Model instance
    :param reverseidx: True to reverse the order of the indexes (essentially transposing the matrix)
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: GMPL dat code for list parameter either as a StringIO
        object or a string.
    """

    param = 'param ' + param_name + ' default 0 :=\n'
    for (i, t, j, w) in inst.TourTypeDay_idx:
        try:
            val = int(round(inst.TourTypeDay[i, t, j, w].value))
        except TypeError:
            val = inst.TourTypeDay[i, t, j, w].value

        if val > 0:
            pos_list = [str(p) for p in (i, t, j, w)]

            if reverseidx:
                pos_list.reverse()

            data_row = ' '.join(pos_list) + ' ' + str(val) + '\n'
            param += data_row

    param += ";\n"
    if isStringIO:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def tourtypedayshift_to_param(param_name, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Phase 1 TourTypeDayShift variable to a GMPL representation of a parameter.

    Used by `solvemwts` in Phase 2 dat file creation.

    :param param_name: name (str) of parameter in GMPL file
    :param inst: Model instance
    :param reverseidx: True to reverse the order of the indexes (essentially transposing the matrix)
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: GMPL dat code for list parameter either as a StringIO
        object or a string.
    """

    param = 'param ' + param_name + ' default 0 :=\n'
    for (i, t, k, j, w) in inst.TourTypeDayShift_idx:
        try:
            val = int(round(inst.TourTypeDayShift[i, t, k, j, w]()))
        except TypeError:
            val = inst.TourTypeDayShift[i, t, k, j, w]()

        if val > 0:
            pos_list = [str(p) for p in (i, t, k, j, w)]

            if reverseidx:
                pos_list.reverse()

            data_row = ' '.join(pos_list) + ' ' + str(val) + '\n'
            param += data_row

    param += ";\n"
    if isStringIO:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def weekenddaysworked_to_param(param_name, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Phase 1 WeekendDaysWorked variable to a GMPL representation of a parameter.

    Used by `solvemwts` in Phase 2 dat file creation.

    :param param_name: name (str) of parameter in GMPL file
    :param inst: Model instance
    :param reverseidx: True to reverse the order of the indexes (essentially transposing the matrix)
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: GMPL dat code for list parameter either as a StringIO
        object or a string.
    """

    param = 'param ' + param_name + ' default 0 :=\n'
    for (i, t, d) in inst.weekend_days_worked_idx:
        try:
            val = int(round(inst.WeekendDaysWorked[i, t, d]()))
        except TypeError:
            val = inst.WeekendDaysWorked[i, t, d]()

        if val > 0:
            pos_list = [str(p) for p in (i, t, d)]

            if reverseidx:
                pos_list.reverse()

            data_row = ' '.join(pos_list) + ' ' + str(val) + '\n'
            param += data_row

    param += ";\n"
    if isStringIO:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def multiweekdaysworked_to_param(param_name, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Phase 1 MWDW variable to a GMPL representation of a parameter.

    Used by `solvemwts` in Phase 2 dat file creation.

    :param param_name: name (str) of parameter in GMPL file
    :param inst: Model instance
    :param reverseidx: True to reverse the order of the indexes (essentially transposing the matrix)
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: GMPL dat code for list parameter either as a StringIO
        object or a string.
    """

    param = 'param ' + param_name + ' default 0 :=\n'
    for (i, t, p1, p2) in inst.multiweekdaysworked_idx:
        try:
            val = int(round(inst.MultiWeekDaysWorked[i, t, p1, p2]()))
        except TypeError:
            val = inst.MultiWeekDaysWorked[i, t, p1, p2]()

        if val > 0:
            pos_list = [str(p) for p in (i, t, p1, p2)]

            if reverseidx:
                pos_list.reverse()

            data_row = ' '.join(pos_list) + ' ' + str(val) + '\n'
            param += data_row

    param += ";\n"
    if isStringIO:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param





def tour_WIN_TT_to_param(inst, isStringIO=True):
    """
    Convert WIN_TT values to a GMPL representation of a parameter.

    Used by `solvemwts` in Phase 2 dat file creation.

    :param inst: Model instance
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: GMPL dat code for list parameter either as a StringIO
        object or a string.
    """

    n_tours = int(round((sum(inst.TourType[i, t]() for (i, t) in inst.okTourType))))
    tour_num = 0

    WIN_x = [0 for _ in range(n_tours + 1)]
    TT_x = [0 for _ in range(n_tours + 1)]
    for (i, t) in inst.okTourType:
        try:
            val = int(round(inst.TourType[i, t]()))
        except TypeError:
            val = inst.TourType[i, t]

        if val > 0:
            tour_num_lower = tour_num + 1
            tour_num_upper = tour_num + val
            tour_num = tour_num_upper

            for idx in range(tour_num_lower, tour_num_upper + 1):
                WIN_x[idx] = i
                TT_x[idx] = t

    param_name = 'WIN_x'
    param = 'param ' + param_name + ':=\n'
    for i in range(1, n_tours + 1):
        data_row = str(i) + ' ' + str(WIN_x[i]) + '\n'
        param += data_row

    param += ";\n"

    param_name = 'TT_x'
    param += 'param ' + param_name + ':=\n'
    for i in range(1, n_tours + 1):
        data_row = str(i) + ' ' + str(TT_x[i]) + '\n'
        param += data_row

    param += ";\n"

    if isStringIO:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def logger(f, msg, ts):
    #    print msg, ts
    msgts = '{}: {}\n'.format(msg, ts)
    f.write(msgts)



