# -------------------------------------------------------------------------------
# Name:        pymwts_utils.py
# Purpose:
#
# Author:      isken
#
# Created:     03/13/2012
# Copyright:   (c) isken 2012
# Licence:     <your licence>
# -------------------------------------------------------------------------------

import io
import pyomo
import numpy as np


# The following lives in mwts_phase1 and mwts_phase2 models project. Need to figure out how to use it from there.
# Just copying it for now so I can get this model working.

def n_prds_per_week_init(M):
    return M.n_days_per_week()*M.n_prds_per_day()


def n_prds_per_cycle_init(M):
    return M.n_weeks()*M.n_days_per_week()*M.n_prds_per_day()


def g_period_init(M, i, j, w):
    return ((w - 1) * M.n_days_per_week() * M.n_prds_per_day() +
            (j - 1) * M.n_prds_per_day() + i)


def g_tuple_to_period(i, j, w, n_prds_per_day, n_days_per_week):
    return ((w - 1) * n_days_per_week * n_prds_per_day +
            (j - 1) * n_prds_per_day + i)


def period_increment(M, i, j, w, incr):
    p = M.g_period[i, j, w]
    if p + incr <= M.n_prds_per_cycle:
        return p + incr
    else:
        return p + incr - M.n_prds_per_cycle


def g_period_increment(M, p, incr):
    if p + incr <= M.n_prds_per_cycle:
        return p + incr
    else:
        return p + incr - M.n_prds_per_cycle


def g_period_difference(M, b_prd, e_prd):
    if e_prd >= b_prd:
        return e_prd - b_prd + 1
    else:
        return M.n_prds_per_cycle + e_prd - b_prd + 1


def g_prd_to_tuple(M, p):
    #    param which_prd{p in 1..(n_days+1)*n_prds_per_day} :=
    #   p-n_prds_per_day*(ceil(p/n_prds_per_day-1));
    #
    # param which_day{p in 1..(n_days+1)*n_prds_per_day} :=
    #   (if p>n_prds_per_day*n_days then 1 else 1+ceil(p/n_prds_per_day-1));

    n_week = ((p - 1) // M.n_prds_per_week.value) + 1
    prds_remainder = p - (n_week - 1) * M.n_prds_per_week
    if prds_remainder == 0:
        n_day = 1
    else:
        n_day = ((prds_remainder - 1) // M.n_prds_per_day.value) + 1

    prds_remainder = prds_remainder - (n_day - 1) * M.n_prds_per_day
    if prds_remainder == 0:
        n_period = 1
    else:
        n_period = prds_remainder

    return n_period, n_day, n_week



# The following lives in  pymtwsio project. Need to figure out how to use it from there.
# Just copying it for now so I can get this model working.
def scalar_to_param(pname, pvalue, isStringIO=True):
    """
    Convert a scalar to a GMPL representation of a parameter.
    Inputs:
        pname - string name of paramter in GMPL file
        pvalue - value of parameter
        isStringIO - true to return StringIO object, false to return string

    Output:
        GMPL dat code for scalar parameter either as a StringIO
        object or a string.

        Example:

            param n_prds_per_day :=  48; 
    """

    param = 'param ' + pname + ' :=  ' + str(pvalue) + ';\n'
    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def list_to_param(pname, plist, reverseidx=False, isStringIO=True):
    """
    Convert a list to a GMPL representation of a parameter.
    Inputs:
        pname - string name of paramter in GMPL file
        plist - list containing parameter (could be N-Dimen list)
        reverseidx - True to reverse the order of the indexes (essentially transposing the matrix)
        isStringIO - True to return StringIO object, False to return string

    Output:
        GMPL dat code for list parameter either as a StringIO
        object or a string.

        Example:
            param midnight_thresh:=
             1 100
             2 100
             3 100
            ; 

    """
    # Convert parameter as list to an ndarray
    parray = np.array(plist)

    # Denumerate the array to get at the index tuple and array value
    paramrows = np.ndenumerate(parray)
    param = 'param ' + pname + ':=\n'
    for pos, val in paramrows:
        poslist = [str(p + 1) for p in pos]
        if reverseidx:
            poslist.reverse()
        datarow = ' '.join(poslist) + ' ' + str(val) + '\n'
        param += datarow

    param += ";\n"
    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def shift_to_param(pname, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Pyomo indexed variable to a GMPL representation of a parameter.
    Inputs:
        pname - string name of paramter in GMPL file
        plist - list containing parameter (could be N-Dimen list)
        reverseidx - True to reverse the order of the indexes (essentially transposing the matrix)
        isStringIO - True to return StringIO object, False to return string

    Output:
        GMPL dat code for list parameter either as a StringIO
        object or a string.

        Example:
            param midnight_thresh:=
             1 100
             2 100
             3 100
            ; 

    """
    # Convert parameter as list to an ndarray

    # Denumerate the array to get at the index tuple and array value

    param = 'param ' + pname + ' default 0 :=\n'
    for (i, j, w, k, t) in inst.okShifts:
        try:
            val = int(round(inst.Shift[i, j, w, k, t]()))
        except:
            val = inst.Shift[i, j, w, k, t]()

        if val > 0:
            poslist = [str(p) for p in (i, j, w, k, t)]

            if reverseidx:
                poslist.reverse()

            datarow = ' '.join(poslist) + ' ' + str(val) + '\n'
            param += datarow

    param += ";\n"
    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def tourtype_to_param(pname, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Pyomo indexed variable to a GMPL representation of a parameter.
    Inputs:
        pname - string name of paramter in GMPL file
        plist - list containing parameter (could be N-Dimen list)
        reverseidx - True to reverse the order of the indexes (essentially transposing the matrix)
        isStringIO - True to return StringIO object, False to return string

    Output:
        GMPL dat code for list parameter either as a StringIO
        object or a string.

        Example:
            param midnight_thresh:=
             1 100
             2 100
             3 100
            ; 

    """
    # Convert parameter as list to an ndarray

    # Denumerate the array to get at the index tuple and array value

    param = 'param ' + pname + ' default 0 :=\n'
    for (i, t) in inst.TourType_idx:
        try:
            val = int(round(inst.TourType[i, t]()))
        except:
            val = inst.TourType[i, t]()

        if val > 0:
            poslist = [str(p) for p in (i, t)]

            if reverseidx:
                poslist.reverse()

            datarow = ' '.join(poslist) + ' ' + str(val) + '\n'
            param += datarow

    param += ";\n"
    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def dailytourtype_to_param(pname, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Pyomo indexed variable to a GMPL representation of a parameter.
    Inputs:
        pname - string name of paramter in GMPL file
        plist - list containing parameter (could be N-Dimen list)
        reverseidx - True to reverse the order of the indexes (essentially transposing the matrix)
        isStringIO - True to return StringIO object, False to return string

    Output:
        GMPL dat code for list parameter either as a StringIO
        object or a string.

        Example:
            param midnight_thresh:=
             1 100
             2 100
             3 100
            ; 

    """
    # Convert parameter as list to an ndarray

    # Denumerate the array to get at the index tuple and array value

    param = 'param ' + pname + ' default 0 :=\n'
    for (i, t, j, w) in inst.DailyTourType_idx:
        try:
            val = int(round(inst.DailyTourType[i, t, j, w].value))
        except:
            val = inst.DailyTourType[i, t, j, w].value

        if val > 0:
            poslist = [str(p) for p in (i, t, j, w)]

            if reverseidx:
                poslist.reverse()

            datarow = ' '.join(poslist) + ' ' + str(val) + '\n'
            param += datarow

    param += ";\n"
    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def dailyshiftworked_to_param(pname, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Pyomo indexed variable to a GMPL representation of a parameter.
    Inputs:
        pname - string name of paramter in GMPL file
        plist - list containing parameter (could be N-Dimen list)
        reverseidx - True to reverse the order of the indexes (essentially transposing the matrix)
        isStringIO - True to return StringIO object, False to return string

    Output:
        GMPL dat code for list parameter either as a StringIO
        object or a string.

        Example:
            param midnight_thresh:=
             1 100
             2 100
             3 100
            ; 

    """
    # Convert parameter as list to an ndarray

    # Denumerate the array to get at the index tuple and array value

    param = 'param ' + pname + ' default 0 :=\n'
    for (i, t, k, j, w) in inst.DailyShiftWorked_idx:
        try:
            val = int(round(inst.DailyShiftWorked[i, t, k, j, w]()))
        except:
            val = inst.DailyShiftWorked[i, t, k, j, w]()

        if val > 0:
            poslist = [str(p) for p in (i, t, k, j, w)]

            if reverseidx:
                poslist.reverse()

            datarow = ' '.join(poslist) + ' ' + str(val) + '\n'
            param += datarow

    param += ";\n"
    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def weekenddaysworked_to_param(pname, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Pyomo indexed variable to a GMPL representation of a parameter.
    Inputs:
        pname - string name of paramter in GMPL file
        plist - list containing parameter (could be N-Dimen list)
        reverseidx - True to reverse the order of the indexes (essentially transposing the matrix)
        isStringIO - True to return StringIO object, False to return string

    Output:
        GMPL dat code for list parameter either as a StringIO
        object or a string.

        Example:
            param midnight_thresh:=
             1 100
             2 100
             3 100
            ; 

    """
    # Convert parameter as list to an ndarray

    # Denumerate the array to get at the index tuple and array value

    param = 'param ' + pname + ' default 0 :=\n'
    for (i, t, d) in inst.weekenddaysworked_idx:
        try:
            val = int(round(inst.WeekendDaysWorked[i, t, d]()))
        except:
            val = inst.WeekendDaysWorked[i, t, d]()

        if val > 0:
            poslist = [str(p) for p in (i, t, d)]

            if reverseidx:
                poslist.reverse()

            datarow = ' '.join(poslist) + ' ' + str(val) + '\n'
            param += datarow

    param += ";\n"
    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def multiweekdaysworked_to_param(pname, inst, reverseidx=False, isStringIO=True):
    """
    Convert a Pyomo indexed variable to a GMPL representation of a parameter.
    Inputs:
        pname - string name of paramter in GMPL file
        plist - list containing parameter (could be N-Dimen list)
        reverseidx - True to reverse the order of the indexes (essentially transposing the matrix)
        isStringIO - True to return StringIO object, False to return string

    Output:
        GMPL dat code for list parameter either as a StringIO
        object or a string.

        Example:
            param midnight_thresh:=
             1 100
             2 100
             3 100
            ;

    """
    # Convert parameter as list to an ndarray

    # Denumerate the array to get at the index tuple and array value

    param = 'param ' + pname + ' default 0 :=\n'
    for (i, t, p1, p2) in inst.multiweekpattern_idx:
        try:
            val = int(round(inst.MultiWeekPattern[i, t, p1, p2]()))
        except:
            val = inst.MultiWeekPattern[i, t, p1, p2]()

        if val > 0:
            poslist = [str(p) for p in (i, t, p1, p2)]

            if reverseidx:
                poslist.reverse()

            datarow = ' '.join(poslist) + ' ' + str(val) + '\n'
            param += datarow

    param += ";\n"
    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def weekenddaysworked_to_tourskeleton(inst, isStringIO=True):
    """
    

    """

    # build the header
    headerlist = ['n', 'i', 't', 'p']
    for w in range(1, inst.n_weeks + 1):
        for d in range(1, inst.n_days_per_week + 1):
            headerlist.append(str(w) + '_' + str(d))

    param = ','.join(headerlist) + '\n'

    tnum = 0
    for (i, t, pattern) in inst.weekenddaysworked_idx:
        try:
            val = int(round(inst.WeekendDaysWorked[i, t, pattern]()))
        except:
            val = inst.WeekendDaysWorked[i, t, pattern]()

        if val > 0:
            poslist = [str(p) for p in (i, t, pattern)]

            for tour in range(val):
                tnum += 1

                datarow = str(tnum) + ',' + ','.join(poslist) + ','
                e = inst.weekend_type[i, t]
                daylist = []
                for w in range(1, inst.n_weeks + 1):
                    for d in range(1, inst.n_days_per_week + 1):
                        daylist.append(str(inst.A[pattern, d, w, t, e]))

                datarow += ','.join(daylist) + '\n'
                param += datarow

    param += '\n'

    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def dailytourtype_to_tourskeleton(inst, isStringIO=True):
    """
    

    """

    # build the header
    headerlist = ['n', 'i', 't', 'dummy']
    for w in range(1, inst.n_weeks + 1):
        for d in range(1, inst.n_days_per_week + 1):
            headerlist.append(str(w) + '_' + str(d))

    param = ','.join(headerlist) + '\n'

    tnum = 0
    for (i, t) in inst.okTourType:
        if inst.TourType[i, t].value > 0:
            poslist = [str(p) for p in (i, t)]
            tnum += 1
            datarow = str(tnum) + ',' + ','.join(poslist) + ', x, '
            daylist = []
            for w in range(1, inst.n_weeks + 1):
                for d in range(1, inst.n_days_per_week + 1):
                    daylist.append(str(inst.DailyTourType[i, t, d, w]()))

            datarow += ','.join(daylist) + '\n'
            param += datarow

    param += "\n"

    # build the header

    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def tour_WIN_TT_to_param(inst, isStringIO=True):
    """
    Convert certain Pyomo indexed variables to a GMPL representation of a parameter.
    
    Operates on WIN and TT variables.
    
    Inputs:
        inst - string name of paramter in GMPL file
        isStringIO - True to return StringIO object, False to return string

    Output:
        GMPL dat code for list parameter either as a StringIO
        object or a string.

        Example:
            param midnight_thresh:=
             1 100
             2 100
             3 100
            ; 

    """

    n_tours = int(round((sum(inst.TourType[i, t]() for (i, t) in inst.okTourType))))
    tour_num = 0

    WIN_x = [0 for j in range(n_tours + 1)]
    TT_x = [0 for j in range(n_tours + 1)]
    for (i, t) in inst.okTourType:
        try:
            val = int(round(inst.TourType[i, t]()))
        except:
            val = inst.TourType[i, t]

        if val > 0:
            tnum_lower = tour_num + 1
            tnum_upper = tour_num + val
            # print tnum_lower, tnum_upper
            # touridxs = range(tour_num+1,tour_num+int(M.TourType[i,t]))
            tour_num = tnum_upper

            for idx in range(tnum_lower, tnum_upper + 1):
                WIN_x[idx] = i
                TT_x[idx] = t

                # Denumerate the array to get at the index tuple and array value

    pname = 'WIN_x'
    param = 'param ' + pname + ':=\n'
    for i in range(1, n_tours + 1):
        datarow = str(i) + ' ' + str(WIN_x[i]) + '\n'
        param += datarow

    param += ";\n"

    pname = 'TT_x'
    param += 'param ' + pname + ':=\n'
    for i in range(1, n_tours + 1):
        datarow = str(i) + ' ' + str(TT_x[i]) + '\n'
        param += datarow

    param += ";\n"

    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def logger(f, msg, ts):
    #    print msg, ts
    msgts = '{}: {}\n'.format(msg, ts)
    f.write(msgts)


def write_phase1_shiftsummary(inst, isStringIO=True):
    """
    Write out a multiweek summary of daily shift worked variables 
    in multi-week format. Hopefully useful for debugging.
    
    Inputs:
        pname - string name of paramter in GMPL file
        plist - list containing parameter (could be N-Dimen list)
        isStringIO - True to return StringIO object, False to return string

    Output:
        Shift summary as StringIO object or a string.

        Example:
            TODO

    """

    param = ''
    for t in inst.activeTT:
        for i in inst.PERIODS:
            if (i, t) in inst.okTourType:
                for k in inst.tt_length_x[t]:
                    datarow = '{}|{}|{}|{}|{}'.format(t, i, k, inst.TourType[i, t], inst.TourType[i, t].value)
                    for w in inst.WEEKS:
                        for j in inst.DAYS:
                            if (i, t, k, j, w) in inst.DailyShiftWorked_idx:
                                datarow += '|{}|{}'.format(inst.DailyShiftWorked[i, t, k, j, w], inst.DailyShiftWorked[i, t, k, j, w].value)
                            else:
                                datarow += '|{}|{}'.format('No shift', 0)
                    if inst.TourType[i, t].value > 0:
                        datarow += '\n'
                        param += datarow

    if isStringIO:
        paramout = io.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param
