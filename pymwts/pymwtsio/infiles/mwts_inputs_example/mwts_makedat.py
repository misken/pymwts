#-------------------------------------------------------------------------------
# Name:        mwts_makedat
# Purpose:     mwts dat file creation
#
# Author:      isken
#
# Created:     28/09/2011
# Copyright:   (c) isken 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import sys
import StringIO
import yaml
from numpy import *
import csv
import json
import itertools


"""
mwts_makedat is a module for reading input files for mwts problems
and creating an AMPL/GMPL data file.
It is a replacement for the ancient createssdat.c program that was
used to create AMPL/GMPL dat files for one week tour scheduling problems.

"""



def create_weekend_base(n_weeks):
    """
    Generate basis for cartesion product of [0,1] lists based
    on number of weeks in scheduling problem. Each list
    element is one week. The tuple of binary values represent the
    first and second day of the weekend for that week. A 1 means
    the day is worked, a 0 means it is off.

    Input:
        n_weeks - number of weeks in scheduling horizon

    Output:
        Result is all the possible n_weeks weekends worked patterns.
        Example: n_weeks = 4 --> 256 possible weekends worked patterns. This
        exhaustive list can later be filtered to only include desirable patterns.

    """


    basis_list = [[0,0],[1,0],[0,1],[1,1]]

    mw_basis_list = []
    for i in range(n_weeks):
        mw_basis_list.append(basis_list)

    # Use itertools to create the n_weeks cartesion product of the basis_list.

    return list(itertools.product(*mw_basis_list))

def filterpatterns(pattern,ttnum,wkendtype,ttspec):
    """
    Creates a sequence of binary values to be used for list filtering. This
    function will contain the various rules used to filter out weekend days
    worked patterns that we don't want to allow.

    For now I'm hard coding in rules but need to develop an approach to
    flexibly specifiying rules to apply to filter out undesirable weekends
    worked patterns.

        Inputs:
          x - list of 2-tuples representing weekend days worked. Each list
            element is one week. The tuple of binary values represent the
            first and second day of the weekend for that week. A 1 means
            the day is worked, a 0 means it is off.
          type - 1 --> weekend consists of Saturday and Sunday
                 2 --> weekend consists of Friday and Saturday


          max_days_worked - max # of weekend days worked over horizon
          max_wkends_worked - max # of weekends in which >= 1 day worked
          half_weekends_ok - True or False
          max_consec_wkends - max consecutive weeks with >= 1 day worked



        Examples:
          (1) Type 1, work every other weekend
            pattern = [(0,1),(1,0),(0,1),(1,0)], type = 1

          (2) Type 1, work every other weekend
            pattern = [(1,1),(0,0),(1,1),(0,0)], type = 2

        Output: True --> keep pattern
                False --> discard pattern

    """
    n_weeks = len(pattern)
    keep = True

    tourtype = [t for t in ttspec['tourtypes'] if t['ttnum'] == ttnum]

    # No more than max_days_worked over the scheduling horizon
    max_days_worked = tourtype[0]['max_days_worked']
    if not (sum(pattern) <= max_days_worked):
        keep = False 

    # No consecutive weekends with one or more days worked
    window = ntuples(pattern,2)
    for pair in window:
        if sum(pair) > 2:
            keep = False

    # No half-weekends
    if not tourtype[0]['half_weekends_ok'] and num_half_weekends(pattern,wkendtype) > 0:
        keep = False

    return keep



def num_full_weekends(x,wkendtype):
    """
    Returns number of full weekends (both days) worked in a given weekends worked pattern.

    Inputs:
          x - list of 2-tuples representing weekend days worked. Each list
            element is one week. The tuple of binary values represent the
            first and second day of the weekend for that week. A 1 means
            the day is worked, a 0 means it is off.

          wkendtype - 1 --> weekend consists of Saturday and Sunday
                 2 --> weekend consists of Friday and Saturday

    Output:
        Number of full weekends worked

    Example:
        n = num_full_weekends([(0,1),(1,0),(0,1),(1,0)],1)
        # n = 2

        n = num_full_weekends([(0,1),(1,0),(0,1),(0,0)],1)
        # n = 1

        n = num_full_weekends([(1,1),(1,0),(1,1),(1,0)],2)
        # n = 2

        n = num_full_weekends([(0,1),(1,0),(0,1),(0,0)],2)
        # n = 0
        
    """
    if wkendtype == 2:
        L1 = [sum(j) for j in x]
        n = sum([(1 if j == 2 else 0) for j in L1])
    else:
        n = 0
        for j in range(len(x)):
            if j < len(x) - 1:
                if x[j][1] == 1 and x[j+1][0] == 1:
                    n += 1
            else:
                if x[j][1] == 1 and x[0][0] == 1:
                    n += 1

    return n


def num_half_weekends(x,wkendtype):
    """
    Returns number of half weekends (one day) worked in a given weekends worked pattern.

    Inputs:
          x - list of 2-tuples representing weekend days worked. Each list
            element is one week. The tuple of binary values represent the
            first and second day of the weekend for that week. A 1 means
            the day is worked, a 0 means it is off.

          wkendtype - 1 --> weekend consists of Saturday and Sunday
                 2 --> weekend consists of Friday and Saturday

    Output:
        Number of half weekends worked

    Example:
        n = num_half_weekends([(0,1),(1,0),(0,1),(1,0)],1)
        # n = 0

        n = num_half_weekends([(0,1),(1,0),(0,1),(0,0)],1)
        # n = 1

        n = num_half_weekends([(1,1),(1,0),(1,1),(1,0)],2)
        # n = 2

        n = num_half_weekends([(0,1),(1,0),(0,1),(0,0)],2)
        # n = 3
        
    """
    if wkendtype == 2:
        L1 = [sum(j) for j in x]
        n = sum([(1 if j == 1 else 0) for j in L1])
    else:
        n = 0
        for j in range(len(x)):
            if j < len(x) - 1:
                if x[j][1] + x[j+1][0] == 1:
                    n += 1
            else:
                if x[j][1] + x[0][0] == 1:
                    n += 1

    return n

def ntuples(lst, n):
    return zip(*[lst[i:]+lst[:i] for i in range(n)])


def dmd_min_to_dat(gmpl_param_name,fn_dmd_or_min,mode='unsliced',isStringIO=True):
    # The input file containing demand by period is assumed to contain
    # each day on a separate row. The number of columns is the same as the
    # number of periods per day. If the file was for a two week problem with
    # half-hour bins, there would be 14 rows and 48 columns. Demand here is
    # really the target staffing level and can be a real number.
    fin = open(fn_dmd_or_min,"r")

    # Read all the lines, strip the trailing spaces, split on the columns
    # and cast the resultant strings to floats. We end up with a 2D array
    # implemented as a list of lists of floats. Done with the input file.
    days = fin.readlines()
    days = [day.rstrip() for day in days]
    days = [day.split() for day in days]
    for day in days:
        day[:] = [float(dmd) for dmd in day]
    fin.close

    # We always assume a 7 day week.
    num_weeks = len(days)/7
    # Not checking for missing or extra columns. Assuming input file
    # creator got it right.
    num_prds = len(days[0])

    # Now it's time to write the GMPL code for this input data element.


    if mode == 'sliced':
        param = 'param ' + gmpl_param_name + ' := '
        for week in range(1,num_weeks + 1):
            # Write the GMPL indexed parameter slice specifier for the week
            weekheader = '\n[*,*,{0}] :'.format(week) + '\n'
            weekheader += ' '.join(map(str, range(1,8)))
            weekheader += ' :=\n'
            param += weekheader

            # Need to transpose the demand values so that days become cols
            # and bins become rows and then write out the GMPL matrix.
            for prd in range(num_prds):
                prd_line = []
                prd_line.append(prd+1)
                prd_line.extend([days[(week-1)*7+day][prd] for day in range(7)])

                prd_line_out = '{0:3d}{1:7.2f}{2:7.2f}{3:7.2f}{4:7.2f}{5:7.2f}{6:7.2f}{7:7.2f}'.format(*prd_line)
                prd_line_out += '\n'
                param += prd_line_out

    else: # Unsliced format
        weeks_of_dmd = []
        for week in range(1,num_weeks + 1):
            week_of_days = []
            for day in range(1,8):
                week_of_days.append(days[7*(week-1) + day - 1])
            weeks_of_dmd.append(week_of_days)

        # Need to reverse the index list so that it is period, day and
        # matches the parameter spec. 
        param = list_to_param(gmpl_param_name, weeks_of_dmd, reverseidx=True)




    # End the GMPL parameter spec and close the file
    if isStringIO:
        paramout = StringIO.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


##param dmd_staff := [*,*,1] :
##        1     2     3     4     5     6     7 :=
##  1   5.0   4.0   4.0   4.0   5.0   5.0   5.0
##  2   5.0   4.0   4.0   4.0   5.0   5.0   5.0
##  3   5.0   4.0   4.0   4.0   5.0   5.0   5.0
##  4   5.0   4.0   4.0   4.0   5.0   5.0   5.0


def scalar_to_param(pname,pvalue,isStringIO=True):
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
        paramout = StringIO.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param

def list_to_param(pname,plist,reverseidx=False,isStringIO=True):
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
    parray = array(plist)

    # Denumerate the array to get at the index tuple and array value
    paramrows = ndenumerate(parray)
    param = 'param ' + pname + ':=\n'
    for pos, val in paramrows:
        poslist = [str(p + 1) for p in pos]
        if reverseidx:
            poslist.reverse()
        datarow = ' '.join(poslist) + ' ' + str(val) + '\n'
        param += datarow

    param += ";\n"
    if isStringIO:
        paramout = StringIO.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param

def shiftlencons_to_param(pname,ttspec,plist,isStringIO=True):
    """
    Convert the shift length specific inputs for the days worked and periods
    worked constraints to a GMPL representation of a parameter.
    Cannot use the generic list_to_param function above since the potentially
    jagged nature of the lists storing these parameters makes it impossible to
    convert to a numpy array for denumeration.

    Inputs:
        pname - string name of paramter in GMPL file
        plist - list containing parameter (could be N-Dimen list)
        isStringIO - true to return StringIO object, false to return string

    Output:
       param tt_shiftlen_min_dys_weeks:=
       1 6 1 3 
       1 6 2 5 
       1 6 3 5 
       1 6 4 5 
       ...

    """


    lengths = get_lengths_from_mix(ttspec)
    param = 'param ' + pname + ':=\n'


    for t in range(0,len(plist)):  # Outer loop is tour types in mix
        t_x = ttspec['tourtypes'][t]['ttnum'] # Get tour type number
        for s in range(0,len(plist[t])):    # Inner loop is shift length
        # Get shift length index
            s_x = lengths.index(ttspec['tourtypes'][t]['shiftlengths'][s]['numbins'])
        # Generate the GMPL rows for this tour type, shift length
            for w in range(0,len(plist[t][s])):
                rowlist = [str(t_x),str(s_x + 1),str(w+1),str(plist[t][s][w])]
                datarow = ' '.join(rowlist) + ' ' + '\n'
                param += datarow

    param += ";\n"
    if isStringIO:
        paramout = StringIO.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param

def list_to_indexedset(sname,slist,isStringIO=True):
    """
    Convert a list to a GMPL representation of a parameter.
    Inputs:
        sname - string name of set in GMPL file
        slist - list containing set (could be N-Dimen list)
        isStringIO - true to return StringIO object, false to return string

    Output:
        set tt_length_x[1] :=
          5 6;

    """
    # Convert set as list to GMPL string rep'n
    gset = ''
    sindex = 0
    for s in slist:
        gset += 'set ' + sname + '[' + str(sindex + 1) + '] :=\n'
        datarow = ' '.join(map(str, s)) + ';\n'
        gset += datarow
        sindex += 1

    if isStringIO:
        gsetout = StringIO.StringIO()
        gsetout.write(gset)
        return gsetout.getvalue()
    else:
        return gset

def mix_days_prds_params(ttspec,pname,nonshiftlen_pname,shiftlen_pname,isStringIO=True):
    """
    Convert the various tour type mix lower and upper bounds (both cumulative
    and non-cumulative and both shift length specific and non-shift length
    specific) to their GMPL parameter representation.

    It's a wrapper function in that it calls list_to_param() for non-shift
    length specific inputs and shiftlencons_to_param() for shift length
    specific inputs.

    Inputs:
        ttspec - the tour type spec object created from the mix file
        pname - string name of paramter in GMPL file
        nonshiftlen_pname - string name of non-shift length specific mix parameter key in YAML file
        shiftlen_pname - string name of shift length specific mix parameter key in YAML file

    Output:
       param tt_shiftlen_min_dys_weeks:=
       1 6 1 3 
       1 6 2 5 
       1 6 3 5 
       1 6 4 5 
       ...

    """

    L = []
    isShiftLen = False
    for m in ttspec['tourtypes']:
        if 'shiftlen' in pname:
            isShiftLen = True
            shiftL = []
            for s in m['shiftlengths']:
                shiftL.append(s[shiftlen_pname])
            L.append(shiftL)
        else:
            if nonshiftlen_pname in m:
                L.append(m[nonshiftlen_pname])
            else:
                L.append(m['shiftlengths'][0][shiftlen_pname])

    if not isShiftLen:
        return list_to_param(pname,L)
    else:
        return shiftlencons_to_param(pname,ttspec,L)



def mix_to_dat(probspec,isStringIO=True):
    """
    Reads a YAML mix file and generates all of the GMPL dat components associated with
    the mix inputs.

    Inputs:
        ttspec - the tour type spec object created from the mix file
        pname - string name of paramter in GMPL file
        nonshiftlen_pname - string name of non-shift length specific mix parameter key in YAML file
        shiftlen_pname - string name of shift length specific mix parameter key in YAML file

    Output:
       param tt_shiftlen_min_dys_weeks:=
       1 6 1 3 
       1 6 2 5 
       1 6 3 5 
       1 6 4 5 
       ...

    """

    # Open the mix file and load it into a YAML object

    fn_mix = probspec['reqd_files']['filename_mix']
    fin = open(fn_mix,"r")
    ttspec = yaml.load(fin)

    mixout = StringIO.StringIO()

##    print ttspec
##    print ttspec['tourtypes']
##    print ttspec['tourtypes'][0]
##    print ttspec['tourtypes'][0]['min_days_week']

    # Get set of shift lengths and order them ascending by length
    lenset = set([])
    for m in ttspec['tourtypes']:
        for s in m['shiftlengths']:
            lenset.add(s['numbins'])
    lengths = list(lenset)
    lengths.sort()
    len_param = list_to_param('lengths', lengths)

    # Number of shift lengths
    n_lengths = size(lengths)
    numlen_param = scalar_to_param('n_lengths', n_lengths)

    # Number of tour types
    n_ttypes = size(ttspec['tourtypes'])
    numttypes_param = scalar_to_param('n_tts', n_ttypes)

    # Tour type length sets
    lenxset = get_length_x_from_mix(ttspec)

    lenxset_set = list_to_indexedset('tt_length_x', lenxset)

    # Midnight threshold for weekend assignments
    midthresholds = [m['midnight_thresh'] for m in ttspec['tourtypes']]
    midthresh_param = list_to_param('midnight_thresh', midthresholds)


    # Parttime flag and bound
    ptflags = [m['is_parttime'] for m in ttspec['tourtypes']]
    ptflags_param = list_to_param('tt_parttime', ptflags)

    ptfrac = ttspec['max_parttime_frac']
    ptfrac_param = scalar_to_param('max_parttime_frac', ptfrac)

    # Global start window width
    width = ttspec['g_start_window_width']
    width_param = scalar_to_param('g_start_window_width', width)

    # Lower and upper bounds on number scheduled
    if 'opt_files' in probspec and 'filename_ttbounds' in probspec['opt_files']:
        fn_ttbnds = probspec['opt_files']['filename_ttbounds']
        fin_ttbnds = open(fn_ttbnds,"r")
        ttbndsspec = yaml.load(fin_ttbnds)
        tt_lb = [m['tt_lb'] for m in ttbndsspec['tourtypes']]
        tt_lb_param = list_to_param('tt_lb', tt_lb)
        tt_ub = [m['tt_ub'] for m in ttbndsspec['tourtypes']]
        tt_ub_param = list_to_param('tt_ub', tt_ub)
    else:
        tt_lb = [m['tt_lb'] for m in ttspec['tourtypes']]
        tt_lb_param = list_to_param('tt_lb', tt_lb)
        tt_ub = [m['tt_ub'] for m in ttspec['tourtypes']]
        tt_ub_param = list_to_param('tt_ub', tt_ub)

    # Cost multiplier
    tt_cost_multiplier = [m['tt_cost_multiplier'] for m in ttspec['tourtypes']]
    tt_cost_multiplier_param = list_to_param('tt_cost_multiplier',
                                                    tt_cost_multiplier)


    # Min and max cumulative days and prds worked over the weeks

    tt_min_dys_weeks_param = mix_days_prds_params(ttspec,
                   'tt_min_dys_weeks','min_days_week',
                                   'min_shiftlen_days_week')

    tt_max_dys_weeks_param = mix_days_prds_params(ttspec,
                   'tt_max_dys_weeks','max_days_week',
                                   'max_shiftlen_days_week')

    tt_min_prds_weeks_param = mix_days_prds_params(ttspec,
                   'tt_min_prds_weeks','min_prds_week',
                                   'min_shiftlen_prds_week')

    tt_max_prds_weeks_param = mix_days_prds_params(ttspec,
                   'tt_max_prds_weeks','max_prds_week',
                                   'max_shiftlen_prds_week')



    # Min and max days and prds worked over the weeks
    # for each shift length workable in the tour type

    tt_shiftlen_min_dys_weeks_param = mix_days_prds_params(ttspec,
                   'tt_shiftlen_min_dys_weeks','min_days_week',
                                   'min_shiftlen_days_week')

    tt_shiftlen_max_dys_weeks_param = mix_days_prds_params(ttspec,
                   'tt_shiftlen_max_dys_weeks','max_days_week',
                                   'max_shiftlen_days_week')

    tt_shiftlen_min_prds_weeks_param = mix_days_prds_params(ttspec,
                   'tt_shiftlen_min_prds_weeks','min_prds_week',
                                   'min_shiftlen_prds_week')

    tt_shiftlen_max_prds_weeks_param = mix_days_prds_params(ttspec,
                   'tt_shiftlen_max_prds_weeks','max_prds_week',
                                   'max_shiftlen_prds_week')




    # Min and max days and prds worked each week

    tt_min_cumul_dys_weeks_param = mix_days_prds_params(ttspec,
                   'tt_min_cumul_dys_weeks','min_cumul_days_week',
                                   'min_shiftlen_cumul_days_week')

    tt_max_cumul_dys_weeks_param = mix_days_prds_params(ttspec,
                   'tt_max_cumul_dys_weeks','max_cumul_days_week',
                                   'max_shiftlen_cumul_days_week')

    tt_min_cumul_prds_weeks_param = mix_days_prds_params(ttspec,
                   'tt_min_cumul_prds_weeks','min_cumul_prds_week',
                                   'min_shiftlen_cumul_prds_week')

    tt_max_cumul_prds_weeks_param = mix_days_prds_params(ttspec,
                   'tt_max_cumul_prds_weeks','max_cumul_prds_week',
                                   'max_shiftlen_cumul_prds_week')


    # Min and max cumulative days and prds worked over the weeks
    # for each shift length workable in the tour type

    tt_shiftlen_min_cumul_dys_weeks_param = mix_days_prds_params(ttspec,
                   'tt_shiftlen_min_cumul_dys_weeks','min_cumul_days_week',
                                   'min_shiftlen_cumul_days_week')

    tt_shiftlen_max_cumul_dys_weeks_param = mix_days_prds_params(ttspec,
                   'tt_shiftlen_max_cumul_dys_weeks','max_cumul_days_week',
                                   'max_shiftlen_cumul_days_week')

    tt_shiftlen_min_cumul_prds_weeks_param = mix_days_prds_params(ttspec,
                   'tt_shiftlen_min_cumul_prds_weeks','min_cumul_prds_week',
                                   'min_shiftlen_cumul_prds_week')

    tt_shiftlen_max_cumul_prds_weeks_param = mix_days_prds_params(ttspec,
                   'tt_shiftlen_max_cumul_prds_weeks','max_cumul_prds_week',
                                   'max_shiftlen_cumul_prds_week')


    # Put the parameter pieces together into a single StringIO object
    print >>mixout, numlen_param
    print >>mixout, len_param
    print >>mixout, numttypes_param
    print >>mixout, lenxset_set
    print >>mixout, midthresh_param
    print >>mixout, tt_lb_param
    print >>mixout, tt_ub_param
    print >>mixout, tt_cost_multiplier_param
    print >>mixout, ptflags_param
    print >>mixout, ptfrac_param
    print >>mixout, width_param

    print >>mixout, tt_min_cumul_dys_weeks_param
    print >>mixout, tt_max_cumul_dys_weeks_param
    print >>mixout, tt_min_cumul_prds_weeks_param
    print >>mixout, tt_max_cumul_prds_weeks_param

    print >>mixout, tt_min_dys_weeks_param
    print >>mixout, tt_max_dys_weeks_param
    print >>mixout, tt_min_prds_weeks_param
    print >>mixout, tt_max_prds_weeks_param

    print >>mixout, tt_shiftlen_min_dys_weeks_param
    print >>mixout, tt_shiftlen_max_dys_weeks_param
    print >>mixout, tt_shiftlen_min_prds_weeks_param
    print >>mixout, tt_shiftlen_max_prds_weeks_param

    print >>mixout, tt_shiftlen_min_cumul_dys_weeks_param
    print >>mixout, tt_shiftlen_max_cumul_dys_weeks_param
    print >>mixout, tt_shiftlen_min_cumul_prds_weeks_param
    print >>mixout, tt_shiftlen_max_cumul_prds_weeks_param

    # print mixout.getvalue()

    if isStringIO:
        return mixout.getvalue()
    else:
        smixout = mixout.read()
        return smixout

def get_length_x_from_mix(ttspec):
    """
    Get list of lists of shift length indexes for each tour type from
    a mix spec.
    Inputs:
        ttspec - yaml representation of tour type mix parameters
    Output:
        A list of lists whose elements are the shift length indexes for
        each tour type.

        Example: [[1,2],[2]]
    """
    # Get set of shift lengths and order them ascending by length
    lenset = set([])
    for m in ttspec['tourtypes']:
        for s in m['shiftlengths']:
            lenset.add(s['numbins'])
    lengths = list(lenset)
    lengths.sort()
    lenxset = []
    for m in ttspec['tourtypes']:
        shifts = [lengths.index(s['numbins']) for s in m['shiftlengths']]
        shifts = [s + 1 for s in shifts]
        shifts.sort()
        lenxset.append(shifts)

    return lenxset

def get_lengths_from_mix(ttspec):
    """
    Get set of shift lengths and order them ascending by length
    Inputs:
        ttspec - yaml representation of tour type mix parameters
    Output:
        A sorted list of shift lengths.

        Example: [8, 16, 20, 24]
    """
    #
    lenset = set([])
    for m in ttspec['tourtypes']:
        for s in m['shiftlengths']:
            lenset.add(s['numbins'])
    lengths = list(lenset)
    lengths.sort()


    return lengths

def csvrow_to_yaml(fn_csv, isStringIO=True):
    """
    Convert a comma delimited row of data into a
    a yaml representation that can be inserted into the yaml mix file.

    This procedure does not not know or care what each row means in the sense
    It's just taking a comma or semicolon delimited row and converts it to yaml.

    Inputs:
        fn_csv - csv filename containing rows of size n_periods_per_day
        isStringIO - true to return StringIO object, false to return string
    Output:
        yaml version of csv row of data either as a StringIO
        object or a string.

        Example:
            Input:     0, 1, 0, 0
            Output:   [0, 1, 0, 0]

    """
    
    fin = open(fn_csv,'r')
    dialect = csv.Sniffer().sniff(fin.read(1024),delimiters=',;')
    fin.seek(0)
    ash_data = csv.reader(fin,dialect)
    
    ash_list = [map(float,row) for row in ash_data]
    fin.close

    yamlstr = ''
    for row in ash_list:
        yamlstr += (' - ' + str(row) + '\n')

    if isStringIO:
        yamlout = StringIO.StringIO()
        yamlout.write(yamlstr)
        return yamlout.getvalue()
    else:
        return yamlstr


def ash_to_dat(fn_yni,fn_mix,isStringIO=True):
    """
    Convert allowable shift start time inputs into GMPL dat form.
    Inputs:
        fn_yni - filename of yaml ini scenario file
        fn_mix - filename of yaml tour type mix file
        isStringIO - true to return StringIO object, false to return string
    Output:
        GMPL dat code for allowable shift start times either as a StringIO
        object or a string.

        Example:
            param allow_start:=
              1 1 1 2 0.0
              2 1 1 2 0.0
              3 1 1 2 0.0
              4 1 1 2 0.0
              ...
              13 1 1 2 1.0
              14 1 1 2 1.0
              15 1 1 2 1.0
    """
    fin_yni = open(fn_yni,"r")
    probspec = yaml.load(fin_yni)
    fin_mix = open(fn_mix,"r")
    ttspec = yaml.load(fin_mix)

# param allow_start[i,j,t,s] = 1 if period i and day j is an allowable
#    shift start time for shift length s of tour type t

    lenxset = get_lengths_from_mix(ttspec)
    ash_rows = []
    for m in ttspec['tourtypes']:
        for s in m['shiftlengths']:
            for j in range(len(s['allowable_starttimes'])):
                for i in range(len(s['allowable_starttimes'][j])):
                    length_x = lenxset.index(s['numbins'])
                    L = [i+1,j+1,length_x+1,m['ttnum'],s['allowable_starttimes'][j][i]]
                    ash_rows.append(L)

    param = 'param allow_start:=\n'
    for val in ash_rows:
        datarow = ' '.join(map(str, val)) + '\n'
        param += datarow

    param += ";\n"
    if isStringIO:
        paramout = StringIO.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param

def wkends_to_dat(fn_yni,fn_mix,isStringIO=True):
    fin_yni = open(fn_yni,"r")
    probspec = yaml.load(fin_yni)
    fin_mix = open(fn_mix,"r")
    ttspec = yaml.load(fin_mix)

    n_weeks = probspec['time']['n_weeks']
    n_ttypes = size(ttspec['tourtypes'])
    patterns_all = create_weekend_base(n_weeks)
    wkend_patterns = []

    wkend_days = [[],[]]
    wkend_days[0] = [1,7]
    wkend_days[1] = [6,7]

    wkend_rows = []
    num_wkend_rows = []

    for m in ttspec['tourtypes']:
        tt = m['ttnum']

        wkend_patterns = [[],[]]
        wkend_patterns[0] = [row for row in patterns_all if filterpatterns(row,tt,1,ttspec)]
        wkend_patterns[1] = [row for row in patterns_all if filterpatterns(row,tt,2,ttspec)]


# param A[p,j,w,t,e] = 1 if weekend pattern p calls for work on day j of week k for tour type t having weekend type e and 0 otherwise


    for i in range(2):
        for t in range(1,n_ttypes+1):
            for p in range(len(wkend_patterns[i])):
                for w in range(n_weeks):
                    for j in range(2):
                        L = [p+1,wkend_days[i][j],w+1,t,i+1,wkend_patterns[i][p][w][j]]                       
                        wkend_rows.append(L)


    for t in range(1,n_ttypes+1):
        for i in range(2):
            L = [i+1,t,len(wkend_patterns[i])]
            num_wkend_rows.append(L)


    param = 'param num_weekend_patterns:=\n'
    for val in num_wkend_rows:
        datarow = ' '.join(map(str, val)) + '\n'
        param += datarow
    param += ";\n"

    param += '\nparam A:=\n'
    for val in wkend_rows:
        datarow = ' '.join(map(str, val)) + '\n'
        param += datarow

    param += ";\n"
    if isStringIO:
        paramout = StringIO.StringIO()
        paramout.write(param)
        return paramout.getvalue()
    else:
        return param


def tester():
    #print csvrow_to_yaml('infiles/oneweekash.csv',False)
    p = create_weekend_base(4)

##    p = [(0,1),(1,1),(0,0),(1,0)]
##    n = num_full_weekends(p,1)

def mwts_createdat(fn_yni,fn_dat):
    """
    Create a GMPL dat file for mwts problems.
    Inputs:
        fn_yni - Name of YAML input file for the mwts problem
    Output:
        fn_dat - Name of GMPL dat file to create
    """
    fin = open(fn_yni,"r")
    probspec = yaml.load(fin)

    # General section
    num_prds_per_day_param = scalar_to_param('n_prds_per_day',
                                      probspec['time']['n_prds_per_day'])

    num_days_per_week_param = scalar_to_param('n_days_per_week',
                                      probspec['time']['n_days_per_week'])

    num_weeks_param = scalar_to_param('n_weeks',
                                      probspec['time']['n_weeks'])

    # Cost related

    labor_budget_param = scalar_to_param('labor_budget',probspec['cost']
                                       ['labor_budget'])

    cu1_param = scalar_to_param('cu1',probspec['cost']
                                       ['understaff_cost_1'])

    cu2_param = scalar_to_param('cu2',probspec['cost']
                                       ['understaff_cost_2'])

    usb_param = scalar_to_param('usb',probspec['cost']
                                       ['understaff_1_lb'])

    # Demand section
    dmd_dat = dmd_min_to_dat('dmd_staff',probspec['reqd_files']['filename_dmd'],mode='unsliced')
    # Min staff section
    min_dat = dmd_min_to_dat('min_staff',probspec['reqd_files']['filename_min'],mode='unsliced')
    # Mix section
    mix_dat = mix_to_dat(probspec)

    # Weekends worked patterns section
    wkends_dat = wkends_to_dat(fn_yni,probspec['reqd_files']['filename_mix'])

    # Allowable shift start time section
    ash_dat = ash_to_dat(fn_yni,probspec['reqd_files']['filename_mix'])



    # Put the pieces together
    dat = StringIO.StringIO()

    print >>dat, num_prds_per_day_param
    print >>dat, num_days_per_week_param
    print >>dat, num_weeks_param

    print >>dat, labor_budget_param
    print >>dat, cu1_param
    print >>dat, cu2_param
    print >>dat, usb_param

    print >>dat, mix_dat
    print >>dat, dmd_dat
    print >>dat, min_dat
    print >>dat, wkends_dat
    print >>dat, ash_dat

    fout = open(fn_dat,"w")
    print >>fout, dat.getvalue()
    fout.close()


def main():
##    dmd_min_to_dat('dmd_staff','infiles/jax_4week.dmd','test.out')
##    testout = dmd_min_to_datstringio('dmd_staff','infiles/jax_4week.dmd')
##    fout = open('teststringio.out',"w")
##    print "Starting to write stringio"
##    print >>fout, testout
##    print "Ending write stringio"
##    fout.close()
    mwts_createdat('/home/mark/Documents/research/MultiWeek/exps/mwts01/inputs/simple.yni','/home/mark/Documents/research/MultiWeek/wsmwts/pymwts/tests/simple.dat')
    #tester()

if __name__ == '__main__':
    main()
