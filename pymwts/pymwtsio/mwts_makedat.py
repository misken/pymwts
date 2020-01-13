"""
Read input files for mwts problems and create a GMPL data file.
"""

# Author: misken
# License: TBD

import io
import csv
import itertools

import yaml
import numpy as np
import pandas as pd


def dmd_min_to_dat(gmpl_param_name, fn_dmd_or_min, sep=',', header=None, comment='#',
                   mode='unsliced', isstringio=True):
    """
    Convert plain text matrices of dmd or min data into their dat equivalents.

    The input file containing demand by period is assumed to contain
    each day on a separate row. The number of columns is the same as the
    number of periods per day. If the file was for a two week problem with
    half-hour epoch_tuples, there would be 14 rows and 48 columns. Demand here is
    really the target staffing level and can be a real number.

    Default is comma separated file with no header line.

    :param gmpl_param_name: name for resulting GMPL param value
    :type gmpl_param_name: str
    :param fn_dmd_or_min: filename for target or minimum staffing levels
    :type fn_dmd_or_min: str
    :param sep: delimiter (default is ',')
    :param header: see Pandas `read_csv()` docs (default is None)
    :param comment: comment character (default is '#')
    :param mode: GMPL dat format; either 'sliced' or 'unsliced' (default)
    :param isstringio: True (default) to return StringIO object, False to return string
    :return:
    """

    # with open(fn_dmd_or_min, "r") as f_in:
    #
    #     # Read all the lines, strip the trailing spaces, split on the columns
    #     # and cast the resultant strings to floats. We end up with a 2D array
    #     # implemented as a list of lists of floats.
    #     lines_raw = f_in.readlines()
    #     lines = [line.rstrip() for line in lines_raw]
    #     days = [day.split() for line in lines]
    #     for day in days:
    #         day[:] = [float(dmd) for dmd in day]

    days_csv = pd.read_csv(fn_dmd_or_min, sep=sep, header=header, comment=comment)
    # Convert to numpy array which we can treat like list of lists
    days = days_csv.values

    # We always assume a 7 day week.
    num_weeks = len(days) // 7
    # Not checking for missing or extra columns. Assuming input file
    # creator got it right.
    num_prds = len(days[0])

    # Check if writing sliced or unsliced format
    if mode == 'sliced':
        param = 'param ' + gmpl_param_name + ' := '
        for week in range(1, num_weeks + 1):
            # Write the GMPL indexed parameter slice specifier for the week
            weekheader = '\n[*,*,{0}] :'.format(week) + '\n'
            weekheader += ' '.join(map(str, range(1, 8)))
            weekheader += ' :=\n'
            param += weekheader

            # Need to transpose the demand values so that days become cols
            # and epoch_tuples become rows and then write out the GMPL matrix.
            for prd in range(num_prds):
                prd_line = [prd + 1]
                prd_line.extend([days[(week - 1) * 7 + day][prd] for day in range(7)])

                prd_line_out = '{0:3d}{1:7.2f}{2:7.2f}{3:7.2f}{4:7.2f}{5:7.2f}{6:7.2f}{7:7.2f}'.format(*prd_line)
                prd_line_out += '\n'
                param += prd_line_out

    else:  # Unsliced format
        weeks_of_dmd = []
        for week in range(1, num_weeks + 1):
            week_of_days = []
            for day in range(1, 8):
                week_of_days.append(days[7 * (week - 1) + day - 1])
            weeks_of_dmd.append(week_of_days)

        # Need to reverse the index list so that it is period, day and
        # matches the parameter spec. 
        param = list_to_param(gmpl_param_name, weeks_of_dmd, reverseidx=True)

    if isstringio:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def scalar_to_param(gmpl_param_name, param_value, isstringio=True):
    """
    Convert a scalar to a GMPL representation of a parameter.

    :param gmpl_param_name: name for resulting GMPL param value
    :type gmpl_param_name: str
    :param param_value:
    :param isstringio: True (default) to return StringIO object, False to return string
    :return: GMPL dat code for scalar parameter either as a StringIO object or a string.
    :rtype: StringIO or str

    Example:

     scalar_to_param('n_prds_per_day', 48) -->  'param n_prds_per_day :=  48;'
    """

    param = 'param ' + gmpl_param_name + ' :=  ' + str(param_value) + ';\n'
    if isstringio:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def list_to_param(gmpl_param_name, param_list, reverseidx=False, isstringio=True):
    """
    Convert a list to a GMPL representation of a parameter.

    :param gmpl_param_name: name for resulting GMPL param value
    :param param_list: list containing parameter (could be N-Dimen list)
    :param reverseidx: True to reverse the order of the indexes (essentially transposing the matrix)
    :param isstringio: True to return StringIO object, False to return string
    :return: GMPL dat code for list parameter
    :rtype:  StringIO object or a string.

    Example:

     scalar_to_param('midnight_thresh', [100, 100, 100]) -->

        param midnight_thresh:=
             1 100
             2 100
             3 100
            ;
    """

    # Convert parameter as list to an ndarray
    param_array = np.array(param_list)

    # Denumerate the array to get at the index tuple and array value
    param_rows = np.ndenumerate(param_array)
    param = 'param ' + gmpl_param_name + ':=\n'
    for pos, val in param_rows:
        pos_list = [str(p + 1) for p in pos]
        if reverseidx:
            pos_list.reverse()
        data_row = ' '.join(pos_list) + ' ' + str(val) + '\n'
        param += data_row

    param += ";\n"
    if isstringio:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def shiftlencons_to_param(gmpl_param_name, ttspec, param_list, isstringio=True):
    """
    Convert the shift length specific inputs for the days worked and periods
    worked constraints to a GMPL parameter representation.

    Cannot use the generic list_to_param function above since the potentially
    jagged nature of the lists storing these parameters makes it impossible to
    convert to a numpy array for denumeration.

    :param gmpl_param_name:
    :param ttspec:
    :param param_list:
    :param isstringio: True to return StringIO object, False to return string
    :return: GMPL dat code for list parameter
    :rtype:  StringIO object or a string.


    Example:

       param tt_shiftlen_min_dys_weeks:=
       1 6 1 3 
       1 6 2 5 
       1 6 3 5 
       1 6 4 5 
       ...

    """

    lengths = get_lengths_from_mix(ttspec)
    param = 'param ' + gmpl_param_name + ':=\n'

    for t in range(0, len(param_list)):  # Outer loop is tour types in mix
        t_x = ttspec['tourtypes'][t]['ttnum']  # Get tour type number
        for s in range(0, len(param_list[t])):  # Inner loop is shift length
            # Get shift length index
            s_x = lengths.index(ttspec['tourtypes'][t]['shiftlengths'][s]['numbins'])
            # Generate the GMPL rows for this tour type, shift length
            for w in range(0, len(param_list[t][s])):
                row_list = [str(t_x), str(s_x + 1), str(w + 1), str(param_list[t][s][w])]
                data_row = ' '.join(row_list) + ' ' + '\n'
                param += data_row

    param += ";\n"
    if isstringio:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def list_to_indexedset(gmpl_set_name, set_list, isstringio=True):
    """
    Convert a list to a GMPL representation of a parameter.

    :param gmpl_set_name: name for resulting GMPL param value
    :param set_list: list containing set (could be N-Dimen list)
    :param isstringio: True to return StringIO object, False to return string
    :return: GMPL dat code for list parameter
    :rtype:  StringIO object or a string.

    Example:

        set tt_length_x[1] :=
          5 6;

    """

    gmpl_set = ''
    set_index = 0
    for s in set_list:
        gmpl_set += 'set ' + gmpl_set_name + '[' + str(set_index + 1) + '] :=\n'
        data_row = ' '.join(map(str, s)) + ';\n'
        gmpl_set += data_row
        set_index += 1

    if isstringio:
        gmpl_set_out = io.StringIO()
        gmpl_set_out.write(gmpl_set)
        return gmpl_set_out.getvalue()
    else:
        return gmpl_set


def mix_days_prds_params(ttspec, param_name, non_shiftlen_param_name, shiftlen_param_name):
    """
    Convert the various tour type mix lower and upper bounds (both cumulative
    and non-cumulative and both shift length specific and non-shift length
    specific) to their GMPL parameter representation.

    It's a wrapper function in that it calls list_to_param() for non-shift
    length specific inputs and shiftlencons_to_param() for shift length
    specific inputs.

    :param ttspec: the tour type spec object created from the mix file
    :param param_name: name for resulting GMPL param value
    :param non_shiftlen_param_name: string name of non-shift length specific mix parameter key in YAML file
    :param shiftlen_param_name: string name of shift length specific mix parameter key in YAML file
    :return: GMPL dat code for parameters
    :rtype:  StringIO object or a string.

    Example:

       param tt_shiftlen_min_dys_weeks:=
       1 6 1 3 
       1 6 2 5 
       1 6 3 5 
       1 6 4 5 
       ...

    """

    param_list = []
    is_shift_len = False
    for m in ttspec['tourtypes']:
        if 'shiftlen' in param_name:
            is_shift_len = True
            shift_lens = []
            for s in m['shiftlengths']:
                shift_lens.append(s[shiftlen_param_name])
            param_list.append(shift_lens)
        else:
            if non_shiftlen_param_name in m:
                param_list.append(m[non_shiftlen_param_name])
            else:
                param_list.append(m['shiftlengths'][0][shiftlen_param_name])

    if not is_shift_len:
        return list_to_param(param_name, param_list)
    else:
        return shiftlencons_to_param(param_name, ttspec, param_list)


def mix_to_dat(prob_spec, isstringio=True):
    """
    Read a YAML mix file and generates all of the GMPL dat components associated with the mix inputs.

    :param prob_spec: Problem spec object resulting from loading YAML project yni file.
    :param isstringio: True to return StringIO object, False to return string
    :return: GMPL dat code for mix related parameters
    :rtype:  StringIO object or a string.
    """

    # Open the mix file and load it into a YAML object

    fn_mix = prob_spec['reqd_files']['filename_mix']
    f_mix_in = open(fn_mix, "r")
    ttspec = yaml.safe_load(f_mix_in)
    f_mix_in.close()

    mix_out = io.StringIO()

    # Get set of shift lengths and order them ascending by length
    len_set = set([])
    for m in ttspec['tourtypes']:
        for s in m['shiftlengths']:
            len_set.add(s['numbins'])
    lengths = list(len_set)
    lengths.sort()
    len_param = list_to_param('lengths', lengths)

    # Number of shift lengths
    n_lengths = np.size(lengths)
    n_len_param = scalar_to_param('n_lengths', n_lengths)

    # Number of tour types
    n_ttypes = np.size(ttspec['tourtypes'])
    n_ttypes_param = scalar_to_param('n_tts', n_ttypes)

    # Tour type length sets
    lenx_set = get_length_x_from_mix(ttspec)
    lenx_set_set = list_to_indexedset('tt_length_x', lenx_set)

    # Midnight threshold for weekend assignments
    # midthresholds = [m['midnight_thresh'] for m in ttspec['tourtypes']]
    # midthresh_param = list_to_param('midnight_thresh', midthresholds)

    # Parttime flag and bound
    pt_flags = [m['is_parttime'] for m in ttspec['tourtypes']]
    pt_flags_param = list_to_param('tt_parttime', pt_flags)

    pt_frac = ttspec['max_parttime_frac']
    pt_frac_param = scalar_to_param('max_parttime_frac', pt_frac)

    # Global start window width
    width = ttspec['g_start_window_width']
    width_param = scalar_to_param('g_start_window_width', width)

    # Lower and upper bounds on number scheduled
    if 'opt_files' in prob_spec and 'filename_ttbounds' in prob_spec['opt_files']:
        fn_ttbnds = prob_spec['opt_files']['filename_ttbounds']
        fin_ttbnds = open(fn_ttbnds, "r")
        ttbnds_spec = yaml.safe_load(fin_ttbnds)
        tt_lb = [m['tt_lb'] for m in ttbnds_spec['tourtypes']]
        tt_lb_param = list_to_param('tt_lb', tt_lb)
        tt_ub = [m['tt_ub'] for m in ttbnds_spec['tourtypes']]
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
                                                  'tt_min_dys_weeks', 'min_days_week',
                                                  'min_shiftlen_days_week')

    tt_max_dys_weeks_param = mix_days_prds_params(ttspec,
                                                  'tt_max_dys_weeks', 'max_days_week',
                                                  'max_shiftlen_days_week')

    tt_min_prds_weeks_param = mix_days_prds_params(ttspec,
                                                   'tt_min_prds_weeks', 'min_prds_week',
                                                   'min_shiftlen_prds_week')

    tt_max_prds_weeks_param = mix_days_prds_params(ttspec,
                                                   'tt_max_prds_weeks', 'max_prds_week',
                                                   'max_shiftlen_prds_week')

    # Min and max days and prds worked over the weeks
    # for each shift length workable in the tour type
    tt_shiftlen_min_dys_weeks_param = mix_days_prds_params(ttspec,
                                                           'tt_shiftlen_min_dys_weeks', 'min_days_week',
                                                           'min_shiftlen_days_week')

    tt_shiftlen_max_dys_weeks_param = mix_days_prds_params(ttspec,
                                                           'tt_shiftlen_max_dys_weeks', 'max_days_week',
                                                           'max_shiftlen_days_week')

    tt_shiftlen_min_prds_weeks_param = mix_days_prds_params(ttspec,
                                                            'tt_shiftlen_min_prds_weeks', 'min_prds_week',
                                                            'min_shiftlen_prds_week')

    tt_shiftlen_max_prds_weeks_param = mix_days_prds_params(ttspec,
                                                            'tt_shiftlen_max_prds_weeks', 'max_prds_week',
                                                            'max_shiftlen_prds_week')

    # Min and max days and prds worked each week
    tt_min_cumul_dys_weeks_param = mix_days_prds_params(ttspec,
                                                        'tt_min_cumul_dys_weeks', 'min_cumul_days_week',
                                                        'min_shiftlen_cumul_days_week')

    tt_max_cumul_dys_weeks_param = mix_days_prds_params(ttspec,
                                                        'tt_max_cumul_dys_weeks', 'max_cumul_days_week',
                                                        'max_shiftlen_cumul_days_week')

    tt_min_cumul_prds_weeks_param = mix_days_prds_params(ttspec,
                                                         'tt_min_cumul_prds_weeks', 'min_cumul_prds_week',
                                                         'min_shiftlen_cumul_prds_week')

    tt_max_cumul_prds_weeks_param = mix_days_prds_params(ttspec,
                                                         'tt_max_cumul_prds_weeks', 'max_cumul_prds_week',
                                                         'max_shiftlen_cumul_prds_week')

    # Min and max cumulative days and prds worked over the weeks
    # for each shift length workable in the tour type
    tt_shiftlen_min_cumul_dys_weeks_param = mix_days_prds_params(ttspec,
                                                                 'tt_shiftlen_min_cumul_dys_weeks',
                                                                 'min_cumul_days_week',
                                                                 'min_shiftlen_cumul_days_week')

    tt_shiftlen_max_cumul_dys_weeks_param = mix_days_prds_params(ttspec,
                                                                 'tt_shiftlen_max_cumul_dys_weeks',
                                                                 'max_cumul_days_week',
                                                                 'max_shiftlen_cumul_days_week')

    tt_shiftlen_min_cumul_prds_weeks_param = mix_days_prds_params(ttspec,
                                                                  'tt_shiftlen_min_cumul_prds_weeks',
                                                                  'min_cumul_prds_week',
                                                                  'min_shiftlen_cumul_prds_week')

    tt_shiftlen_max_cumul_prds_weeks_param = mix_days_prds_params(ttspec,
                                                                  'tt_shiftlen_max_cumul_prds_weeks',
                                                                  'max_cumul_prds_week',
                                                                  'max_shiftlen_cumul_prds_week')

    # Put the parameter pieces together into a single StringIO object
    print(mix_out, n_len_param)
    print(mix_out, len_param)
    print(mix_out, n_ttypes_param)
    print(mix_out, lenx_set_set)
    # print >>mixout, midthresh_param
    print(mix_out, tt_lb_param)
    print(mix_out, tt_ub_param)
    print(mix_out, tt_cost_multiplier_param)
    print(mix_out, pt_flags_param)
    print(mix_out, pt_frac_param)
    print(mix_out, width_param)

    print(mix_out, tt_min_cumul_dys_weeks_param)
    print(mix_out, tt_max_cumul_dys_weeks_param)
    print(mix_out, tt_min_cumul_prds_weeks_param)
    print(mix_out, tt_max_cumul_prds_weeks_param)

    print(mix_out, tt_min_dys_weeks_param)
    print(mix_out, tt_max_dys_weeks_param)
    print(mix_out, tt_min_prds_weeks_param)
    print(mix_out, tt_max_prds_weeks_param)

    print(mix_out, tt_shiftlen_min_dys_weeks_param)
    print(mix_out, tt_shiftlen_max_dys_weeks_param)
    print(mix_out, tt_shiftlen_min_prds_weeks_param)
    print(mix_out, tt_shiftlen_max_prds_weeks_param)

    print(mix_out, tt_shiftlen_min_cumul_dys_weeks_param)
    print(mix_out, tt_shiftlen_max_cumul_dys_weeks_param)
    print(mix_out, tt_shiftlen_min_cumul_prds_weeks_param)
    print(mix_out, tt_shiftlen_max_cumul_prds_weeks_param)

    if isstringio:
        return mix_out.getvalue()
    else:
        return mix_out.read()


def get_length_x_from_mix(ttspec):
    """
    Get list of lists of shift length indexes for each tour type from a mix spec.
    
    :param ttspec: YAML representation of tour type mix parameters
    :return: A list of lists whose elements are the shift length indexes for each tour type.

    Example: [[1,2],[2]]
    """

    # Get set of shift lengths and order them ascending by length
    len_set = set([])
    for m in ttspec['tourtypes']:
        for s in m['shiftlengths']:
            len_set.add(s['numbins'])
    lengths = list(len_set)
    lengths.sort()
    lenx_set = []
    for m in ttspec['tourtypes']:
        shifts = [lengths.index(s['numbins']) for s in m['shiftlengths']]
        shifts = [s + 1 for s in shifts]
        shifts.sort()
        lenx_set.append(shifts)

    return lenx_set


def get_lengths_from_mix(ttspec):
    """
    Get set of shift lengths and order them ascending by length

    :param ttspec: YAML representation of tour type mix parameters
    :return: A sorted list of shift lengths.

    Example: [8, 16, 20, 24]
    """

    len_set = set([])
    for m in ttspec['tourtypes']:
        for s in m['shiftlengths']:
            len_set.add(s['numbins'])
    lengths = list(len_set)
    lengths.sort()

    return lengths


def csvrow_to_yaml(fn_csv, isstringio=True):
    """
    Convert a comma delimited row of data into a
    a yaml representation that can be inserted into the yaml mix file.

    This procedure does not not know or care what each row means in the sense
    that it's just taking a comma or semicolon delimited row and converts it to YAML.

    :param fn_csv: csv filename
    :param isstringio: True to return StringIO object, False to return string
    :return: yaml version of csv row of data either as a StringIO object or a string.

    Example:
            Input:     0, 1, 0, 0
            Output:   [0, 1, 0, 0]
    """

    f_csv_in = open(fn_csv, 'r')
    dialect = csv.Sniffer().sniff(f_csv_in.read(1024), delimiters=',;')
    f_csv_in.seek(0)
    ash_data = csv.reader(f_csv_in, dialect)

    ash_list = [map(float, row) for row in ash_data]
    f_csv_in.close()

    yaml_str = ''
    for row in ash_list:
        yaml_str += (' - ' + str(row) + '\n')

    if isstringio:
        yaml_out = io.StringIO()
        yaml_out.write(yaml_str)
        return yaml_out.getvalue()
    else:
        return yaml_str


def ash_to_dat(fn_mix, isstringio=True):
    """
    Convert allowable shift start time inputs into GMPL dat form.

    :param fn_mix: filename of yaml tour type mix file
    :param isstringio: true to return StringIO object, false to return string
    :return: GMPL dat code for allowable shift start times either as a StringIO
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

    f_mix_in = open(fn_mix, "r")
    ttspec = yaml.safe_load(f_mix_in)
    f_mix_in.close()

    # param allow_start[i,j,t,s] = 1 if period i and day j is an allowable
    #    shift start time for shift length s of tour type t

    lenx_set = get_lengths_from_mix(ttspec)
    ash_rows = []
    for m in ttspec['tourtypes']:
        for s in m['shiftlengths']:
            for j in range(len(s['allowable_starttimes'])):
                for i in range(len(s['allowable_starttimes'][j])):
                    length_x = lenx_set.index(s['numbins'])
                    row = [i + 1, j + 1, length_x + 1, m['ttnum'], s['allowable_starttimes'][j][i]]
                    ash_rows.append(row)

    param = 'param allow_start:=\n'
    for val in ash_rows:
        data_row = ' '.join(map(str, val)) + '\n'
        param += data_row

    param += ";\n"
    if isstringio:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def wkends_to_dat(fn_yni, fn_mix, fn_wkd, isstringio=True):
    """

    :param fn_yni:
    :param fn_mix:
    :param fn_wkd:
    :param isstringio:
    :return:
    """

    fin_yni = open(fn_yni, "r")
    probspec = yaml.safe_load(fin_yni)
    fin_mix = open(fn_mix, "r")
    ttspec = yaml.safe_load(fin_mix)

    fin_wkd = open(fn_wkd, "r")
    wkdspec = yaml.safe_load(fin_wkd)

    n_weeks = probspec['time']['n_weeks']
    n_ttypes = np.size(ttspec['tourtypes'])
    patterns_all = create_weekend_base(n_weeks)
    wkend_patterns = []

    wkend_days = [[], []]
    wkend_days[0] = [1, 7]
    wkend_days[1] = [6, 7]

    wkend_rows = []
    num_wkend_rows = []

    for m in ttspec['tourtypes']:
        tt = m['ttnum']

        wkend_patterns = [[], []]
        wkend_patterns[0] = [row for row in patterns_all if filter_patterns(row, tt, 1, wkdspec)]
        wkend_patterns[1] = [row for row in patterns_all if filter_patterns(row, tt, 2, wkdspec)]

    # param A[p,j,w,t,e] = 1 if weekend pattern p calls for work on day j of week k
    # for tour type t having weekend type e and 0 otherwise

    for i in range(2):
        for t in range(1, n_ttypes + 1):
            for p in range(len(wkend_patterns[i])):
                for w in range(n_weeks):
                    for j in range(2):
                        temp_list = [p + 1, wkend_days[i][j], w + 1, t, i + 1, wkend_patterns[i][p][w][j]]
                        wkend_rows.append(temp_list)

    for t in range(1, n_ttypes + 1):
        for i in range(2):
            temp_list = [i + 1, t, len(wkend_patterns[i])]
            num_wkend_rows.append(temp_list)

    param = 'param num_weekend_patterns:=\n'
    for val in num_wkend_rows:
        data_row = ' '.join(map(str, val)) + '\n'
        param += data_row
    param += ";\n"

    param += '\nparam A:=\n'
    for val in wkend_rows:
        data_row = ' '.join(map(str, val)) + '\n'
        param += data_row

    param += ";\n"

    # Midnight threshold for weekend assignments
    # TODO - generalize for tour types, right now assuming global used
    mid_thresholds = [wkdspec['wkd_global']['midnight_thresh'] for _ in ttspec['tourtypes']]
    mid_thresh_param = list_to_param('midnight_thresh', mid_thresholds)
    param += mid_thresh_param

    if isstringio:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def tester():
    p = create_weekend_base(4)
    print(p)

    # #    p = [(0,1),(1,1),(0,0),(1,0)]
    # #    n = num_full_weekends(p,1)


def mwts_createdat(fn_yni, fn_dat):
    """
    Create a GMPL dat file for mwts problems.

    :param fn_yni: Name of YAML input file for the mwts problem
    :param fn_dat: Name of GMPL dat file to create
    :return:
    """

    f_yni_in = open(fn_yni, "r")
    prob_spec = yaml.safe_load(f_yni_in)

    # Scheduling cycle
    num_prds_per_day_param = scalar_to_param('n_prds_per_day',
                                             prob_spec['time']['n_prds_per_day'],
                                             isstringio=False)

    num_days_per_week_param = scalar_to_param('n_days_per_week',
                                              prob_spec['time']['n_days_per_week'],
                                              isstringio=False)

    num_weeks_param = scalar_to_param('n_weeks',
                                      prob_spec['time']['n_weeks'],
                                      isstringio=False)

    # Cost related
    labor_budget_param = scalar_to_param('labor_budget', prob_spec['cost']['labor_budget'],
                                         isstringio=False)

    cu1_param = scalar_to_param('cu1', prob_spec['cost']['understaff_cost_1'],
                                isstringio=False)

    cu2_param = scalar_to_param('cu2', prob_spec['cost']['understaff_cost_2'],
                                isstringio=False)

    usb_param = scalar_to_param('usb', prob_spec['cost']['understaff_1_ub'],
                                isstringio=False)

    # Demand and min staff
    dmd_dat = dmd_min_to_dat('dmd_staff', prob_spec['reqd_files']['filename_dmd'],
                             mode='unsliced', isstringio=False)
    min_dat = dmd_min_to_dat('min_staff', prob_spec['reqd_files']['filename_min'],
                             mode='unsliced', isstringio=False)

    # Tour type mix
    mix_dat = mix_to_dat(prob_spec, isstringio=False)

    # Weekends worked patterns section
    wkends_dat = wkends_to_dat(fn_yni,
                               prob_spec['reqd_files']['filename_mix'],
                               prob_spec['reqd_files']['filename_wkd'],
                               isstringio=False)

    # Allowable shift start time section
    ash_dat = ash_to_dat(prob_spec['reqd_files']['filename_mix'], isstringio=False)

    # Put the pieces together
    dat_list = [num_prds_per_day_param,
                num_days_per_week_param,
                num_weeks_param,
                labor_budget_param,
                cu1_param,
                cu2_param,
                usb_param,
                mix_dat,
                dmd_dat,
                min_dat,
                wkends_dat,
                ash_dat]

    dat_str = '\n'.join(dat_list)

    with open(fn_dat, "w") as f_dat_out:
        f_dat_out.write(dat_str)


def create_yaml_ash():
    """
    Guessing this was a utility used while creating ash files

    :return:
    """
    shift_lens = ['8', '10', '12', '14', '16', '20', '24']

    for shift_len in shift_lens:
        fn_csv = '../exps/mwts02/inputs/ash/ash_' + shift_len + '.fn_csv'
        fn_yml = '../exps/mwts02/inputs/ash/ash_' + shift_len + '.yaml'

        ash = io.StringIO()
        print(ash, csvrow_to_yaml(fn_csv))
        f_out = open(fn_yml, "w")
        print(f_out, ash.getvalue())
        f_out.close()


def create_weekend_base(n_weeks):
    """
    Generate basis for cartesion product of [0,1] lists based
    on number of weeks in scheduling problem.

    Returns a list in which each list
    element is one week. The tuple of binary values represent the
    first and second day of the weekend for that week. A 1 means
    the day is worked, a 0 means it is off. This
    exhaustive list can later be filtered to only include desirable patterns.

    :param n_weeks: number of weeks in scheduling horizon
    :type n_weeks: int
    :return: all of the possible n_weeks weekends worked patterns
    :rtype: list of tuples (of size n_weeks) of lists (of size 2)
    """

    basis_list = [[0, 0], [1, 0], [0, 1], [1, 1]]

    mw_basis_list = []
    for i in range(n_weeks):
        mw_basis_list.append(basis_list)

    return list(itertools.product(*mw_basis_list))


def filter_patterns(pattern, ttnum, wkend_type, wkdspec):
    """
    TODO: Create a sequence of binary values to be used for list filtering.

    This function will contain the various rules used to filter out weekend days
    worked patterns that we don't want to allow.

    For now I'm hard coding in rules but need to develop an approach to
    flexibly specifying rules to apply to filter out undesirable weekends
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

          (2) Type 2, work every other weekend
            pattern = [(1,1),(0,0),(1,1),(0,0)], type = 2

        Output: True --> keep pattern
                False --> discard pattern

    :param pattern:
    :param ttnum:
    :param wkend_type:
    :param wkdspec:
    :return:
    """

    # n_weeks = len(pattern)
    keep = True

    if 'wkd_global' in wkdspec:
        max_days_worked = wkdspec['wkd_global']['max_days_worked']
        is_halfweekend_ok = wkdspec['wkd_global']['half_weekends_ok']
        max_wkends_worked = wkdspec['wkd_global']['max_wkends_worked']
        max_consec_wkends = wkdspec['wkd_global']['max_consec_wkends']

    else:
        # TODO - implement tour type specific global overrides
        # tourtype = [t for t in wkdspec['tourtypes'] if t['ttnum'] == ttnum]

        # PLACEHOLDERS
        max_days_worked = 0
        is_halfweekend_ok = False
        max_wkends_worked = 0
        max_consec_wkends = 0

    # No more than max_days_worked over the scheduling horizon
    # max_days_worked = wkdspec[0]['max_days_worked']
    if not (sum(pattern) <= max_days_worked):
        keep = False

        # Half-weekends
    if not is_halfweekend_ok and num_half_weekends(pattern, wkend_type) > 0:
        keep = False

        # Max weekends (full or half) worked
    tot_wkends_worked = num_half_weekends(pattern, wkend_type) + num_full_weekends(pattern, wkend_type)
    if tot_wkends_worked > max_wkends_worked:
        keep = False

    # Limit on consecutive weekends with one or more days worked
    consec_wkends_worked = num_consecutive_weekends(pattern, wkend_type)
    if consec_wkends_worked > max_consec_wkends:
        keep = False

    return keep


def num_full_weekends(x, wkend_type):
    """
    Compute number of full weekends (both days) worked in a given weekends worked pattern.

    :param x: list of 2-tuples representing weekend days worked. Each list
            element is one week. The tuple of binary values represent the
            first and second day of the weekend for that week. A 1 means
            the day is worked, a 0 means it is off.
    :param wkend_type: 1 --> weekend consists of Saturday and Sunday
                       2 --> weekend consists of Friday and Saturday
    :return: Number of full weekends worked

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

    if wkend_type == 2:
        n_wkend_days_worked = [sum(j) for j in x]
        n = sum([(1 if j == 2 else 0) for j in n_wkend_days_worked])
    else:
        n = 0
        for j in range(len(x)):
            if j < len(x) - 1:
                if x[j][1] == 1 and x[j + 1][0] == 1:
                    n += 1
            else:
                if x[j][1] == 1 and x[0][0] == 1:
                    n += 1

    return n


def num_half_weekends(x, wkend_type):
    """
    Compute number of full weekends (both days) worked in a given weekends worked pattern.

    :param x: list of 2-tuples representing weekend days worked. Each list
            element is one week. The tuple of binary values represent the
            first and second day of the weekend for that week. A 1 means
            the day is worked, a 0 means it is off.
    :param wkend_type: 1 --> weekend consists of Saturday and Sunday
                       2 --> weekend consists of Friday and Saturday
    :return: Number of half weekends worked

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
    if wkend_type == 2:
        n_wkend_days_worked = [sum(j) for j in x]
        n = sum([(1 if j == 1 else 0) for j in n_wkend_days_worked])
    else:
        n = 0
        for j in range(len(x)):
            if j < len(x) - 1:
                if x[j][1] + x[j + 1][0] == 1:
                    n += 1
            else:
                if x[j][1] + x[0][0] == 1:
                    n += 1

    return n


def n_tuples(lst, n):
    """
    No idea why I created this and it is not used.

    :param lst:
    :param n:
    :return:
    """
    return zip(*[lst[i:] + lst[:i] for i in range(n)])


def is_weekend_worked(x, wkend_type, i):
    """
    Determine if i'th weekend is worked (full or half) for a given
    weekend pattern x and weekend type,

    :param x: list of 2-tuples representing weekend days worked. Each list
            element is one week. The tuple of binary values represent the
            first and second day of the weekend for that week. A 1 means
            the day is worked, a 0 means it is off.
    :param wkend_type: 1 --> weekend consists of Saturday and Sunday
                       2 --> weekend consists of Friday and Saturday
    :type wkend_type: int
    :param i: which weekend to check
    :type i: int
    :return: True if weekend i is worked, False otherwise

    Example:
        b = is_weekend_worked([(0,1),(1,0),(0,1),(1,0)],1,1)
        # b = True

        b = is_weekend_worked([(0,1),(1,0),(0,1),(1,0)],1,2)
        # b = False

        b = is_weekend_worked([(1,0),(1,0),(1,1),(1,0)],2,1)
        # b = True

        b = is_weekend_worked([(0,1),(1,0),(0,1),(0,0)],2,4)
        # b = False

    """
    if wkend_type == 2:
        if sum(x[i - 1]):
            return True
        else:
            return False
    else:
        n_weeks = len(x)
        if i < n_weeks:
            sunday_idx = i
        else:
            sunday_idx = 0

        if (x[i - 1][1] + x[sunday_idx][0]) > 0:
            return True
        else:
            return False


def num_consecutive_weekends(x, wkend_type, circular=True):
    """
    Returns largest number of consecutive weekends worked in a pattern.

    :param x: list of 2-tuples representing weekend days worked. Each list
            element is one week. The tuple of binary values represent the
            first and second day of the weekend for that week. A 1 means
            the day is worked, a 0 means it is off.
    :param wkend_type: 1 --> weekend consists of Saturday and Sunday
                       2 --> weekend consists of Friday and Saturday
    :param circular: True  --> week n wraps to week 1
                     False --> weeks do not wrap
    :return: Number of consecutive weekends worked

    Example:
        n = num_consecutive_weekends([(0,1),(1,0),(0,1),(1,0)],1)
        # n = 1

        n = num_consecutive_weekends([(0,1),(1,1),(0,1),(0,0)],1)
        # n = 3

        n = num_consecutive_weekends([(1,1),(1,0),(1,1),(1,0)],2)
        # n = 4

        n = num_consecutive_weekends([(0,1),(1,0),(0,1),(0,0)],2)
        # n = 3

    """

    n = 0
    n_consec = 0
    if circular:
        pattern = x + x
    else:
        pattern = x

    n_weeks = len(pattern)
    for i in range(1, n_weeks + 1):
        if is_weekend_worked(pattern, wkend_type, i):
            n += 1
            if n > n_consec:
                n_consec = n
        else:
            n = 0

    return n_consec


def main():
    pass


if __name__ == '__main__':
    main()
