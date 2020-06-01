"""
Read input files for mwts problems and create a GMPL data file.
"""

# Author: misken
# License: TBD

import io
import re
import json
from datetime import time

import pandas as pd
import numpy as np


def shift_label(start_prd, shift_len, prds_per_day):
    mins_per_prd = 1440.0 / prds_per_day
    end_prd = start_prd + shift_len
    if end_prd > prds_per_day:
        end_prd = - (prds_per_day - end_prd)

    start_hours = (start_prd - 1) * mins_per_prd / 60.0
    start_whole_hours = int(start_hours)
    start_whole_minutes = int((start_hours - start_whole_hours) * 60)
    start_time = time(start_whole_hours, start_whole_minutes).strftime('%H%M')

    end_hours = (end_prd - 1) * mins_per_prd / 60.0
    end_whole_hours = int(end_hours)
    end_whole_minutes = int((end_hours - end_whole_hours) * 60)
    try:
        end_time = time(end_whole_hours, end_whole_minutes).strftime('%H%M')
    except:
        pass

    label = '{}-{}'.format(start_time, end_time)

    return label

def make_tours(tur_df, prds_per_fte, nweeks):
    tour_df = tur_df.groupby('tournum').agg(tourtype=('tourtype', 'min'),
                                            tot_shifts=('tourtype', 'size'), tot_periods=('shiftlength', 'sum'),
                                            startwin=('startwin', 'min'))

    tour_df['tot_hours'] = tour_df['tot_periods'] * 40.0 / prds_per_fte
    tour_df['tot_ftes'] = tour_df['tot_periods'] / nweeks / prds_per_fte
    tour_df.reset_index(inplace=True, drop=False)

    tour_df = tour_df.astype(dtype={'tournum': np.int32,
                                'tourtype': np.int32,
                                'tot_shifts': np.int32,
                                'tot_periods': np.int32,
                                'startwin': np.int32})

    return tour_df


def make_tourtype_summary(tour_df):
    ttype_sum_df = tour_df.groupby('tourtype').agg(num_tours=('tourtype', 'count'),
                                                   tot_periods=('tot_periods', 'sum'),
                                                   tot_shifts=('tot_shifts', 'sum'),
                                                   tot_hours=('tot_hours', 'sum'),
                                                   tot_ftes=('tot_ftes', 'sum'), )

    ttype_sum_df.reset_index(inplace=True)
    return ttype_sum_df


def make_summary(tours_df, prds_per_fte, nweeks):
    # get_fte = lambda x: x.sum() / prds_per_fte / nweeks
    # https://stackoverflow.com/questions/38179212/custom-describe-or-aggregate-without-groupby/41363399
    summary_df = tours_df.groupby(lambda _: 0).agg(num_tours=('tourtype', 'count'),
                                                   tot_periods=('tot_periods', 'sum'),
                                                   tot_shifts=('tot_shifts', 'sum'),
                                                   tot_hours=('tot_hours', 'sum'),
                                                   tot_ftes=('tot_ftes', 'sum'))

    return summary_df


def write_phase1_capsummary(inst, isStringIO=True):
    """
    Write out a multiweek summary of capacity, demand, understaffing.

    :param inst: Model instance
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: capacity summary as StringIO object or a string.
    """
    param = 'period day week dmd cap us1 us2 ustot\n'
    rows = [(i, j, w,
             inst.dmd_staff[i, j, w],
             inst.cov[i, j, w].value,
             inst.under1[i, j, w].value,
             inst.under2[i, j, w].value,
             inst.under1[i, j, w].value + inst.under2[i, j, w].value)
            for i in inst.PERIODS
            for j in inst.DAYS
            for w in inst.WEEKS
            ]

    for row in rows:
        row = [str(r) for r in row]
        data_row = ' '.join(row)
        data_row += '\n'
        param += data_row

    if isStringIO:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def write_phase1_shiftsummary(inst, isStringIO=True):
    """
    Write out a multiweek summary of daily shift worked variables
    in multi-week format - hopefully useful for debugging.

    :param inst: Model instance
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: Shift summary as StringIO object or a string.
    """

    param = ''
    for t in inst.activeTT:
        for i in inst.PERIODS:
            if (i, t) in inst.okTourType:
                for k in inst.tt_length_x[t]:
                    data_row = '{}|{}|{}|{}|{}'.format(
                        t, i, k, inst.TourType[i, t], inst.TourType[i, t].value)
                    for w in inst.WEEKS:
                        for j in inst.DAYS:
                            if (i, t, k, j, w) in inst.TourTypeDayShift_idx:
                                data_row += '|{}|{}'.format(
                                    inst.TourTypeDayShift[i, t, k, j, w],
                                            inst.TourTypeDayShift[i, t, k, j, w].value)
                            else:
                                data_row += '|{}|{}'.format('No shift', 0)
                    if inst.TourType[i, t].value > 0:
                        data_row += '\n'
                        param += data_row

    if isStringIO:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def write_phase2_toursummary(fn_tur, prds_per_fte, nweeks, isStringIO=True):
    """
    Write out a summary of tours by tour type

    :param inst: Model instance
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: tour summary as StringIO object or a string.
    """

    colnames = ['tournum', 'tourtype', 'week', 'period', 'day', 'shiftlength', 'startwin']

    tur_df = pd.read_csv(fn_tur,
                         skiprows=3,
                         names=colnames,
                         sep='\s+',
                         dtype={'tournum': np.int32,
                                'tourtype': np.int32,
                                'week': np.int32,
                                'period': np.int32,
                                'day': np.int32,
                                'shiftlength': np.int32,
                                'startwin': np.int32})

    tours_df = make_tours(tur_df, prds_per_fte, nweeks)
    ttype_sum_df = make_tourtype_summary(tours_df)
    fte_sum_df = make_summary(tours_df, prds_per_fte, nweeks)

    ttype_sum = ttype_sum_df.to_string(index=False)
    fte_sum = fte_sum_df.to_string(index=False)

    if isStringIO:
        ttype_sum_out = io.StringIO()
        ttype_sum_out.write(ttype_sum)

        fte_sum_out = io.StringIO()
        fte_sum_out.write(fte_sum)

        return ttype_sum_out.getvalue(), fte_sum_out.getvalue()
    else:
        return ttype_sum, fte_sum


def create_mwt(fn_tur, prds_per_day, nweeks, prds_per_fte, scenario, output_path):

    """
    Write out tour schedule

    :param inst: Model instance
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: tour summary as StringIO object or a string.
    """

    phase2_mwt_file = output_path + scenario + '_phase2_mwt.csv'

    colnames = ['tournum', 'tourtype', 'week', 'period', 'day', 'shiftlength', 'startwin']

    tur_df = pd.read_csv(fn_tur,
                         skiprows=3,
                         names=colnames,
                         sep='\s+',
                         dtype={'tournum': np.int32,
                                'tourtype': np.int32,
                                'week': np.int32,
                                'period': np.int32,
                                'day': np.int32,
                                'shiftlength': np.int32,
                                'startwin': np.int32})

    tur_df['shift_label'] = tur_df.apply(lambda  x: shift_label(x['period'], x['shiftlength'], prds_per_day), axis=1)

    tours_df = make_tours(tur_df, prds_per_fte, nweeks)
    ntours = len(tours_df.index)

    mwtours = []
    for t in range(1, ntours + 1):
        mwtourspec = []
        for w in range(1, nweeks + 1):
            for _ in range(1, 8):
                mwtourspec.append('x')
        mwtours.append(mwtourspec)

    for index, row in tur_df.iterrows():
        tournum = row['tournum']
        week = row['week']
        dow = row['day']
        shift = row['shift_label']
        mwtours[tournum - 1][(week - 1) * 7 + dow - 1] = shift

    tour_shifts_df = pd.DataFrame(mwtours)
    # tour_shifts_df.index.rename('tournum', inplace=True)
    # tour_shifts_df.reset_index(inplace=True, drop=False)

    tour_schedule_df = tours_df.join(tour_shifts_df)
    tour_schedule_df.to_csv(phase2_mwt_file)

    pass


def create_mwt_old(fn_tur, prds_per_day, output_stub, output_path):
    pattLineType = re.compile('^(TTS|PP4|n_weeks|n_tours)')    # The regex to match the file section headers

    #pattTourShift = re.compile('^tourshift\[([0-9]+),([0-9]+),([0-9]+),([0-9]+),([0-9]+),([0-9]+),([0-9]+)')
    outFiles = {}
    whichFile = None


#===============================================================================
# PP4
#  1  2  1 11  4 16 11  4
#  1  2  1 11  6 16 11  6
#  1  2  2 11  3 16 11  3
#  1  2  2 11  4 16 11  4
#  1  2  2 11  5 16 11  5
#  1  2  2 11  6 16 11  6
#  2  3  1 11  1 16 11  1
#  2  3  1 11  3 24 11  3
#  2  3  1 11  5 24 11  5
#  2  3  1 11  6 24 11  6
#  2  3  2 11  1 24 11  1
#  2  3  2 11  2 24 11  2
#  2  3  2 11  6 24 11  6
#===============================================================================


    # outFiles['TTS'] = open(output_path+stubOutput+'.tts','w')
    outFiles['MWT'] = open(output_path + output_stub + '.mwt', 'w')
    # pp2key = []



    #outFiles['PP2'] = open(stubOutput+'.pp2','w')
    #outFiles['PP3'] = open(stubOutput+'.pp3','w')

    inFile = open(fn_tur)           # Open the file
    pp4mode = False
    toursInitialized = False
    num_weeks = 0
    num_tours = 0

    while 1:
        line = inFile.readline()           # Read a single line
        if not line: break                 # Bail if EOF

        matchobj = pattLineType.match(line)    # See if this line is start of new section
        if matchobj:
            if matchobj.group(1) == 'PP4':
                pp4mode = True
            else:
                pp4mode = False

            # If this is the num_weeks section, create the keys and files
            # for each week
            if matchobj.group(1) == 'n_weeks':
                num_weeks = int(line.split()[1])
                # for w in range(1,num_weeks+1):
                #     pp2key.append("pp2_%d" % (w))
                #     key = pp2key[w-1]
                #     outFiles[key] = open("%s%s_w%d.pp2" % (output_path, stubOutput, w),'w')

            if matchobj.group(1) == 'n_tours':
                num_tours = int(line.split()[1])

            # Initialize the tour data structure if not done yet
            if num_weeks>0 and num_tours>0 and not toursInitialized:
                tours = []
                mwtours = []
                tourtypes = []
                for w in range(1,num_weeks+1):
                    weeklytours = []
                    for t in range(1,num_tours+1):
                        tourspec = []

                        for _ in range(1,8):
                            tourspec.append(0)
                        tourspec.append(1)
                        for _ in range(1,8):
                            tourspec.append(0)
                        weeklytours.append(tourspec)
                    tours.append(weeklytours)

                for t in range(1,num_tours+1):
                    mwtourspec = []
                    for w in range(1,num_weeks+1):
                        for _ in range(1,8):
                            mwtourspec.append('x')
                    mwtours.append(mwtourspec)

                toursInitialized = True
                
        else:
            if pp4mode:
                tour = line.split()
                if len(tour) > 5:
                    tournum = int(tour[0])
                    ttype = int(tour[1])
                    week = int(tour[2])
                    startprd = int(tour[3])
                    dow = int(tour[4])
                    shiftlen = int(tour[5])


                    tours[week-1][tournum-1][dow-1] = startprd
                    tours[week-1][tournum-1][dow-1+8] = shiftlen

                    shift = shift_label(startprd, shiftlen, prds_per_day)
                    #'(' + str(startprd) + ':' + str(shiftlen) + ')'
                    mwtours[tournum-1][(week-1)*7+dow-1] = shift


                #whichFile.write(line)         # No match, just write the line to active output file
            else:
                # outFiles['TTS'].write(line)   # Write TTS line
                tourtypes.append(int(line))


    # for w in range(1,num_weeks+1):
    #     key = pp2key[w-1]
    #     for t in range(1,num_tours+1):
    #         for i in tours[w-1][t-1]:
    #              outFiles[key].write("%5d" % i)
    #         outFiles[key].write('\n')



    for t in range(1,num_tours+1):
        print (mwtours[t-1])
        # outFiles['MWT'].write(mwtours[t])
        # outFiles['MWT'].write(str(tourtypes[t-1]))
        json.dump(mwtours[t-1],outFiles['MWT'])
        outFiles['MWT'].write('\n')

    # Close all the files
    # for k in pp2key:
    #     outFiles[k].close

    outFiles['MWT'].close
    inFile.close

def weekenddaysworked_to_tourskeleton(inst, isStringIO=True):
    """
    Use values of WeekendDaysWorked variables to create a
    "weekend tour skeleton" which indicates days worked and days off.

    :param inst: Model instance
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: Tour skeleton as either as a StringIO object or a string.
    """

    header_list = ['n', 'i', 't', 'p']
    for w in range(1, inst.n_weeks + 1):
        for d in range(1, inst.n_days_per_week + 1):
            header_list.append(str(w) + '_' + str(d))

    param = ','.join(header_list) + '\n'

    tour_num = 0
    for (i, t, pattern) in inst.weekend_days_worked_idx:
        try:
            val = int(round(inst.WeekendDaysWorked[i, t, pattern]()))
        except TypeError:
            val = inst.WeekendDaysWorked[i, t, pattern]()

        if val > 0:
            pos_list = [str(p) for p in (i, t, pattern)]

            for tour in range(val):
                tour_num += 1

                data_row = str(tour_num) + ',' + ','.join(pos_list) + ','
                e = inst.weekend_type[i, t]
                day_list = []
                for w in range(1, inst.n_weeks + 1):
                    for d in range(1, inst.n_days_per_week + 1):
                        day_list.append(str(inst.A_wkend_days[pattern, d, w, t, e]))

                data_row += ','.join(day_list) + '\n'
                param += data_row

    param += '\n'

    if isStringIO:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param


def tourtypeday_to_tourskeleton(inst, isStringIO=True):
    """
    Use values of TourTypeDay variables to create a "tour skeleton"
    which indicates days worked and days off.

    :param inst: Model instance
    :param isStringIO: True (default) to return StringIO object, False to return string
    :return: Tour skeleton as either as a StringIO object or a string.
    """

    # 'dummy' entry in list is for column heading lineup in csv import
    header_list = ['n', 'i', 't', 'dummy']
    for w in range(1, inst.n_weeks + 1):
        for d in range(1, inst.n_days_per_week + 1):
            header_list.append(str(w) + '_' + str(d))

    param = ','.join(header_list) + '\n'

    tour_num = 0
    for (i, t) in inst.okTourType:
        if inst.TourType[i, t].value > 0:
            pos_list = [str(p) for p in (i, t)]
            tour_num += 1
            data_row = str(tour_num) + ',' + ','.join(pos_list) + ', x, '
            day_list = []
            for w in range(1, inst.n_weeks + 1):
                for d in range(1, inst.n_days_per_week + 1):
                    day_list.append(str(inst.TourTypeDay[i, t, d, w]()))

            data_row += ','.join(day_list) + '\n'
            param += data_row

    param += "\n"

    # build the header

    if isStringIO:
        param_out = io.StringIO()
        param_out.write(param)
        return param_out.getvalue()
    else:
        return param

def main():
    #extract_outfiles(sys.argv[1],sys.argv[2])
    pass

if __name__ == '__main__':
    #extract_outfiles(sys.argv[1],sys.argv[2])
    pass

