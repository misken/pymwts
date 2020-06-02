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


def write_phase2_tours(phase2_inst, prds_per_fte, tot_dmd, scenario, output_path):

    """
    Write out tour schedule files

    :param phase2_inst: Phase 2 model instance
    :param prds_per_fte: Number of periods which is equal to 1 FTE
    :param scenario: string used in filenames
    :output_path: output files get written here
    :return: nothing
    """

    phase2_mwt_file = output_path + scenario + '_phase2_mwt.csv'
    phase2_tur_file = output_path + scenario + '_phase2_tur.csv'
    phase2_tourtypesum_file = output_path + scenario + '_phase2_tourtypesum.csv'
    phase2_ftesum_file = output_path + scenario + '_phase2_ftesum.csv'

    # Create the tur file
    idx_copy = []
    for idx in phase2_inst.TourShift:
        idx_copy.append(idx)
    idx_sorted = sorted(idx_copy)

    prds_per_day = phase2_inst.n_prds_per_day.value
    n_weeks = phase2_inst.n_weeks.value
    colnames = ['tournum', 'tourtype', 'week', 'period', 'day', 'shiftlength', 'startwin']
    header = ' '.join(colnames)

    with open(phase2_tur_file, "w") as f2_tour:
        f2_tour.write('{}\n'.format(header))
        for idx in idx_sorted:
            if phase2_inst.TourShift[idx].value > 0:
                # tournum, tt, week, prd, day, shiftlenprds, startwin
                f2_tour.write(
                    '{} {} {} {} {} {} {}\n'.format(idx[0], idx[5], idx[3], idx[1], idx[2],
                                                    phase2_inst.lengths[idx[4]],
                                                    phase2_inst.WIN_x[idx[0]]))

    tur_df = pd.read_csv(phase2_tur_file,
                         header=0,
                         sep='\s+',
                         dtype={'tournum': np.int32,
                                'tourtype': np.int32,
                                'week': np.int32,
                                'period': np.int32,
                                'day': np.int32,
                                'shiftlength': np.int32,
                                'startwin': np.int32})

    tur_df['shift_label'] = tur_df.apply(lambda  x: shift_label(x['period'], x['shiftlength'], prds_per_day), axis=1)

    tours_df = make_tours(tur_df, prds_per_fte, n_weeks)
    ntours = len(tours_df.index)

    mwtours = []
    for t in range(1, ntours + 1):
        mwtourspec = []
        for w in range(1, n_weeks + 1):
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

    # Create the column headers
    dow_abbrevs = ['Su', 'Mo', 'Tu', 'We', 'Th', 'Fr', 'Sa']
    header = ['{}-{}'.format(d, w) for w in range(1, n_weeks + 1) for d in dow_abbrevs]
    tour_shifts_df.set_axis(header, axis=1,inplace=True)

    tour_schedule_df = tours_df.join(tour_shifts_df)
    tour_schedule_df.to_csv(phase2_mwt_file, index=False)

    # Create the tour summary files
    ttype_sum_df = make_tourtype_summary(tours_df)
    fte_sum_df = make_summary(tours_df, prds_per_fte, n_weeks)
    fte_sum_df['tot_dmd'] = tot_dmd
    fte_sum_df['sched_eff'] = tot_dmd / fte_sum_df['tot_periods']

    ttype_sum_df.to_csv(phase2_tourtypesum_file, index=False)
    fte_sum_df.to_csv(phase2_ftesum_file, index=False)


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

