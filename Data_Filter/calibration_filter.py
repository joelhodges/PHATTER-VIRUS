# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 14:42:45 2023

@author: joeli
"""

from astropy.table import Table
from astropy.time import Time, TimeDelta
import argparse as ap
import warnings

warnings.filterwarnings("ignore")

parser = ap.ArgumentParser(add_help=True)

parser.add_argument('input_file_full', type=str,
                    help='''name of the input file for the full dataset''')

parser.add_argument('input_file_m33', type=str,
                    help='''name of the input file for the M33 observations''')

parser.add_argument('output_file', type=str,
                    help='''name of the output file''')

parser.add_argument("label", type=str,
                    help='''Label to be filtered for, e.g. skyflat''')

parser.add_argument("input_row", type=int,
                    help='''row of the M33 input table to match''')

parser.add_argument("-dt", "--dawn_time",
                    help='''The time (in hh:mm:ss format) for the last observation time in the morning (defaults to 10 am)''',
                    type=str, default='10:00:00')

parser.add_argument("-tt", "--twilight_time",
                    help='''The time (in hh:mm:ss format) for the first observation time in the evening (defaults to 4 pm)''',
                    type=str, default='16:00:00')

args = parser.parse_args(args=None)

T = Table.read("./"+args.input_file_full, format='ascii.fixed_width_two_line')
m33_observations = Table.read("./"+args.input_file_m33, format='ascii.fixed_width_two_line')

# Find time difference between twilight and dawn
input_row = m33_observations[args.input_row]
time_range = TimeDelta(1)-((Time(str(input_row['Date'])[:11] + args.twilight_time)) - (Time(str(input_row['Date'])[:11] + args.dawn_time)))

def find_calibration_exposures(data_table, input_row, desired_label, date_constraint = 0):
    exposure_list = []
    dawn_cutoff = Time(str(input_row['Date'])[:11] + args.dawn_time, format='fits', out_subfmt='date_hms')
    
    # Adjust dawn to occur on the correct day of observation
    if(dawn_cutoff - Time(str(input_row['Date']), format='fits', out_subfmt='date_hms') < 0):
        dawn_cutoff = dawn_cutoff + TimeDelta(1)
       
    # Calculate twilight time
    twilight_cutoff = dawn_cutoff - time_range
    
    # Find frames corresponding to the correct label that fall within the correct time range
    for row in data_table:
        if desired_label in str(row['Description']):
            calibration_time = Time(str(row['Date']))
            
            if (twilight_cutoff - calibration_time < 0) and (dawn_cutoff - calibration_time > 0):
                exposure_list.append(row)
    return exposure_list

calibrations = find_calibration_exposures(T, input_row, args.label)
    
T = Table(rows = calibrations, names=['Exposure', 'Description', 'Date', 'Temp'])

T.write(args.output_file, format='ascii.fixed_width_two_line', overwrite=True)
