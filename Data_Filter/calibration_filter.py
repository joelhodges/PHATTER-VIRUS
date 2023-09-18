# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 14:42:45 2023

@author: joeli
"""

from astropy.table import Table
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

args = parser.parse_args(args=None)

T = Table.read("./"+args.input_file_full, format='ascii.fixed_width_two_line')
m33_observations = Table.read("./"+args.input_file_m33, format='ascii.fixed_width_two_line')


def find_calibration_exposures(data_table, input_row, desired_label, date_constraint = 0):
    exposure_list = []
    
    for row in data_table:
        if desired_label in str(row['Description']):
            if str(row['Date'])[:10] == str(input_row['Date'])[:10]:
                exposure_list.append(row)
    return exposure_list

calibrations = find_calibration_exposures(T, m33_observations[args.input_row], args.label)
    
T = Table(rows = calibrations, names=['Exposure', 'Description', 'Date', 'Temp'])

T.write(args.output_file, format='ascii.fixed_width_two_line', overwrite=True)
