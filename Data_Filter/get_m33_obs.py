# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 13:43:12 2023

@author: joeli
"""

from astropy.table import Table
import argparse as ap
import warnings

warnings.filterwarnings("ignore")

parser = ap.ArgumentParser(add_help=True)

parser.add_argument('input_file', type=str,
                    help='''name of the input file''')

parser.add_argument('output_file_m33', type=str,
                    help='''name of the output file for M33 observations''')

parser.add_argument('output_file_dex', type=str,
                    help='''name of the output file for DEX observations''')

parser.add_argument("-tc", "--time_constraint",
                    help='''The time constraint of the DEX observations in hours (defaults to 2)''',
                    type=float, default=2)

args = parser.parse_args(args=None)

T = Table.read("./"+args.input_file, format='ascii.fixed_width_two_line')

def get_observations(data_table, label):
    output_list = []
    
    for row in data_table:
        if label in str(row['Description']):
            output_list.append(row)
    return output_list

def get_close_dex(m33_obs, dex_obs, time_constraint = args.time_constraint):
    output_list = []
    
    for dex in dex_obs:
        dex_time = str(dex['Date'])
        dex_time_hours = int(dex_time[11:13]) + int(dex_time[14:16])/60 + int(dex_time[17:19])/3600
        
        for m33 in m33_obs:
            m33_time = str(m33['Date'])
            m33_time_hours = int(m33_time[11:13]) + int(m33_time[14:16])/60 + int(m33_time[17:19])/3600
            
            # Checks if the two observations occurred on the same day
            if dex_time[:10] == m33_time[:10]:
                        
                # Checks if the two observations occured within the time constraint
                # Accurate to the second (but not tenth of second)
                if (abs(dex_time_hours - m33_time_hours)%24 < time_constraint) or (abs(dex_time_hours + m33_time_hours)%24 < time_constraint):
                    output_list.append(dex)
                    break
                        
    return output_list
            

m33_observations = get_observations(T, 'M33')
dex_observations = get_observations(T, 'DEX')

close_dex_observations = get_close_dex(m33_observations, dex_observations)

print(str(len(m33_observations))+' M33 observations found')
print(str(len(close_dex_observations))+' DEX observations found')
    
T = Table(rows = m33_observations, names=['Exposure', 'Description', 'Date', 'Temp'])
T.write(args.output_file_m33, format='ascii.fixed_width_two_line', overwrite=True)

T = Table(rows = close_dex_observations, names=['Exposure', 'Description', 'Date', 'Temp'])
T.write(args.output_file_dex, format='ascii.fixed_width_two_line', overwrite=True)
