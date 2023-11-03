#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:22:48 2023

@author: gregz
"""
import numpy as np
import argparse as ap
import logging
import warnings
import tables
import sys
from astropy.table import Table
sys.path.append("..") 

import matplotlib.pyplot as plt
import seaborn as sns
from astropy.time import Time, TimeDelta
from astropy.io import fits
from multiprocessing import Pool

from virusraw import VIRUSRaw
from virusobs import VIRUSObs

# =============================================================================
# Warning Suppression: good for clarity, bad for debugging
# =============================================================================
warnings.filterwarnings("ignore")

# =============================================================================
# Plotting hacks
# =============================================================================
sns.set_context('talk')
sns.set_style('ticks')
plt.rcParams["font.family"] = "Times New Roman"


def setup_logging(logname='input_utils'):
    '''Set up a logger for shuffle with a name ``input_utils``.

    Use a StreamHandler to write to stdout and set the level to DEBUG if
    verbose is set from the command line
    '''
    log = logging.getLogger('input_utils')
    if not len(log.handlers):
        fmt = '[%(levelname)s - %(asctime)s] %(message)s'
        fmt = logging.Formatter(fmt)

        level = logging.INFO

        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        handler.setLevel(level)

        log = logging.getLogger('input_utils')
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log

def find_calibration_exposures(Obj_row, Obj_Table, cal='Hg', 
                               time_constraint=2.):
    '''
    

    Parameters
    ----------
    Obj_row : astropy.table.Table row
        Row for the science observation
    Obj_Table : astropy.table.Table
        Table of the observations over a date range
    cal : str, optional
        Calibration or observation to link to science. The default is 'Hg'.
    time_constraint : float, optional
        Time window in hours. The default is 2..

    Returns
    -------
    output_list : TYPE
        DESCRIPTION.

    '''
    keys = list([str(t) for t in Obj_Table['Exposure']])
    values = list(Obj_Table['Description'])
    
    # time_obs = Time(str(Obj_row['Date']))
    date_list = [Time(str(row['Date'])) for row in Obj_Table]
    
    twi_beg, twi_end = get_twi_times(Obj_row)
    
    # time_diffs = [cal_time - time_obs for cal_time in date_list]
    
    cal_obs = []
    for key, value, time in zip(keys, values, date_list):
        if value == cal:
            if (time - twi_beg >=0) and (twi_end - time >= 0):
                cal_obs.append(key)
    
    
    return cal_obs

def get_twi_times(Obj_row):
    time_range = Time('2023-01-02T10:00:00') - Time('2023-01-01T16:00:00')
    
    twilight_end = Time(str(Obj_row['Date'])[:11] + '10:00:00', format='fits', out_subfmt='date_hms')
    
    # Adjust dawn to occur on the correct day of observation
    if(twilight_end - Time(str(Obj_row['Date']), format='fits', out_subfmt='date_hms') < 0):
        twilight_end = twilight_end + TimeDelta(1)
       
    # Calculate twilight time
    twilight_beginning = twilight_end - time_range
    
    return (twilight_beginning, twilight_end)

def reduce(key, h5table, basepath=None, ifuslots=None):
    date = key[:8]
    obs = int(key[8:15])
    exp = int(key[15:])
    virus = VIRUSRaw(date, obs, h5table, exposure_number=exp, basepath=basepath, ifuslots=ifuslots)
    return virus

def get_science(sciind):
    twi_obs = find_calibration_exposures(T[sciind], T, cal='skyflat', time_constraint=12.)
    twi_list = []
    for twio in twi_obs:
        twi_list.append(reduce(twio, h5table, basepath=basedir, ifuslots=ifuslots))
    CdA_obs = find_calibration_exposures(T[sciind], T, cal='CdA', time_constraint=12.)
    CdA_list = []
    for CdAo in CdA_obs:
        CdA_list.append(reduce(CdAo, h5table, basepath=basedir, ifuslots=ifuslots))
    Hg_obs = find_calibration_exposures(T[sciind], T, cal='Hg', time_constraint=12.)
    Hg_list = []
    for Hgo in Hg_obs:
        Hg_list.append(reduce(Hgo, h5table, basepath=basedir, ifuslots=ifuslots))
    Dark_obs = find_calibration_exposures(T[sciind], T, cal='dark', time_constraint=12.)
    Dark_list = []
    for Darko in Dark_obs:
        Dark_list.append(reduce(Darko, h5table, basepath=basedir, ifuslots=ifuslots))
    LDLS_obs = find_calibration_exposures(T[sciind], T, cal='ldls_long', time_constraint=12.)
    LDLS_list = []
    for LDLSo in LDLS_obs:
        LDLS_list.append(reduce(LDLSo, h5table, basepath=basedir, ifuslots=ifuslots))
    sky_obs = find_calibration_exposures(T[sciind], T, cal='dex', time_constraint=12.)
    sky_list = []
    for skyo in sky_obs:
        sky_list.append(reduce(skyo, h5table, basepath=basedir, ifuslots=ifuslots))
        
    VIRUS1 = reduce(keys[sciind], h5table, basepath=basedir, ifuslots=ifuslots)
    VIRUS2 = reduce(keys[sciind+1], h5table, basepath=basedir, ifuslots=ifuslots)
    VIRUS3 = reduce(keys[sciind+2], h5table, basepath=basedir, ifuslots=ifuslots)
    SCIENCE = VIRUSObs([VIRUS1, VIRUS2, VIRUS3],
                       arcRaw_list=Hg_list+CdA_list, DarkRaw_list=Dark_list,
                       twiRaw_list=twi_list, LDLSRaw_list=LDLS_list, dither_index=[0, 1, 2])
    return SCIENCE
    

parser = ap.ArgumentParser(add_help=True)

parser.add_argument('object_table', type=str,
                    help='''name of the output file''')

parser.add_argument('hdf5file', type=str,
                    help='''name of the hdf5 file''')

parser.add_argument('outname', type=str,
                    help='''name appended to shift output''')

parser.add_argument('target', type=str,
                    help='''name of the target for reduction''')

parser.add_argument('--basedir', type=str,
                    help='''name appended to shift output''',
                    default='/work/03946/hetdex/maverick')

if __name__ == '__main__':
    # =============================================================================
    # Get Arguments from parser and setup logging
    # =============================================================================
    args = parser.parse_args(args=None)
    args.log = setup_logging('get_object_table')
    
    # =============================================================================
    # Hard coding for the VIRUS rectified wavelengths
    # =============================================================================
    def_wave = np.linspace(3470, 5540, 1036)
    
    
    # =============================================================================
    # Load the *.h5 calibration table
    # =============================================================================
    basedir = args.basedir
    hdf5file = args.hdf5file
    
    h5file = tables.open_file(hdf5file, mode='r')
    h5table = h5file.root.Cals
    
    # =============================================================================
    # Get all of the ifuslot numbers from the calibration table
    # =============================================================================
    ifuslots = list(np.unique(['%03d' % i for i in h5table.cols.ifuslot[:]]))
    
    # =============================================================================
    # Get the table of observations
    # =============================================================================
    T = Table.read(args.object_table, format='ascii.fixed_width_two_line')
    
    keys = list([str(t) for t in T['Exposure']])
    values = list(T['Description'])
    
    target_name = args.target
    
    sci_obs = [key for key, value in zip(keys, values) if target_name in value.lower()]
    sci_inds = [j for key, value, j in zip(keys, values, np.arange(len(keys))) if target_name in value.lower()]
    sci_unique_obs, uinds = np.unique([sci_o[:-2] for sci_o in sci_obs], return_index=True)
    sci_unique_inds = np.array(sci_inds)[uinds]
    
    P = Pool(16)
    res = P.map(get_science, sci_unique_inds)
    P.close()
