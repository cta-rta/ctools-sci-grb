# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import argparse
import warnings
import yaml
import os
import logging
from os.path import join, isfile, isdir
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
from lib.visibility_v2 import Visibility, complete_irf_name
from astropy.coordinates import solar_system_ephemeris
#solar_system_ephemeris.set('jpl') 

# parse command line inputs
parser = argparse.ArgumentParser(description='This software runs the CTA visibility for a given configuration pass via YAML configuration file. The output is saved as NPY binary file.')
parser.add_argument('-f', '--config', required=True, type=str, help='configuration yaml file')
# configuration file
cf = parser.parse_args().config
# load params configuration from cf
with open(cf) as f:
    cfg = yaml.load(f, Loader=yaml.FullLoader)

# ----------------------------------------------------------------------------- catalog
if '$' in cfg['path']['catalog']:
    catalog = os.path.expandvars(cfg['path']['catalog'])
else:
    catalog = cfg['path']['catalog']
if cfg['path']['output'] == None:
    output = 'visibility_output.npy'
elif '$' in cfg['path']['output']:
    output = os.path.expandvars(cfg['path']['output'])
else:
    output = cfg['path']['output']

if not isdir(catalog):
    raise ValueError('Please correctly select a catalog folder')    

if cfg['path']['filename'] == None:
    runids = [f for f in os.listdir(catalog) if '.fits' in f and isfile(join(catalog, f))]
    if len(runids) == 0:
        raise ValueError('No valid FITS file found')    
elif type(cfg['path']['filename']) == str:
    if not isfile(join(catalog, cfg['path']['filename'])):
            raise ValueError(f'Specified template {runid} does not exist in catalog')
    runids = [cfg['path']['filename']]
else:
    runids = cfg['path']['filename']
    for runid in runids:
        if not isfile(join(catalog, runid)):
            raise ValueError(f'Specified template {runid} does not exist in catalog')
runids = sorted(runids)

# -------------------------------------------------------------------------log the configuration

logname = output.replace('.npy','.log')
logging.basicConfig(filename=logname, filemode='w+', level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info('#################')
logging.info('# CONFIGURATION #')
logging.info(f'#################\n\n{yaml.dump(cfg)}')
logging.info('##############')
logging.info('# VISIBILITY #')
logging.info('##############')

# ----------------------------------------------------------------------------- loop runid

data = {}
# ignore warnings
with warnings.catch_warnings():
    warnings.filterwarnings('ignore')
    for runid in runids:
        logging.info('------------------------------------------------------------------ #')
        print(f'Processing {runid}')
        logging.info(f'Processing {runid}')
        logging.info('----------')
        data[f'{runid.replace(".fits", "")}'] = {}
        # load template
        with fits.open(join(catalog, runid)) as hdul:
            hdr = hdul[0].header
            # source coordinates
            source_radec = SkyCoord(ra=hdr['RA'] * u.deg, dec=hdr['DEC'] * u.deg, frame='icrs')
            # source trigger time and afterglow duration
            try:
                trigger = Time(hdr['GRBJD'] * u.day, format='jd')
            except KeyError:
                raise ValueError('This catalog cannot be processed. The headers do not contain a "GRBJD" trigger time keyword.')

            try:
                times = np.array(hdul['TIMES (AFTERGLOW)'].data.tolist())
            except KeyError:
                times = np.array(hdul['TIMES'].data.tolist())
            afterglow_duration = Time((times[-1] - times[1])[0] / 86400, format='jd')

        # --------------------------------------------------------------------------- loop sites
        for site in cfg['sites_list']:
            logging.info(f'{site} site')
            # initialise
            visibility = Visibility()
            # visibility points in JD and AltAz
            visibility.visibility_points(trigger, afterglow_duration, cfg['total_points'])
            visibility.visibility_altaz(source_radec, cfg['sites_list'][site])
            # find nights account for Moon (use default Moon thresholds)
            if cfg['setup']['moon_sep'] == None:
                nights = visibility.get_nighttime(twilight=cfg['setup']['twilight'])
            else: 
                nights = visibility.get_nighttime_moonlight(twilight=cfg['setup']['twilight'], moon_sep=cfg['setup']['moon_sep'], fov_rad=cfg['setup']['fov_rad'], moonpha=0, max_moonpha=cfg['setup']['moon_pha'])
            #del visibility
            # within each night find IRFs
            #print(nights)
            if len(nights['start']) == 1 and nights['start'] < 0:
                logging.info('................Not visible either to moonlight or daylight conditions')
                #del nights, irfs, site_coords
                   
            logging.info('Observability windows:') 
            data[f'{runid.replace(".fits", "")}'][f'{site}'] = {}
            for i in range(len(nights['start'])):
                #print('Night', i+1, 'of', len(nights['start']))
                logging.info(f'................Night {i+1} of {len(nights["start"])} in [{nights["start"][i]}, {nights["stop"][i]}]')
                data[f'{runid.replace(".fits", "")}'][f'{site}'][f'night{i+1:02d}'] = {'start': nights["start"][i], 'stop': nights["stop"][i]}
                #print(nights['start'][i], nights['stop'][i])
                t_start = Time(nights['start'][i], format='jd')
                night_duration = Time(nights['stop'][i] - nights['start'][i], format='jd')
                # initialise
                visibility = Visibility()
                # visibility points in JD and AltAz
                visibility.visibility_points(t_start, night_duration, cfg['window_points'])
                visibility.visibility_altaz(source_radec, cfg['sites_list'][site])
                # IRFs and relative time intervals
                irfs = visibility.associate_irf_zenith_angle(cfg['setup']['thresholds'], cfg['setup']['zenith_angles'])
                if type(irfs['start']) == float and irfs['start'] < 0:
                    logging.info('................Not visible due to low altitude')
                    data[f'{runid.replace(".fits", "")}'][f'{site}'][f'night{i+1:02d}']['irfs'] = irfs
                else:
                    logging.info('................Altitude intervals:')
                    for n in range(len(irfs['zref'])):
                        logging.info(f'................Zenith Ref. {irfs["zref"][n]} in [{irfs["start"][n]}, {irfs["stop"][n]}]')
                    data[f'{runid.replace(".fits", "")}'][f'{site}'][f'night{i+1:02d}']['irfs'] = irfs
                #del visibility
                #print(irfs)
                #print(irfs['start'][0], irfs['stop'][-1])
        #del nights, irfs, site_coords, night_duration
        #del afterglow_duration
np.save(output, data)


print("\n\nExit\n\n")

#os.system(f"cat {logname}")