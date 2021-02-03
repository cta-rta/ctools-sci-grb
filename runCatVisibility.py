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
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
from lib.visibility import Visibility, complete_irf_name

# parse command line inputs
parser = argparse.ArgumentParser(description='This software runs the CTA visibility for a given configuration pass via YAML configuration file. The output is saved as NPY binary file.')
parser.add_argument('-f', '--config', required=True, type=str, help='configuration yaml file')
# configuration file
cf = parser.parse_args().config
# load params configuration from cf
with open(cf) as f:
    cfg = yaml.load(f, Loader=yaml.FullLoader)

logging.basicConfig(filename=__file__.replace('.py','.log'), filemode='w+', level=logging.DEBUG, format='%(asctime)s %(message)s')

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

logging.info('#################')
logging.info('# CONFIGURATION #')
logging.info(f'#################\n\n{yaml.dump(cfg)}')
logging.info('##############')
logging.info('# VISIBILITY #')
logging.info('##############')

# ----------------------------------------------------------------------------- loop runid

data = {}
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
            t_start = Time(hdr['GRBJD'] * u.day, format='jd')
        except KeyError:
            raise ValueError('This catalog cannot be processed. The headers do not contain a "GRBJD" trigger time keyword.')

        try:
            times = np.array(hdul['TIMES (AFTERGLOW)'].data.tolist())
        except KeyError:
            times = np.array(hdul['TIMES'].data.tolist())
        duration = Time(((times[-1] + times[1]) / 2)[0] / 86400, format='jd')

    # ignore warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        # --------------------------------------------------------------------------- loop sites
        for site in cfg['sites_list']:
            logging.info(f'{site} site')
            # initialise site coordinates
            site_coords = EarthLocation.of_site(cfg['sites_list'][site])
            # initialise
            visibility = Visibility()
            # visibility points in JD and AltAz
            visibility.visibility_points(t_start, duration, cfg['total_points'])
            visibility.visibility_altaz(source_radec, cfg['sites_list'][site])
            # find nights account for Moon (use default Moon thresholds)
            nights = visibility.get_nighttime_moonlight(twilight=cfg['setup']['twilight'], moon_sep=cfg['setup']['moon_sep'], fov_rad=cfg['setup']['fov_rad'], moonpha=0, max_moonpha=cfg['setup']['moon_pha'])
            del visibility
            # within each night find IRFs
            if type(nights['start']) == float:
                logging.info('................Not visible either to moonlight or daylight conditions')
                irfs = nights
                irfs['zref'] = np.nan
                #print(f'irfs: {irfs}')
                data[f'{runid.replace(".fits", "")}'][f'{site}'] = irfs
                continue
                         
            logging.info('Observability windows:') 
            for i in range(len(nights['start'])):
                logging.info(f'................Night {i+1} of {len(nights["start"])} in [{nights["start"][i]}, {nights["stop"][i]}]')
                t_start = Time(nights['start'][i], format='jd')
                duration = Time(nights['stop'][i] - nights['start'][i], format='jd')
                # initialise
                visibility = Visibility()
                # visibility points in JD and AltAz
                visibility.visibility_points(t_start, duration, cfg['window_points'])
                visibility.visibility_altaz(source_radec, cfg['sites_list'][site])
                # IRFs and relative time intervals
                irfs = visibility.associate_irf_zenith_angle(cfg['setup']['thresholds'], cfg['setup']['zenith_angles'])
                if type(irfs['zref']) == float:
                    logging.info('................Not visible due to low altitude')
                else:
                    logging.info('................Altitude intervals:')
                    for n in range(len(irfs['zref'])):
                        logging.info(f'................Zenith Ref. {irfs["zref"][n]} in [{irfs["start"][n]}, {irfs["stop"][n]}]')
                data[f'{runid.replace(".fits", "")}'][f'{site}'] = irfs
                del visibility

#print(data)

np.save(output, data)




