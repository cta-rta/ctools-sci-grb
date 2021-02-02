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
parser = argparse.ArgumentParser(description='Simple test for source visibility.')
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
if '$' in cfg['path']['output']:
    output = os.path.expandvars(cfg['path']['output'])
else:
    output = cfg['path']['output']
if '$' in cfg['path']['filename']:
    filename = os.path.expandvars(cfg['path']['filename'])
else:
    filename = cfg['path']['filename']

logging.info(f'catalog: {catalog}')
logging.info(f'output: {output}')
if not isdir(catalog):
    raise ValueError('Please correctly select a catalog folder')

if cfg['path']['filename'] == 'null':
    runids = [f for f in os.listdir(catalog) if '.fits' in f and isfile(join(catalog, f))]
else:
    runids = [filename]

runs = {'north': dict(), 'south': dict()}
data = {'run': 'name', 'visibility': runs}
data = np.array(data.items())
print(data)


# ----------------------------------------------------------------------------- loop runid

for runid in runids:
    logging.info('------------------------------------------------------------------ #')
    logging.info(f'Processing {runid}')
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
            logging.info(f'Processing {site} site visibility')
            # initialise site coordinates
            site_coords = EarthLocation.of_site(cfg['sites_list'][site])
            # initialise
            visibility = Visibility()
            # visibility points in JD and AltAz
            visibility.visibility_points(t_start, duration, cfg['total_points'])
            visibility.visibility_altaz(source_radec, cfg['sites_list'][site])
            # find nights account for Moon (use default Moon thresholds)
            nights = visibility.get_nighttime_moonlight(twilight=cfg['setup']['twilight'], moon_sep=cfg['setup']['moon_sep'], fov_rad=cfg['setup']['fov_rad'])
            print(f'nights: {nights}')
            del visibility
            # within each night find IRFs
            if type(nights['start']) == type(float()):
                logging.info('Not visible')
                irfs = nights
                irfs['zref'] = np.nan
                print(f'irfs: {irfs}')
                break
            for i in range(len(nights['start'])):
                logging.info(f'Night {i+1} of {len(nights["start"])}')
                t_start = Time(nights['start'][i], format='jd')
                duration = Time(nights['stop'][i] - nights['start'][i], format='jd')
                # initialise
                visibility = Visibility()
                # visibility points in JD and AltAz
                visibility.visibility_points(t_start, duration, cfg['window_points'])
                visibility.visibility_altaz(source_radec, cfg['sites_list'][site])
                # IRFs and relative time intervals
                irfs = visibility.associate_irf_zenith_angle(cfg['setup']['thresholds'], cfg['setup']['zenith_angles'])
                print(f'irfs {irfs}')
                del visibility





