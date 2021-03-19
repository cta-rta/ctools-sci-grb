import sys
import os
from os.path import join, isfile, isdir
import gammalib
import ctools
import cscripts
import gammapy
from astropy.io import fits
import numpy as np
import math as m
import astropy
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.utils import iers
from astropy.io import fits
import warnings
import argparse
import yaml
import logging
import numpy as np
import glob


parser = argparse.ArgumentParser(description='TThe significance of a GRB observation is computed at different IRFs according to a visibility tabele, created with runCatVisibility.py. A configuation YAML file is required, the output is saved as NPY binary file.')
parser.add_argument('-f', '--config', required=True, type=str, help='configuration yaml file')
# configuration file
cf = parser.parse_args().config
# load params configuration from cf
with open(cf) as f:
    cfg = yaml.load(f, Loader=yaml.FullLoader)


# ----------------------------------------------------------------------------- catalog
if '$' in cfg['path']['xmlcatalog']:
    catalog = os.path.expandvars(cfg['path']['xmlcatalog'])
else:
    catalog = cfg['path']['xmlcatalog']

if '$' in cfg['path']['output']:
    output = os.path.expandvars(cfg['path']['output'])
else:
    output = cfg['path']['output']

if cfg['path']['xmlfilename'] != None and '$' in cfg['path']['xmlfilename']:
    filename = os.path.expandvars(cfg['path']['xmlfilename'])
else:
    filename = cfg['path']['xmlfilename']

if not isdir(catalog):
    raise ValueError('Please correctly select a catalog folder')

if cfg['path']['xmlfilename'] == None:
    runids = glob.glob(cfg['path']['xmlcatalog'] + '/**/*.xml', recursive=True)
    if len(runids) == 0:
        raise ValueError('No valid XML file found')
else:
    runids = [filename]
runids = sorted(runids)


#--------------------------------------------------------------reading visibility table

data= np.load(output, allow_pickle=True, encoding='latin1', fix_imports=True).flat[0]
events = list(data.keys())
sites = list(data[events[0]].keys())
#data = table.copy()
#------------------------------#from lib.visibility import Visibility, complete_irf_name---------------------------------importing ctools variables

caldb=cfg['ctools']['caldb']
sim_rad = cfg['ctools']['rad']
sim_e_max=cfg['ctools']['Emax']
seed = cfg['ctools']['seed']
fitmodel=cfg['ctools']['fitmodel']


nn=0.
ss=0.

for n, event in enumerate(events):
    print(f'\nProcessing {event}')

    #mxl model for the event
    inmodel = runids[n]

    with fits.open(cfg['path']['catalog']+f'/{event}.fits') as hdul:
        hdr = hdul[0].header
        t_trigger = Time(hdr['GRBJD'] * u.day, format='jd')
        trigger=float(t_trigger.jd)

#-------------the following parameters can be useful for plots, otherwise are useless

        #z = hdr['z']
        #Eiso = hdr ['EISO']
    #adding reshift and Eiso to data dictionary
    #data[event]['z'] = z
    #data[event]['Eiso'] = Eiso


    for site in sites:
        print(f'\nProcessing site {site}')

        if type(data[event][site]['zref']) == float:
            print(f'\tThis contains NaNs event---> the source is not observable at the site.')

            data[event][site]['sigma_ON/OFF'] = np.nan
            if site == 'North':
                nn=nn+1
            else:
                ss=ss+1
        else:
            sigma= np.empty(shape= len(data[event][site]['zref']))

            for i in range(len(data[event][site]['zref'])):
                t_min = data[event][site]["start"][i]
                t_max= data[event][site]["stop"][i]
                name_irf = (f'{site}_z{data[event][site]["zref"][i]}_0.5h')


                #converting times from jd to seconds from trigger
                sim_t_min=(t_min- trigger)*86400
                sim_t_max=(t_max-trigger)*86400

                #selection of e_min according visibility tablesding to the irf
                if 'z60' in name_irf:
                    sim_e_min=0.101
                else:
                    sim_e_min=0.03


#-----------------------------------------------------------------------Simulation
                sim = ctools.ctobssim()
                sim['inmodel'] = inmodel
                sim['caldb'] = caldb
                sim['irf'] = name_irf
                sim['ra'] = 0.
                sim['dec'] = 0.
                sim['rad'] = sim_rad
                sim['tmin'] = sim_t_min
                sim['tmax'] = sim_t_max
                sim['emin'] = sim_e_min
                sim['emax'] = sim_e_max
                sim['seed'] = seed
                sim['outevents'] ='events_full_GRB.fits'
                sim.execute()


                onoff_sim = cscripts.csphagen()
                onoff_sim['inobs'] =  'events_full_GRB.fits'
                onoff_sim['inmodel'] = fitmodel
                onoff_sim['srcname'] = 'Crab'
                onoff_sim['ebinalg'] = 'LOG'
                onoff_sim['emin'] = sim_e_min
                onoff_sim['emax'] = sim_e_max
                onoff_sim['enumbins'] = 20
                onoff_sim['coordsys'] = 'CEL'

                onoff_sim['ra'] = 0.
                onoff_sim['dec'] = 0.5
                onoff_sim['rad'] = 0.2
                onoff_sim['caldb'] = caldb
                onoff_sim['irf'] = name_irf
                onoff_sim['bkgmethod'] = 'REFLECTED'
                onoff_sim['use_model_bkg'] = False
                onoff_sim['srcregfile'] =  'regioni_on.reg'
                onoff_sim['bkgregfile'] =  'regioni_off.reg'
                onoff_sim['outobs'] =  'GRBobs.xml'
                onoff_sim['outmodel'] =  'GRBmodel.xml'
                onoff_sim['stack'] = False
                onoff_sim.execute()

                a = 0
                with open( 'onoff_off.reg', 'r') as regioni_off:
                    for line in regioni_off:
                        if line.startswith('fk5'):
                            a += 1
                regioni_off.close()

                on = fits.open( 'onoff_pha_on.fits')
                tbdata_on = on[1].data
                conteggi_on = tbdata_on.field('counts')

                off = fits.open( 'onoff_pha_off.fits')
                tbdata_off = off[1].data
                conteggi_off = tbdata_off.field('counts')

                somma_on = 0
                somma_off = 0
                for valore_on in conteggi_on:
                    somma_on += valore_on

                for valore_off in conteggi_off:
                    somma_off += valore_off

                    on.close()
                    off.close()

                somma_off = somma_off/(a*1.0)
                a = 1.0/a

                valore = 2*(somma_on * m.log((1+a)/a*(somma_on/(somma_off+somma_on))) + somma_off * m.log((1+a)*(somma_off/(somma_off+somma_on))))
                if valore < 0:
                    valore = 0


                sigma[i] = str(m.sqrt(valore))
                print(f'\tSigma ON/OFF_{site}_z{data[event][site]["zref"][i]}: DONE')

            data[event][site]['sigma_ON/OFF'] = sigma

print (f"{nn} events can't be detected by CTA North")
print (f"{ss} events can't be detected by CTA South")
#print (data)
np.save(cfg['path']['sigmaoutput'] , data)
