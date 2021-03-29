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

if cfg['path']['output'] == None:
    output = 'visibility_output.npy'
elif '$' in cfg['path']['output']:
    output = os.path.expandvars(cfg['path']['output'])
else:
    output = cfg['path']['output']

if cfg['path']['xmlfilename'] == None:
    runids = glob.glob(cfg['path']['xmlcatalog'] + '/**/*.xml', recursive=True)
    if len(runids) == 0:
        raise ValueError('No valid XML file found')

elif type(cfg['path']['xmlfilename']) == str:
    if not isfile(join(catalog, cfg['path']['xmlfilename'])):
            raise ValueError(f'Specified template {runid} does not exist in catalog')
    runids = [cfg['path']['xmlfilename']]
else:
    runids = cfg['path']['xmlfilename']
    for runid in runids:
        if not isfile(join(catalog, runid)):
            raise ValueError(f'Specified template {runid} does not exist in catalog')
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
seeds = np.random.randint(1,1000,size=cfg['ctools']['iterations'])  #seeds randomly genrated to be used for the different simulations
fitmodel=cfg['ctools']['fitmodel']

for runid in runids:
    if type(cfg['path']['xmlfilename']) == str:
        inmodel = join(catalog, runid)
        src=f'{runid.replace(".xml", "")}'
        src=f'{src.replace("model_", "")}'

    elif cfg['path']['xmlfilename'] == None:
        inmodel = runid
        if len(runids) == len(data[events[0]]):
            for event in events:
                src = event
        else:
            print( "the possibility to run subcatalogs hasn't been implemented yet: you can either run one event per time or generate a visibility table for the subcataalog")
        break

    for event in events:
        if event == src:
            #print (event)
            print(f'\nProcessing {event}')
            with fits.open(cfg['path']['catalog']+f'/{event}.fits') as hdul:
                    hdr = hdul[0].header
                    t_trigger = Time(hdr['GRBJD'] * u.day, format='jd')
                    trigger=float(t_trigger.jd)
            for site in sites:
                    print(f'\nProcessing site {site}')
                    for night in data[event][site]:
                        print(f'\nProcessing {night}')
                        if type(data[event][site]) == float:
                            print(f'\tThis contains NaNs ---> the source is not observable due to daylight or moon.')
                        for i in range(len(data[event][site][night]['irfs']['zref'])):
                            if type(data[event][site][night]['irfs']['zref'][i]) == float:
                                print(f'\tThis contains NaNs event---> the source is not observable at the site.')
                                data[event][site][night]['sigma_ON/OFF'] = float(-9)
                            elif data[event][site][night]['irfs']['zref'][i] == -9.0:
                                print(f'\tThis contains NaNs event---> the source is not observable at the site.')
                                data[event][site][night]['sigma_ON/OFF'] = float(-9)

                            else:
                                somma_on=np.zeros(shape=len(data[event][site][night]['irfs']['zref']))
                                somma_off=np.zeros(shape=len(data[event][site][night]['irfs']['zref']))
                                valore=np.zeros(shape=len(seeds))
                                for j, seed in enumerate(seeds):
                                    seed = int(seed)
                                    print(f'seed number {j+1}: {seed}')
                                    night_start=data[event][site][night]["start"]
                                    for i in range(len(data[event][site][night]['irfs']['zref'])):
                                        t_min = data[event][site][night]['irfs']["start"][i]
                                        t_max= data[event][site][night]['irfs']["stop"][i]

                                        name_irf = (f'{site}_z{data[event][site][night]["irfs"]["zref"][i]}_0.5h')

                                        #converting times from jd to seconds from trigger

                                        sim_t_min=(t_min - trigger)*86400
                                        sim_t_max= (t_max-trigger)*86400

                                        print(f'\tsim_t_start (seconds from trigger) {site}-{night}-{name_irf}: {sim_t_min}')
                                        print(f'\tsim_t_stop (seconds from trigger){site}-{night}-{name_irf}:{sim_t_max}')

                                        #selection of e_min according visibility tablesding to the irf
                                        if 'z60' in name_irf:
                                            sim_e_min=0.101
                                        else:
                                            sim_e_min=0.03
                                        print (f'Irf : {name_irf} - Minimum energy: {sim_e_min}')

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
                                        sim['outevents'] =f'events_full_GRB.fits'
                                        sim.execute()

                                        print(f'Simulation site {site} - {night}: DONE')

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
                                        onoff_sim['irf'] =name_irf
                                        onoff_sim['bkgmethod'] = 'REFLECTED'
                                        onoff_sim['use_model_bkg'] = False
                                        onoff_sim['srcregfile'] =  'regioni_on.reg'
                                        onoff_sim['bkgregfile'] =  'regioni_off.reg'
                                        onoff_sim['outobs'] =  'GRBobs.xml'
                                        onoff_sim['outmodel'] =  'GRBmodel.xml'
                                        onoff_sim['stack'] =False
                                        onoff_sim.execute()

                                        a = 0
                                        with open( 'onoff_off.reg', 'r') as regioni_off:
                                            for line in regioni_off:
                                                if line.startswith('fk5'):
                                                    a += 1
                                        regioni_off.close()

                                        on = fits.open('onoff_pha_on.fits')
                                        tbdata_on = on[1].data
                                        conteggi_on = tbdata_on.field('counts')

                                        off = fits.open('onoff_pha_off.fits')
                                        tbdata_off = off[1].data
                                        conteggi_off = tbdata_off.field('counts')

                                        somma_on[i] = 0
                                        somma_off[i]= 0
                                        for valore_on in conteggi_on:
                                            somma_on[i] += valore_on

                                        for valore_off in conteggi_off:
                                            somma_off[i] += valore_off

                                            on.close()
                                            off.close()

                                        somma_off[i] = somma_off[i]/(a*1.0)
                                        a = 1.0/a
#------------------Adding up all the counts for the single night 
                                    conteggi_on=np.sum(somma_on)
                                    conteggi_off=np.sum(somma_off)

                                    print (f'Number of counts (per {night}) in the on region: {conteggi_on}')
                                    print (f'Number of counts (per {night}) in the off region: {conteggi_off}')
                                                    
 #----------------  Computing significance for the total counts of the night                              
                                    try:
                                        valore[j] = 2*(conteggi_on * m.log((1+a)/a*(conteggi_on/(conteggi_off+conteggi_on))) + conteggi_off * m.log((1+a)*(conteggi_off/(conteggi_off+conteggi_on))))
                                        if valore[j] < 0:
                                            valore[j] = 0
                                        valore[j]=np.sqrt(valore[j])
                                    except ValueError:
                                        continue
                                    print (f'significance site {site},{night}, seed {seed} :{valore[j]}')


                                sigma = np.mean(valore)
                                var = np.var(valore)
                                data[event][site][night]['sigma_ON/OFF'] = sigma
                                data[event][site][night]['sigma_var'] = var
                                print (f'mean significance, site {site}, {night}: {sigma}, variance: {var}')

    
    np.save(cfg['path']['sigmaoutput'] , data)
