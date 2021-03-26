import numpy as np
import argparse

# parse command line inputs
parser = argparse.ArgumentParser(description='This script is a simple example on how to read a NPY binary file.')
parser.add_argument('-f', '--file', required=True, type=str, help='configuration yaml file')
# configuration file
filename = parser.parse_args().file

data = np.load(filename, allow_pickle=True, encoding='latin1', fix_imports=True).flat[0]

print(f'These are the parent keywords: \n\t{sorted(data.keys())}')
events = list(data.keys())
print(f'Inside each of the parent keyword you will find the following keys: \n\t{sorted(data[events[0]].keys())}')
sites = list(data[events[0]].keys())
print(f'Each key is a dictionary: \n\t{data[events[0]][sites[0]].keys()}')

print(f'In example you can access data with a nested cycle:')
for n, event in enumerate(events):
    for site in sites:
        print(f'Template {event} - Site {site} - Visibility Intervals:')
        if type(data[event][site]) == float:
            print(f'\tThis contains NaNs ---> the source is not observable due to daylight or moon.')
        else:
            print(f"\tThis contains the start and stop of each night.")
            for night in data[event][site]:
                print(f"\t{night}")
                print(f"\tstart: {data[event][site][night]['start']}")
                print(f"\tstop: {data[event][site][night]['stop']}")
                print(f"\n\tEach night contains the irfs intervals under the 'irfs' keyword.")
                if type(data[event][site][night]['irfs']['zref']) == float:
                    print(f'\tThis contains NaNs ---> the source is not observable at the site.')
                else:
                    for i in range(len(data[event][site][night]['irfs']['zref'])):
                        print(f'\t\tstart: {data[event][site][night]["irfs"]["start"][i]}')
                        print(f'\t\tstop: {data[event][site][night]["irfs"]["stop"][i]}')
                        print(f'\t\tzref: {data[event][site][night]["irfs"]["zref"][i]}')

    if n > 2:
        print("Let's stop after 3 templates.")
        break