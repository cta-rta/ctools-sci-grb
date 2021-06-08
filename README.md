# ctools-sci-grb

This repository hosts the refactoring of the ctools pipeline for CTA-GRB-WG: https://github.com/thomasgas/ctools_pipe.

## Environment **

To create a virtual environment with all required dependencies:

```bash
conda env create --name <envname> --file=environment.yaml
```

Note that you should already have anaconda installed: https://www.anaconda.com/

## Calibration database

To complete the environment be sure to download and install the correct IRFs (only prod2 comes with ctools installation). Public ones can be found here: https://www.cta-observatory.org/science/cta-performance/


## Configuration file

Under cfg you can find a sample configuration file. Description of each parameter is commented within. This file will serve as input when running the code.

## Compiling the visibility table

After adjusting the configuration file to your needs, you can run the code as follows:

```bash
python runCatVisibility.py -f cfg/config.yaml
```

### Reading the visibility table

The output is saved via numpy as a binary NPY file. You can run an example of how to access data like this:

```bash
python readVisTable -f path/to/output.npy
```

### Notebook: plot the visibility
A notebook for useful plot and checks on the visibility is provided in the *notebooks* folder.


## GRB significance
After adjusting the configuration file to your needs, you can run the code as follows:

```bash
python GRB_significance.py -f cfg/config.yaml
```
Visibility tables are used as an input to compute significance and they are read directly from teh directory where you generated them with runCatVisibility.py so, please, do not move or rename them. 
As an input you will also need fits files for the events and xml files. These are generated with ctools_pipe (https://github.com/thomasgas/ctools_pipe). Be aware that the source must be simulated in (ra,dec)=(0, 0.75), in order for the counts of On and Off regions to be extracted correctly. 
The output is saved as a binary NPY file. A notebook (plotSignificance.ipynb in *notebooks* folder) is provided with some ready-to-plot example files (stored into the *notebooks/examples* directory). The notebook allows to display sigma evolution with time. Soon a script to display the npy table content will be provided.
<HR>
[**] subsceptible to changes 
