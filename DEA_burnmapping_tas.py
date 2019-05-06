'''
Based on Fang's batch burn mapping script and adapted for Tas.
A TILELIST (text) file can be used to nominate the DEA tiles
for processing (e.g. tiles covering Tasmania). Incude functions
to show processing times for each step
'''

import pandas as pd
import geopandas as gpd
import xarray as xr
import sys, os
import time
import multiprocessing
ncpus = multiprocessing.cpu_count()

#set path to burncube module
sys.path.append(os.path.abspath('/g/data/r78/DPIPWE_lm/repos/burn-mapping/notebooks/handover'))
from BurnCube import BurnCube #including burn mapping main functions


# Map burnscars from time period specified

########################################################
# location identified from Australian Albers tiles

outputdir = '/g/data/r78/DPIPWE_lm/test_burn_mapping/output_data/2017_2018'
if not os.path.exists(outputdir):
    print("output directory doesn't exist")
    exit()

albers = gpd.read_file('/g/data/r78/DPIPWE_lm/test_burn_mapping/reference_data/Albers_Australia_Coast_Islands_Reefs.shp')

subset = True
label = None
if len(sys.argv)==2:
    label = sys.argv[1]
elif len(sys.argv)==3:
    label = "{},{}".format(sys.argv[1], sys.argv[2])

if label:
    index = albers[albers['label']==label].index[0]
    x = (albers.loc[index]['X_MIN'], albers.loc[index]['X_MAX'])
    y = (albers.loc[index]['Y_MIN'], albers.loc[index]['Y_MAX'])
    output_filename = outputdir + '/BurnScarMap_2017_2018_'+'_'.join(label.split(','))+'.nc'
    print("Working on tile {}...".format(label))
else:
    x, y = (1385000.0, 1375000.0), (-4570000.0, -4580000.0) # location used if TILELIST not provided
    if subset:
        output_filename = 'BurnScarMap_2017_2018_test_subset.nc'
    else:
        output_filename = 'BurnScarMap_2017_2018_test_one.nc'

if os.path.exists(output_filename):
    print("output file already exists.")
    exit()

########################################################

sensor = [7] # sensor to use
datatime = ('2014-01-01', '2018-06-30') # period to retrieve data
referenceperiod = ('2014-01-01', '2017-06-30') # period used for the calculation of geometric median (pre fire season)
mappingperiod = ('2017-07-01', '2018-06-30') # period of interest for change/severity mapping

def burnmap(x, y, sensor=sensor, datatime=datatime,
            referenceperiod=referenceperiod, mappingperiod=mappingperiod,
            method = "NBR", n_procs = ncpus, res = (25, 25)):
    bc = BurnCube()
    print("-- Using {} threads.---".format(n_procs))
#step1: load data and filtering
    start_time = time.monotonic()
    try:
        bc.load_cube(x, y, res, datatime, [sensor])
    except:
        print("Problem loading data")
        exit()
    print("---{} minutes for loading data.---".format((time.monotonic()-start_time)/60))
#step2: calculate geometric median
    start_time = time.monotonic()
    bc.geomedian(referenceperiod)
    print("---{} minutes for geomedian calculation.---".format((time.monotonic()-start_time)/60))
#step3: calculate cosine distance and nbr distance for reference period
    start_time = time.monotonic()
    bc.distances(referenceperiod, n_procs=n_procs)
    print("---{} minutes for cos dist.---".format((time.monotonic()-start_time)/60))
#step4: determine the threshold values
    start_time = time.monotonic()
    bc.outliers()
    print("---{} minutes for outliers calculation.---".format((time.monotonic()-start_time)/60))
#step5: calculate the distances to the reference
    start_time = time.monotonic()
    bc.distances(mappingperiod,n_procs=n_procs)
#step6: burn mapping for the given period
    start_time = time.monotonic()
    out = bc.severitymapping(mappingperiod, n_procs=n_procs, method=method, growing=True)
    print("---{} minutes for burn scar mapping.---".format((time.monotonic()-start_time)/60))
    return out.copy()

########################################################
# split into 2 subsets

xm, ym = (x[0]+x[1])/2, (y[0]+y[1])/2
x1, x2 = (x[0], xm), (xm, x[1])
y1, y2 = (y[0], ym), (ym, y[1])
if subset:
    out1 = burnmap(x1, y, method = 'NBRdist')
    out2 = burnmap(x2, y, method = 'NBRdist')
    out = xr.concat([out1, out2], dim='x')
else:
    out = burnmap(x, y, method ='NBRdist')

# save the output
out.to_netcdf(output_filename)
