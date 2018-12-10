''' Produces a geomedian image from a defined epoch but limited to
    the months of interest e.g.November - March (nominal dry seasonon).
    Takes a list of tiles in TILELIST. Uses the DEADataHandling
    to load the Landsat scenes.
'''

import sys, os
import pandas as pd
import geopandas as gpd
import xarray as xr
import time
import multiprocessing
ncpus = multiprocessing.cpu_count()
#from BurnCube import BurnCube #including burn mapping main functions
#bc = BurnCube()

from datacube_stats.statistics import GeoMedian
from datacube.helpers import ga_pq_fuser
from datacube.storage import masking
from datacube.helpers import write_geotiff

#get the DEA version of the plotting functions
sys.path.append(os.path.abspath('/g/data/r78/DPIPWE_lm/dea-notebooks/10_Scripts'))
import DEAPlotting
import DEADataHandling


import datacube
dc = datacube.Datacube(app='multi_landsat_geomedian')

outputdir = '/g/data/r78/DPIPWE_lm/test_burn_mapping/output_data'
if not os.path.exists(outputdir):
    print("output directory doesn't exist")
    exit()

#########################################################

subset = True
label = None
albers = gpd.read_file('/g/data/r78/DPIPWE_lm/test_burn_mapping/reference_data/Albers_Australia_Coast_Islands_Reefs.shp')

if len(sys.argv)==2:
    label = sys.argv[1]
elif len(sys.argv)==3:
    label = "{},{}".format(sys.argv[1], sys.argv[2])

if label:
    index = albers[albers['label']==label].index[0]
    x = (albers.loc[index]['X_MIN'], albers.loc[index]['X_MAX'])
    y = (albers.loc[index]['Y_MIN'], albers.loc[index]['Y_MAX'])
    output_filename = outputdir + '2015_2016/LS_gm_dry_2015_2016_'+'_'.join(label.split(','))+'.nc'
    print("Working on tile {}...".format(label))
else:
    x, y = (1385000.0, 1375000.0), (-4570000.0, -4580000.0)
    if subset:
        output_filename = 'multi_month_gm_2016-2017_test_subset.nc'
    else:
        output_filename = 'multi_month_gm_2016-2017_test_one.nc'

if os.path.exists(output_filename):
    print("output file already exists.")
    exit()

#####################################################
sens_list = ['ls8', 'ls7']
deriv = 'nbart'
resolution = (-25,25)
bands = ['red', 'green', 'blue', 'nir', 'swir1', 'swir2']

epoch = ('2015-11-01', '2016-04-30') # time query for datacube function can be just years
cmonths = [11,12,1,2,3,4] # a list of months for which you want results

####################################################
# function to return months of interest
def is_cm(month):
    return (month >= 11) | (month <= 4)

# function to load cube of data (for epoch),
# extract months of interst and compute the geomedian
def cm_multigm(x,y):
    query = {'x': x,
             'y': y,
             'time': epoch,
             'resolution': resolution,
             'crs': 'EPSG:3577'}

    ds_m = DEADataHandling.load_clearlandsat(dc=dc, query=query,
                                            product=deriv,
                                            masked_prop=0,
                                            sensors = sens_list, # maybe not needed
                                            bands_of_interest = bands,
                                            mask_pixel_quality=True,
                                            ls7_slc_off=True)


    # extract just the months of interest
    ds_cm = ds_m.sel(time=is_cm(ds_m['time.month']))


    # compute geomedian
    ds_cmgm = GeoMedian().compute(ds_cm)
    return ds_cmgm.copy()
####################################################
xm, ym = (x[0]+x[1])/2, (y[0]+y[1])/2
x1, x2 = (x[0], xm), (xm, x[1])
y1, y2 = (y[0], ym), (ym, y[1])
if subset:
    out1 = cm_multigm(x1, y)
    out2 = cm_multigm(x2, y)
    out = xr.concat([out1, out2], dim='x')
else:
    out = cm_multigm(x, y)

# Output to netcdf
datacube.storage.storage.write_dataset_to_netcdf(out, output_filename)
