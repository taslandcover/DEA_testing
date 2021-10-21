''' Adapted from the burn mapping script. Will produce
    geomedian image for multiple Landsat senors.
    Takes a list of tiles in TILELIST
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
dc = datacube.Datacube(app='nbarx_geomedian')

''' delete if not needed

from datacube_stats.statistics import GeoMedian

import datacube
from datacube.storage import masking
#from datacube.storage.masking import mask_to_dict
from datacube.storage.storage import write_dataset_to_netcdf
from datacube.helpers import ga_pq_fuser
dc = datacube.Datacube(app='dc_burnmap_comp')
'''

outputdir = '/g/data/r78/DPIPWE_lm/output_data'
if not os.path.exists(outputdir):
    print("output directory doesn't exist")
    exit()

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
    output_filename = outputdir + '/ls_multigm_2010_2015_'+'_'.join(label.split(','))+'.nc'
    print("Working on tile {}...".format(label))
else:
    x, y = (1385000.0, 1375000.0), (-4570000.0, -4580000.0)
    if subset:
        output_filename = 'multigm_2016-2017_test_subset.nc'
    else:
        output_filename = 'multigm_2016-2017_test_one.nc'

if os.path.exists(output_filename):
    print("output file already exists.")
    exit()

#####################################################

product = 'nbart' #can be 'nbar', 'nbart' or 'fc'. Defaults to 'nbart'
sensors = ['ls5', 'ls7', 'ls8'] #take or remove as needed
time = ('2010-01-01', '2015-12-31')
resolution = (25,25)

####################################################

def multigm(x, y):
    query = {'x': x,
             'y': y,
             'time': time,
             'resolution': resolution,
             'crs': 'EPSG:3577'}
    dsm = DEADataHandling.load_clearlandsat(dc=dc, query=query,
                                           product=product,
                                           masked_prop=0,
                                           sensors = sensors,
                                           bands_of_interest = ['blue', 'green', 'nir', 'red', 'swir1', 'swir2'],
                                           mask_pixel_quality=True,
                                           ls7_slc_off=True)


    # compute geomedian
    dsm_gm = GeoMedian().compute(dsm)
    return dsm_gm.copy()

####################################################

xm, ym = (x[0]+x[1])/2, (y[0]+y[1])/2
x1, x2 = (x[0], xm), (xm, x[1])
y1, y2 = (y[0], ym), (ym, y[1])
if subset:
    out1 = multigm(x1, y)
    out2 = multigm(x2, y)
    out = xr.concat([out1, out2], dim='x')
else:
    out = multigm(x, y)

# Output to netcdf
datacube.storage.storage.write_dataset_to_netcdf(out, output_filename)
#out.to_netcdf(output_filename)
