''' Adapted from the burn mapping script. Will produce geomedian image
    for multiple Sentinel senors. Takes a list of tiles in TILELIST.
    Uses DEADataHandling to load the Sentinel scenes.
'''

import sys, os
import pandas as pd
import geopandas as gpd
import xarray as xr
import time
import multiprocessing
ncpus = multiprocessing.cpu_count()

from datacube_stats.statistics import GeoMedian
from datacube.helpers import ga_pq_fuser
from datacube.storage import masking
from datacube.helpers import write_geotiff

#get the DEA version of the plotting functions
sys.path.append(os.path.abspath('/g/data/r78/DPIPWE_lm/dea-notebooks/10_Scripts'))
import DEAPlotting
import DEADataHandling


import datacube
dc = datacube.Datacube(app='load_nbarx_geomedian')

###############################################################################
outputdir = '/g/data/r78/DPIPWE_lm/test_burn_mapping/output_data'
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
    output_filename = outputdir + '/sentgm_masked_2016_'+'_'.join(label.split(','))+'.nc'
    print("Working on tile {}...".format(label))
else:
    x, y = (1385000.0, 1375000.0), (-4570000.0, -4580000.0)
    if subset:
        output_filename = 'sentgm_xxxx_xxxx_test_subset.nc'
    else:
        output_filename = 'sentgm_xxxx_xxxx_test_one.nc'

if os.path.exists(output_filename):
    print("output file already exists.")
    exit()

###############################################################################
#custom parameters for query
#product = 'ard' #can specify other products?
time = ('2016-01-01', '2016-12-31')
resolution = (-10, 10)
bands = ['nbart_blue', 'nbart_green', 'nbart_red', 'nbart_nir_1']

####################################################

def multi_sentgm(x,y):
    query = {
            'x': x,
            'y': y,
            #'lat': (-35.27, -35.33),
            #'lon': (149.07, 149.15),
            'crs': 'EPSG:3577',
            'output_crs': 'EPSG:3577',
            'resolution': resolution, #(-20, 20),
            'time': time #('2017-01-01', '2017-02-01')
            }

    ds_sent = DEADataHandling.load_clearsentinel2(dc=dc, query=query,
                                                    masked_prop=0,
                                                    bands_of_interest = bands,
                                                    mask_pixel_quality=True)


    # compute geomedian
    dss_gm = GeoMedian().compute(ds_sent)
    return dss_gm.copy()

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
