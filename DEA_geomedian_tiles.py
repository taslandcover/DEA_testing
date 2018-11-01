''' Adapted from the burn mapping script but only to produce the
    geomedian image. Takes a list of tiles in TILELIST
'''

import pandas as pd
import geopandas as gpd
import xarray as xr
import sys, os
import time
import multiprocessing
ncpus = multiprocessing.cpu_count()
#from BurnCube import BurnCube #including burn mapping main functions
#bc = BurnCube()

from datacube_stats.statistics import GeoMedian

import datacube
from datacube.storage import masking
#from datacube.storage.masking import mask_to_dict
from datacube.storage.storage import write_dataset_to_netcdf
from datacube.helpers import ga_pq_fuser
dc = datacube.Datacube(app='dc_burnmap_comp')

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
    output_filename = outputdir + '/composite_2016-2017_'+'_'.join(label.split(','))+'.nc'
    print("Working on tile {}...".format(label))
else:
    x, y = (1385000.0, 1375000.0), (-4570000.0, -4580000.0)
    if subset:
        output_filename = 'composite_2016-2017_test_subset.nc'
    else:
        output_filename = 'composite_2016-2017_test_one.nc'
        
if os.path.exists(output_filename):
    print("output file already exists.")
    exit() 

#####################################################

# Edit these variables as needed

sensor = 'ls8'
time = ('2016-07-01', '2017-06-30') # period of interest for change/severity mapping
resolution = (25, 25)

####################################################

def gm_comp(x, y):
    query = {'x': x,
             'y': y,
             'time': time,
             'resolution': resolution,
             'crs': 'EPSG:3577'}
    
    ds = dc.load(product=sensor+'_nbart_albers',
                 group_by = 'solar_day',
                 dask_chunks={'time': 1},
                 **query)

    # Load PQ data for same query used to load Landsat data
    pq_ds = dc.load(product = sensor+'_pq_albers',
                    group_by = 'solar_day',
                    fuse_func=ga_pq_fuser,
                    dask_chunks={'time': 1},
                    **query)

    # Use PQ to create mask that is True for pixels that are not affected by clouds, cloud shadow or saturation
    good_quality_ds = masking.make_mask(pq_ds.pixelquality,
                                    cloud_acca='no_cloud',
                                    cloud_fmask='no_cloud',
                                    cloud_shadow_acca='no_cloud_shadow',
                                    cloud_shadow_fmask='no_cloud_shadow',
                                    blue_saturated=False,
                                    green_saturated=False,
                                    red_saturated=False,
                                    nir_saturated=False,
                                    swir1_saturated=False,
                                    swir2_saturated=False,
                                    contiguous=True)

    # Remove -999 nodata values prior to analysing or plotting Landsat imagery by setting all nodata values to `NaN`
    ds = masking.mask_invalid_data(ds)

    # Apply the mask to preserve only the good data
    ds = ds.where(good_quality_ds)

    # compute geomedian
    ds_gm = GeoMedian().compute(ds)
    return ds_gm.copy()

####################################################
xm, ym = (x[0]+x[1])/2, (y[0]+y[1])/2
x1, x2 = (x[0], xm), (xm, x[1])
y1, y2 = (y[0], ym), (ym, y[1])
if subset:
    out1 = gm_comp(x1, y)
    out2 = gm_comp(x2, y)
    out = xr.concat([out1, out2], dim='x')
else:
    out = gm_comp(x, y)

# Output to netcdf
datacube.storage.storage.write_dataset_to_netcdf(out, output_filename)
