''' Produces a geomedian image from a defined epoch but limited
    to the months of interest e.g.November - March (nominal dry seasonon).
    Takes a list of tiles in TILELIST. 
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
    output_filename = outputdir + '/month_gm_2016-2017_'+'_'.join(label.split(','))+'.nc'
    print("Working on tile {}...".format(label))
else:
    x, y = (1385000.0, 1375000.0), (-4570000.0, -4580000.0)
    if subset:
        output_filename = 'month_gm_2016-2017_test_subset.nc'
    else:
        output_filename = 'month_gm_2016-2017_test_one.nc'

if os.path.exists(output_filename):
    print("output file already exists.")
    exit()

#####################################################

sensor = 'ls8' # make list to iterate over
deriv = 'nbart'
prod = sensor + '_'+deriv+'_albers'
prod_pq = sensor+'_pq_albers'

epoch = ('2016-12-01', '2017-01-30') # time query for datacube function can be just years
cmonths = [11,12,1,2,3,4] # a list of months for which you want results

####################################################
'''
# make a list of all the datasets (and pq datasets) available for area and epoch
scenes = dc.find_datasets(product=prod, time=epoch, **query)
pq_scenes = dc.find_datasets(product=prod_pq, time=epoch, **query)

# make new lists for just the months of interest
cm_ds = []
cm_pq_ds = []

for scene in scenes:
    if scene.center_time.month in cmonths:
        cm_ds.append(scene)
    else:
        print('No custom months found')

for pq_scene in pq_scenes:
    if pq_scene.center_time.month in cmonths:
        cm_pq_ds.append(pq_scene)
    else:
        print('No custom pq months found')
'''
####################################################

def load_ds(x, y):
    query = {'x': x,
             'y': y,
             'crs': 'EPSG:3577'}

    # make a list of all the datasets (and pq datasets) available for area and epoch
    scenes = dc.find_datasets(product=prod, time=epoch, **query)
    pq_scenes = dc.find_datasets(product=prod_pq, time=epoch, **query)

    # make new lists for just the months of interest
    cm_ds = []
    cm_pq_ds = []

    for scene in scenes:
        if scene.center_time.month in cmonths:
            cm_ds.append(scene)
        else:
            print('No custom months found')

    for pq_scene in pq_scenes:
        if pq_scene.center_time.month in cmonths:
            cm_pq_ds.append(pq_scene)
        else:
            print('No custom pq months found')


    ds = dc.load(product = prod,
                 datasets = cm_ds
                 group_by = 'solar_day',
                 dask_chunks={'time': 1},
                 **query)

    # Load PQ data for same query used to load Landsat data
    pq_ds = dc.load(product = prod_pq,
                    group_by = 'solar_day',
                    datasets = cm_pq_ds
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
    ds_mgm = GeoMedian().compute(ds)
    return ds_mgm.copy()

####################################################
xm, ym = (x[0]+x[1])/2, (y[0]+y[1])/2
x1, x2 = (x[0], xm), (xm, x[1])
y1, y2 = (y[0], ym), (ym, y[1])
if subset:
    out1 = load_ds(x1, y)
    out2 = load_ds(x2, y)
    out = xr.concat([out1, out2], dim='x')
else:
    out = load_ds(x, y)

# Output to netcdf
datacube.storage.storage.write_dataset_to_netcdf(out, output_filename)
