''' Adapted from the burn mapping script. Will produce
    geomedian image for multiple Landsat senors.
    Takes a list of tiles in TILELIST
'''

import sys, os
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
dc = datacube.Datacube(app='load_nbarx_geomedian')

''' delete if not needed

from datacube_stats.statistics import GeoMedian

import datacube
from datacube.storage import masking
#from datacube.storage.masking import mask_to_dict
from datacube.storage.storage import write_dataset_to_netcdf
from datacube.helpers import ga_pq_fuser
dc = datacube.Datacube(app='dc_burnmap_comp')
'''

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
    output_filename = outputdir + '/multigm_2016-2017_'+'_'.join(label.split(','))+'.nc'
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
#datatime = ('2017-01-01', '2017-01-30') # period to retrieve data
#referenceperiod = ('2013-01-01', '2016-06-30') # period used for the calculation of geometric median
#mappingperiod = ('2016-07-01', '2017-06-30') # period of interest for change/severity mapping
#res = (25, 25)

product = 'nbart' #can be 'nbar', 'nbart' or 'fc'. Defaults to 'nbart'
sensors = ['ls5', 'ls7', 'ls8'] #take or remove as needed
query = {'x': x,
         'y': y,
         'time': ('2016-12-01', '2017-01-30'),
         'resolution': (25,25),
         'crs': 'EPSG:3577'}

####################################################

def multigm(x, y):
    dsm = DEADataHandling.load_clearlandsat(dc=dc, query=query,
                                           product=product,
                                           masked_prop=0,
                                           sensors = sensors,
                                           ls7_slc_off=True)

    
    ''' # Load PQ data for same query used to load Landsat data
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
    '''

    # compute geomedian
    out = GeoMedian().compute(dsm)
    return out.copy()

####################################################

xm, ym = (x[0]+x[1])/2, (y[0]+y[1])/2
x1, x2 = (x[0], xm), (xm, x[1])
y1, y2 = (y[0], ym), (ym, y[1])
if subset:
    out1 = burncomp(x1, y)
    out2 = burncomp(x2, y)
    out = xr.concat([out1, out2], dim='x')
else:
    out = burnmap(x, y)

# Output to netcdf
dc.storage.storage.write_dataset_to_netcdf(out, output_filename)
