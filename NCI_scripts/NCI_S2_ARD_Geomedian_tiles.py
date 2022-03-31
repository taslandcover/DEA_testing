'''
Sentinel 2 Geomedian generation using the DEA
tiles processing method.
Takes a list of tiles in TILELIST.
'''

#What can be deleted here?
#import warnings; warnings.simplefilter('ignore')
import datacube
#import fnmatch
import os
import pandas as pd
import geopandas as gpd

#get the DEA version of the plotting functions
import sys
#sys.path.insert(1, '../Tools/')
#from dea_tools.datahandling import load_ard
sys.path.append(os.path.abspath('/g/data/r78/DPIPWE_lm/repos/dea-notebooks/Scripts'))
from dea_datahandling import load_ard
from datacube_stats.statistics import GeoMedian
from datacube.utils.cog import write_cog
from datacube.drivers.netcdf import write_dataset_to_netcdf
import xarray as xr

#Specify output directory
outputdir = '/g/data/r78/DPIPWE_lm/output_data/'
#outputdir = './'
if not os.path.exists(outputdir):
    print("output directory doesn't exist")
    exit()
    
# Connect to datacube containing Sentinel 2 data
dc = datacube.Datacube(app='load_ard_and_geomedian')

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
    #output_filename = outputdir + '/month_gm_2016-2017_'+'_'.join(label.split(','))+'.nc'
    output_filename = outputdir + '/month_gm_2016-2017_'+'_'.join(label.split(','))
    print("Working on tile {}...".format(label))
else:
    x, y = (1385000.0, 1375000.0), (-4570000.0, -4580000.0)
    if subset:
        #output_filename = 'S2_ARD_gm_2021_test_subset.nc'
        output_filename = 'S2_ARD_gm_2021_test_subset'
    else:
        #output_filename = 'S2_ARD_gm_2021_test_one.nc'
        output_filename = 'S2_ARD_gm_2021_test_one'

if os.path.exists(output_filename):
    print("output file already exists.")
    exit()

#####################################################

def load_ds(x, y):
    query = {'x': x,
             'y': y,
             'crs': 'EPSG:3577',
             'time': ('2021'),
             'measurements': ['nbart_blue', 'nbart_green', 'nbart_red', 'nbart_nir_1'], # Can add nbart_swir2 for true flase colour but change res to 20
             'resolution': (-20, 20),
             'group_by': 'solar_day',
             'output_crs': 'EPSG:3577'}
        
    # Load available data from both Sentinel 2 satellites
    ds = load_ard(dc=dc,
                  products=['s2a_ard_granule', 's2b_ard_granule'],
                  dask_chunks={'time':1},
                  **query)

    '''
    # function to return months of interest (cm = composite month)
    def is_cm(month):
        return (month >= 11) | (month <= 4)

    # extract just the months of interest
    ds_cm = ds.sel(time=is_cm(ds['time.month']))

    ds_cm['time.month'] #take a look at what months we have...
    '''

    '''
    Alternate method for extracting months via list
    c_months = [11,12,1,2,3]
    ds_cm = ds.sel(time=(ds['time.month'].isin(c_months)).dropna(dim='time'))
    '''

    # Compute geomedian here is necessary - either for dataset or subset months
    ds_gm = GeoMedian().compute(ds)
    return ds_gm.copy()
  
#####################################################
xm, ym = (x[0]+x[1])/2, (y[0]+y[1])/2
x1, x2 = (x[0], xm), (xm, x[1])
y1, y2 = (y[0], ym), (ym, y[1])
if subset:
    out1 = load_ds(x1, y)
    out2 = load_ds(x2, y)
    out = xr.concat([out1, out2], dim='x')
else:
    out = load_ds(x, y)

 
''' 
# Output to netcdf
datacube.storage.storage.write_dataset_to_netcdf(out, output_filename)
'''

# Here we can export the geomedian
# for COG we need an array not a dataset
out_da = out.to_array()

# Write multi-band GeoTIFF to a location
write_cog(geo_im=out_da,
          fname=output_filename,
          overwrite=True)
