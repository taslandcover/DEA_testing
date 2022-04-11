import datacube
import os
#import pandas as pd
import geopandas as gpd 
from odc.algo import to_f32, xr_geomedian, int_geomedian

#get the DEA version of the plotting functions
import sys
#sys.path.insert(1, '../Tools/')
sys.path.append(os.path.abspath('/g/data/r78/DPIPWE_lm/repos/dea-notebooks/Tools'))
#from dea_datahandling import load_ard
from dea_tools.datahandling import load_ard
#from dea_tools.dask import create_local_dask_cluster
from datacube_stats.statistics import GeoMedian
from datacube.utils.cog import write_cog
#from datacube.drivers.netcdf import write_dataset_to_netcdf
#import xarray as xr

outputdir = '/g/data/r78/DPIPWE_lm/output_data/'
if not os.path.exists(outputdir):
    print("output directory doesn't exist")
    exit()

# Connect to datacube containing Sentinel 2 data
dc = datacube.Datacube(app='Dask_load_ard_and_geomedian')

subset = True
label = '12,-47'
albers = gpd.read_file('/g/data/r78/DPIPWE_lm/test_burn_mapping/reference_data/Albers_Australia_Coast_Islands_Reefs.shp')

if label:
    index = albers[albers['label']==label].index[0]
    x = (albers.loc[index]['X_MIN'], albers.loc[index]['X_MAX'])
    y = (albers.loc[index]['Y_MIN'], albers.loc[index]['Y_MAX'])
    output_filename = outputdir + '/S2gm_2021-2022_'+'_'.join(label.split(','))+'.tif'
    print("Working on tile {}...".format(label))
else:
    x, y = (1385000.0, 1375000.0), (-4570000.0, -4580000.0)
    if subset:
        output_filename = 'S2_ARD_gm_2021_test_subset.tif'
    else:
        output_filename = 'S2_ARD_gm_2021_test_one.tif'
        

def load_ds(x, y):
    #This query from NCI_Geomedian_Tiles
    query = {'x': x,
          'y': y,
          'crs': 'EPSG:3577',
          'time': ('2021-10', '2022-04'),
          'measurements': ['nbart_blue', 'nbart_green', 'nbart_red', 'nbart_nir_1', 'nbart_swir_2'], # Can add nbart_swir_2 for true flase colour but change res to 20
          'resolution': (-10, 10),
          'group_by': 'solar_day',
          'output_crs': 'EPSG:28355'} #Albers is 'output_crs': 'EPSG:3577'
    
        
    # Load available data from both Sentinel 2 satellites
    ds = load_ard(dc=dc,
                  products=['s2a_ard_granule', 's2b_ard_granule'],
                  dtype='native',
                  dask_chunks={'time':1, "x": 1000, "y": 1000},
                  **query)

  
    # Compute geomedian here is necessary - either for dataset or subset months
    geomedian = int_geomedian(ds)
    myGeomed = geomedian.compute()
    return myGeomed
    
    # compute geomedian
    #ds_gm = GeoMedian().compute(ds)
    #return ds_gm.copy()       
        
#Compute and then concat subsets  
xm, ym = (x[0]+x[1])/2, (y[0]+y[1])/2
x1, x2 = (x[0], xm), (xm, x[1])
y1, y2 = (y[0], ym), (ym, y[1])
if subset:
    out1 = load_ds(x1, y)
    out2 = load_ds(x2, y)
    out = xr.concat([out1, out2], dim='x')
else:
    out = load_ds(x, y)

# Here we can export the geomedian
# for COG we need an array not a dataset
out_da = out.to_array()

# Write multi-band GeoTIFF to a location
print('Saving COG...')
write_cog(geo_im=out_da,
          fname=output_filename,
          overwrite=True)

print('Finished processing')
