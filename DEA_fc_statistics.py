import pandas as pd
import geopandas as gpd
import xarray as xr
import sys, os
import time
import multiprocessing
ncpus = multiprocessing.cpu_count()
#from BurnCube import BurnCube #including burn mapping main functions
#bc = BurnCube()

from datacube_stats.statistics import Percentile

import datacube
from datacube.storage import masking
from datacube.storage.storage import write_dataset_to_netcdf
from datacube.helpers import write_geotiff

sys.path.append(os.path.abspath('/g/data/r78/DPIPWE_lm/dea-notebooks/10_Scripts'))
import DEAPlotting
import DEADataHandling

outputdir = '/g/data/r78/DPIPWE_lm/test_burn_mapping/output_data'

dc = datacube.Datacube(app='dc-FC')

################################################################################
query_PL = {
                'time': ('1990', '2018'),
                'lat': (-40.9, -41.7),
                'long': (146.5, 147.5),
                'resolution': (-25,25)
                }

ds_tam = DEADataHandling.load_clearlandsat(dc=dc, query=query_PL, product='fc', ls7_slc_off=True, masked_prop=0.3)

################################################################################

ds_tam = ds_tam.drop(['UE', 'data_perc'])
percentiles = [0, 5, 20, 50, 80, 95, 100]
FC_percents = Percentile(percentiles)
FC_percentiles = FC_percents.compute(ds_tam)
FC_percentiles.attrs = ds_tam.attrs

FC_percentiles.attrs['units'] = 'fractional_cover_percentage_percentile'

################################################################################

try:
    ds = FC_percentiles
    write_geotiff(filename='/g/data/r78/DPIPWE_lm/output_data/FC_percentiles_tamar.tif', dataset=ds)
    print('wrote to GeoTiff' )
#         DEADataHandling.write_your_netcdf(FC_quantiles.isel(quantile=quant), 'FC_Q'+str(quantiles[quant]), savefilepath+'FC_Q_'+(str(quantiles[quant]).replace('.','_'))+'.nc', crs = ds.crs)
except RuntimeError as err:
    print("RuntimeError: {0}".format(err))
