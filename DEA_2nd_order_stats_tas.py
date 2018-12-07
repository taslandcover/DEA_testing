import warnings; warnings.simplefilter('ignore')

import fnmatch
import os
import pandas as pd
import geopandas as gpd

from datacube_stats.statistics import GeoMedian
from datacube.helpers import ga_pq_fuser
from datacube.storage import masking
from datacube.helpers import write_geotiff
import xarray as xr

#get the DEA version of the plotting functions
import sys
sys.path.append(os.path.abspath('/g/data/r78/DPIPWE_lm/dea-notebooks/10_Scripts'))
sys.path.append(os.path.abspath('/g/data/r78/DPIPWE_lm/datacube-2nd-order-stats'))
import DEAPlotting
import DEADataHandling
from model import SMAD, BCMAD, EMAD, TernaryMAD


import datacube
dc = datacube.Datacube(app='stats_2nd_testing')

###############################################################################

#x, y = (1385000.0, 1390000.0), (-4570000.0, -4575000.0)
sensors = ['ls8', 'ls7'] #take or remove as needed
deriv = 'nbart'
#product = 'nbart' #
time = ('2016-01-01', '2016-12-31')
resolution = (-25,25)
bands = ['red', 'green', 'blue', 'nir', 'swir1', 'swir2']
#epoch = ('2016', '2017') # time query for datacube function can be just years

query = {'x': (1300000.0, 1400000.0),
         'y': (-4700000.0, -4800000.0),
         'time': time,
         'resolution': resolution,
         'crs': 'EPSG:3577'}

###############################################################################

dsma = DEADataHandling.load_clearlandsat(dc=dc, query=query,
                                          #product=product,
                                          masked_prop=0,
                                          sensors = sensors,
                                          bands_of_interest = bands,
                                          mask_pixel_quality=True,
                                          ls7_slc_off=True)

dsma = dsma.drop('data_perc')

##############################################################################

dsma_smad = TernaryMAD().compute(dsma)

ds=xr.Dataset({'smad': (['y','x'], dsma_smad.sdev),
               'emad': (['y','x'], dsma_smad.edev),
               'bcmad': (['y','x'], dsma_smad.bcdev)},
                coords={'x': dsma.x, 'y':dsma.y}, attrs=dsma.attrs)

#datacube.storage.storage.write_dataset_to_netcdf(dsma_smad, '/g/data/r78/DPIPWE_LM/output_data/ls8_smad_test.nc')
datacube.helpers.write_geotiff(filename='/g/data/r78/DPIPWE_lm/output_data/ls8_TMAD_lscomb2016.tif', dataset=ds)
#DEADataHandling.dataset_to_geotiff('dsma_smad_netcdf_test.nc', dsma_smad)
