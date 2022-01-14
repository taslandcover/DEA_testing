'''
 Adapted from the burn mapping script but to produce the S2
 geomedian image using ARD. Takes a list of tiles in TILELIST
'''

import warnings

import datacube
import numpy as np
from odc.algo import int_geomedian, to_f32, xr_geomedian

warnings.filterwarnings("ignore")

import sys, os

sys.path.append(os.path.abspath('/g/data/r78/DPIPWE_lm/dea-notebooks/Tools'))
#from dea_tools.dask import create_local_dask_cluster
from dea_tools.datahandling import load_ard
#from dea_tools.plotting import rgb

outputdir = '/g/data/r78/DPIPWE_lm/test_geomedian_mapping/output_data'
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
    output_filename = outputdir + '/gm_2016-2017_'+'_'.join(label.split(','))+'.nc'
    print("Working on tile {}...".format(label))
else:
    x, y = (1385000.0, 1375000.0), (-4570000.0, -4580000.0)
    if subset:
        output_filename = 'gm_2019-2020_test_subset.nc'
    else:
        output_filename = 'gm_2019-2020_test_one.nc'
        
if os.path.exists(output_filename):
    print("output file already exists.")
    exit() 

""
# Edit these variables as needed

#sensor = 'ls8'
time = ('2019-07', '2020-06') # period of interest for change/severity mapping
resolution = (10, 10)

""
def gm_comp(x, y):
    # Create a reusable query (can change resolution and/or bands here if needed)
    query = {
        "x": x,
        "y": y,
        "time": ("2019-07", "2020-06"),  # Months or day can be specified here
        "measurements": [
            "nbart_blue",
            "nbart_green",
            "nbart_red",
            "nbart_nir_1",
            "nbart_swir_2",
        ],  # Can add nbart_swir_2 for true flase colour but change res to 20 (or resample 10m)
        "output_crs": "EPSG:3577",
        "resolution": (-20, 20),
        "group_by": "solar_day",
    }
    
    # Load available data from both Sentinel 2 satellites
    ds = load_ard(dc=dc,
                  products=['s2a_ard_granule', 's2b_ard_granule'],
                  min_gooddata=0.1, # At least 10% good data
                  **query)
    
    # compute geomedian
    ds_gm = GeoMedian().compute(ds)
    return ds_gm.copy()

""
xm, ym = (x[0]+x[1])/2, (y[0]+y[1])/2
x1, x2 = (x[0], xm), (xm, x[1])
y1, y2 = (y[0], ym), (ym, y[1])
if subset:
    out1 = gm_comp(x1, y)
    out2 = gm_comp(x2, y)
    out = xr.concat([out1, out2], dim='x')
else:
    out = gm_comp(x, y)


# Export the geomedian
# for COG we need an array not a dataset
da = ds_cm_gm.to_array()

# Write multi-band GeoTIFF to a location
write_cog(geo_im=out,
          fname=output_filename,
          overwrite=True)
