import xarray as xr
import sys
from math import nan
import numpy as np
import pandas as pd

# thermal conductivities (W/m)
lamice=2.29 ; lamair=0.023

expname=['SCAM_CLM5_CLM5F_01','SCAM_CLM5_CLM5F_02','SCAM_SNOWD_SNOWDF_01','SCAM_SNOWD_SNOWDF_02']

datpath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM/new_lowrelax/"
outpath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM/new_lowrelax/BULKSNOW/"

for iexp in expname:
    print(iexp)


