import importlib
import pandas as pd
import xarray as xr
import numpy as np
import xesmf as xe
import sys
import warnings
from scipy import signal
import datetime

from ecpaper_utils import readdata_utils as read
from ecpaper_utils import calendar_utils as cal
from ecpaper_utils import shapefile_utils as shp
from ecpaper_utils import averaging_utils as avg
from ecpaper_utils import linfit_utils as linfit

cesmdat = xr.open_dataset(
'/project/cas02/islas/CLM5_CLM4/raw/Isla_CAM6_CLM5/day/TREFHT/TREFHT_1979.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat)}, {'lon': (['lon'], cesmdat.lon)})

pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/era5_ps/"
wgtfile=pathout+"wgtfile.nc"

ps = read.read_sfc("/project/haggis/ERA5/mon/PS/*.nc","1979-01","2014-12")

regridder = xe.Regridder(ps.ps, grid_out, 'bilinear', periodic=True, reuse_weights=False,
  filename=wgtfile)

ps_rg = regridder(ps.ps)

psseason = ps_rg.groupby('time.season').mean('time')

psdjf = psseason.sel(season='DJF') 

psdjf.to_netcdf(pathout+"PS_ERA5_1979_2014_DJF.nc")
