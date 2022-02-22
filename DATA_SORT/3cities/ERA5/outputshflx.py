import importlib
import xarray as xr
import numpy as np
import pandas as pd

from CASutils import readdata_utils as read

path="/project/mojave/observations/ERA5_daily/SHFLX/"

dat = read.read_sfc(path+"*.nc","1979-01-01","2014-12-31")
dat = dat.sel(time=~((dat.time.dt.month==2) & (dat.time.dt.day == 29)))

pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/ERA5/"

cityname=['Saskatoon','Toronto','Siderovsk']
citylon=[253.330, 280.617, 82.3139]
citylat=[52.1579, 43.6532, 66.5973]

era5shflx = xr.DataArray(np.zeros([dat.time.size,3]), coords=[dat.time, cityname], 
                    dims=['time','city'], name='shflx')

for icity in np.arange(0,len(cityname),1):
    era5shflx[:,icity] = dat.shflx.sel(lon=citylon[icity], lat=citylat[icity], method='nearest')

era5shflx.to_netcdf(path=pathout+"ERA5_SHFLX.nc")
