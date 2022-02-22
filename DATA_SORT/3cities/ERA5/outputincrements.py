import importlib
import xarray as xr
import numpy as np
import pandas as pd

from CASutils import readdata_utils as read

path="/project/haggis/ERA5/forincrements/increments/t2m/"

dat = read.read_sfc(path+"*.nc","1979-01-01","2014-12-31")
dat = dat.sel(time=~((dat.time.dt.month==2) & (dat.time.dt.day == 29)))

pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/ERA5/"

cityname=['Saskatoon','Toronto','Siderovsk']
citylon=[253.330, 280.617, 82.3139]
citylat=[52.1579, 43.6532, 66.5973]

era5increments = xr.DataArray(np.zeros([dat.time.size,3]), coords=[dat.time, cityname], 
                    dims=['time','city'], name='increments')
era5forecast = xr.DataArray(np.zeros([dat.time.size,3]), coords=[dat.time, cityname],
                    dims=['time','city'], name='forecast')
era5analysis = xr.DataArray(np.zeros([dat.time.size,3]), coords=[dat.time, cityname],
                    dims=['time','city'], name='analysis')

for icity in np.arange(0,len(cityname),1):
    era5increments[:,icity] = dat.increment.sel(lon=citylon[icity], lat=citylat[icity], method='nearest')
    era5forecast[:,icity] = dat.forecast.sel(lon=citylon[icity], lat=citylat[icity], method='nearest')
    era5analysis[:,icity] = dat.analysis.sel(lon=citylon[icity], lat=citylat[icity], method='nearest')

era5increments.to_netcdf(path=pathout+"ERA5_increments.nc")
era5forecast.to_netcdf(path=pathout+"ERA5_increments.nc", mode="a")
era5analysis.to_netcdf(path=pathout+"ERA5_increments.nc", mode="a")
