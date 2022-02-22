import importlib
import xarray as xr
import numpy as np
import pandas as pd
import sys

from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal
from CASutils import filter_utils as filt

from scipy.stats import linregress

importlib.reload(read)
importlib.reload(cal)
importlib.reload(filt)

pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/flux_trefht_gradients/"

def deseasonalized(data):
    daystr = xr.DataArray(data.indexes['time'].strftime('%m-%d'), coords = data.time.coords, name="daystr")
    datseason = data.groupby(daystr).mean('time')
    dat4harm = filt.calc_season_nharm(datseason,4,dimtime=0)
    datanoms = data.groupby(daystr)-dat4harm
    datanoms = cal.group_season_daily(datanoms,'DJF')
    datmean = datanoms.mean('day')
    datanoms = datanoms - datmean
    return datanoms

print("T2m")
fpath = "/project/mojave/observations/ERA5_daily/T2m/*.nc"
trefht = read.read_sfc(fpath,"1979-01-01","2014-12-31")
trefht = trefht.sel(time=~((trefht.time.dt.month==2) & (trefht.time.dt.day == 29)))
trefhtanoms = deseasonalized(trefht.t2m)

print("SHFLX")
fpath = "/project/cas02/islas/FluxCom/SHFLX/SHFLX_*.nc"
shflx = read.read_sfc(fpath,"1979-01-01","2014-12-31")
shflx = shflx.sel(time=~((shflx.time.dt.month==2) & (shflx.time.dt.day == 29)))
shflx = shflx*(1000000./86400.)
shflxanoms = deseasonalized(shflx.shflx)

print("netrad")
fpath = "/project/cas02/islas/FluxCom/netrad/netrad_*.nc"
netrad = read.read_sfc(fpath,"1979-01-01","2014-12-31")
netrad = netrad.sel(time=~((netrad.time.dt.month==2) & (netrad.time.dt.day == 29)))
netrad = netrad*(1000000./86400.)
netradanoms = deseasonalized(netrad.netrad)

trefhtstacked = trefhtanoms.stack(time=('year','day'))
shflxstacked = shflxanoms.stack(time=('year','day'))
netradstacked = netradanoms.stack(time=('year','day'))

bshflx, ashflx, r, p, stderr = xr.apply_ufunc(
    linregress, trefhtstacked, shflxstacked, input_core_dims=[['time'],['time']],
    output_core_dims=[[],[],[],[],[]], vectorize=True)

bnetrad, anetrad, r, p, stderr = xr.apply_ufunc(
    linregress, trefhtstacked, netradstacked, input_core_dims=[['time'],['time']],
    output_core_dims=[[],[],[],[],[]], vectorize=True)

bshflx = bshflx.rename("bshflx")
bnetrad = bnetrad.rename("bnetrad")

bshflx.to_netcdf(path=pathout+"gradients_fluxcom.nc")
bnetrad.to_netcdf(path=pathout+"gradients_fluxcom.nc",mode="a")

trefhtanoms.to_netcdf(path=pathout+"anomalydata_fluxcom.nc")
shflxanoms.to_netcdf(path=pathout+"anomalydata_fluxcom.nc", mode="a")
netradanoms.to_netcdf(path=pathout+"anomalydata_fluxcom.nc", mode="a")

shflx.to_netcdf(path=pathout+"actualshflx_fluxcom.nc")


