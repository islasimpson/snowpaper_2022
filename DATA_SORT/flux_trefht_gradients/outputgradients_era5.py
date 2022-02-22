import importlib
import xarray as xr
import numpy as np
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
trefht = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
trefht = trefht.sel(time=~((trefht.time.dt.month==2) & (trefht.time.dt.day == 29)))
trefhtanoms = deseasonalized(trefht.t2m)
del(trefht)

print("SHFLX")
fpath="/project/mojave/observations/ERA5_daily/SHFLX/*.nc"
shflx = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
shflx = shflx.sel(time=~((shflx.time.dt.month==2) & (shflx.time.dt.day == 29)))
shflxanoms = deseasonalized(shflx.shflx)
del(shflx)

print("LHFLX")
fpath="/project/mojave/observations/ERA5_daily/LHFLX/*.nc"
lhflx = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
lhflx = lhflx.sel(time=~((lhflx.time.dt.month==2) & (lhflx.time.dt.day == 29)))
lhflxanoms = deseasonalized(lhflx.lhflx)
del(lhflx)

print("FLNS")
fpath="/project/mojave/observations/ERA5_daily/STR/*.nc"
flns = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
flns = flns.sel(time=~((flns.time.dt.month==2) & (flns.time.dt.day == 29)))
flnsanoms = deseasonalized(flns.STR)
del(flns)

print("FSNS")
fpath="/project/mojave/observations/ERA5_daily/SSR/*.nc"
fsns = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
fsns = fsns.sel(time=~((fsns.time.dt.month==2) & (fsns.time.dt.day == 29)))
fsnsanoms = deseasonalized(fsns.SSR)
del(fsns)

print("T850")
fpath="/project/haggis/ERA5/day/T850/t850*.nc"
t850 = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
t850 = t850.sel(time=~((t850.time.dt.month==2) & (t850.time.dt.day == 29)))
t850anoms = deseasonalized(t850.t850)
del(t850)

print("increments")
fpath = "/project/haggis/ERA5/forincrements/increments/t2m/*.nc"
dat = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
increment = dat.increment
increment = increment.sel(time=~((increment.time.dt.month==2) & (increment.time.dt.day == 29)))
incrementanoms = deseasonalized(increment)
del increment

sumflux = -1.*fsnsanoms - flnsanoms - shflxanoms - lhflxanoms
netrad = -1.*fsnsanoms - flnsanoms

trefhtstacked = trefhtanoms.stack(time=('year','day'))
shflxstacked = shflxanoms.stack(time=('year','day'))
lhflxstacked = lhflxanoms.stack(time=('year','day'))
flnsstacked = flnsanoms.stack(time=('year','day'))
fsnsstacked = fsnsanoms.stack(time=('year','day'))
netradstacked = netrad.stack(time=('year','day'))
sumfluxstacked = sumflux.stack(time=('year','day'))
t850stacked = t850anoms.stack(time=('year','day'))
incrementstacked = incrementanoms.stack(time=('year','day'))


bshflx, ashflx, r, p, stderr = xr.apply_ufunc(
    linregress, trefhtstacked, shflxstacked, input_core_dims=[['time'],['time']],
    output_core_dims=[[],[],[],[],[]], vectorize=True)

blhflx, alhflx, r, p, stderr = xr.apply_ufunc(
    linregress, trefhtstacked, lhflxstacked, input_core_dims=[['time'],['time']],
    output_core_dims=[[],[],[],[],[]], vectorize=True)

bflns, aflns, r, p, stderr = xr.apply_ufunc(
    linregress, trefhtstacked, flnsstacked, input_core_dims=[['time'],['time']],
    output_core_dims=[[],[],[],[],[]], vectorize=True)

bfsns, afsns, r, p, stderr = xr.apply_ufunc(
    linregress, trefhtstacked, fsnsstacked, input_core_dims=[['time'],['time']],
    output_core_dims=[[],[],[],[],[]], vectorize=True)

bt850, at850, r, p, stderr = xr.apply_ufunc(
    linregress, trefhtstacked, t850stacked, input_core_dims=[['time'],['time']],
    output_core_dims=[[],[],[],[],[]], vectorize=True)

bsumflux, asumflux, r, p, stderr = xr.apply_ufunc(
    linregress, trefhtstacked, sumfluxstacked, input_core_dims=[['time'],['time']],
    output_core_dims=[[],[],[],[],[]], vectorize=True)

bnetrad, anetrad, r, p, stderr = xr.apply_ufunc(
    linregress, trefhtstacked, netradstacked, input_core_dims=[['time'],['time']],
    output_core_dims=[[],[],[],[],[]], vectorize=True)

bincrement, aincrement, r, p, stderr = xr.apply_ufunc(
    linregress, trefhtstacked, incrementstacked, input_core_dims=[['time'],['time']],
    output_core_dims=[[],[],[],[],[]], vectorize=True)

bshflx = bshflx.rename("bshflx")
blhflx = blhflx.rename("blhflx")
bflns = bflns.rename("bflns")
bfsns = bfsns.rename("bfsns")
bt850 = bt850.rename("bt850")
bsumflux = bsumflux.rename("bsumflux")
bnetrad = bnetrad.rename('bnetrad')
bincrement = bincrement.rename("bincrement")

bshflx.to_netcdf(path=pathout+"gradients_ERA5.nc")
blhflx.to_netcdf(path=pathout+"gradients_ERA5.nc", mode="a")
bflns.to_netcdf(path=pathout+"gradients_ERA5.nc",mode="a")
bfsns.to_netcdf(path=pathout+"gradients_ERA5.nc",mode="a")
bt850.to_netcdf(path=pathout+"gradients_ERA5.nc", mode="a")
bnetrad.to_netcdf(path=pathout+"gradients_ERA5.nc", mode="a")
bsumflux.to_netcdf(path=pathout+"gradients_ERA5.nc", mode="a")
bincrement.to_netcdf(path=pathout+"gradients_ERA5.nc", mode="a")


#trefhtanoms.to_netcdf(path=pathout+"anomalydata_ERA5.nc")
#shflxanoms.to_netcdf(path=pathout+"anomalydata_ERA5.nc", mode="a")
#flnsanoms.to_netcdf(path=pathout+"anomalydata_ERA5.nc", mode="a")
#
#shflx.to_netcdf(path=pathout+"actualshflx_ERA5.nc")



