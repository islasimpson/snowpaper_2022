import importlib
import xarray as xr
import numpy as np
import pandas as pd
import sys
import os

from CASutils import filter_utils as filt
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read

importlib.reload(filt)
importlib.reload(cal)

outpath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/trefhtptile_composites/pos_neg_extremes/"

def calcdeseas(da):
    datseas = da.groupby('time.dayofyear').mean('time', skipna=True)
    dat4harm = filt.calc_season_nharm(datseas, 4, dimtime=0)
    anoms = da.groupby('time.dayofyear') - dat4harm
    datdeseas = cal.group_season_daily(anoms, 'DJF')
    seasmean = datdeseas.mean('day', skipna=True)
    datdeseas = datdeseas - seasmean
    #datdeseas = np.array(datdeseas).flatten()
    return datdeseas

###read in CESM data to determine the temperature bins.

fpath="/project/cas02/islas/CLM5_CLM4/raw/Isla_CAM6_CLM5_002/day/TREFHT/"
trefht = read.read_sfc_cesm(fpath+"*.nc","1979-01-01","2014-12-31")
trefhtanoms = calcdeseas(trefht.TREFHT)
trefhtanoms_stacked = trefhtanoms.stack(time=("year","day"))

binlow_min = trefhtanoms_stacked.quantile(1./100., dim="time")
binlow_max = trefhtanoms_stacked.quantile(10./100., dim="time")

binhigh_min = trefhtanoms_stacked.quantile(90./100., dim="time")
binhigh_max = trefhtanoms_stacked.quantile(99./100., dim="time")

del(trefhtanoms)

###read in ERA5 data and composite

from math import nan


##T2m
fpath = "/project/mojave/observations/ERA5_daily/T2m/*.nc"
trefht = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
trefht = trefht.sel(time=~((trefht.time.dt.month==2) & (trefht.time.dt.day == 29)))
trefhtanoms = calcdeseas(trefht.t2m)
del(trefht)
trefhtanoms = trefhtanoms.stack(time=("year","day"))

mintrefht = trefhtanoms.where( (trefhtanoms >= binlow_min) & (trefhtanoms < binlow_max), nan).mean('time', skipna=True)
maxtrefht = trefhtanoms.where( (trefhtanoms >= binhigh_min) & (trefhtanoms < binhigh_max), nan).mean('time', skipna=True)
##End T2m

##SHFLX
fpath = "/project/mojave/observations/ERA5_daily/SHFLX/*.nc"
shflx = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
shflx = shflx.sel(time=~((shflx.time.dt.month==2) & (shflx.time.dt.day == 29)))
shflxanoms = calcdeseas(shflx.shflx)
del(shflx)
shflxanoms = shflxanoms.stack(time=("year","day"))

minshflx = shflxanoms.where( (trefhtanoms >= binlow_min) & (trefhtanoms < binlow_max), nan).mean('time', skipna=True)
maxshflx = shflxanoms.where( (trefhtanoms >= binhigh_min) & (trefhtanoms < binhigh_max), nan).mean('time', skipna=True)
del(shflxanoms)
##End SHFLX 

##LHFLX
fpath = "/project/mojave/observations/ERA5_daily/LHFLX/*.nc"
lhflx = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
lhflx = lhflx.sel(time=~((lhflx.time.dt.month==2) & (lhflx.time.dt.day == 29)))
lhflxanoms = calcdeseas(lhflx.lhflx)
del(lhflx)
lhflxanoms = lhflxanoms.stack(time=("year","day"))

minlhflx = lhflxanoms.where( (trefhtanoms >= binlow_min) & (trefhtanoms < binlow_max), nan).mean('time', skipna=True)
maxlhflx = lhflxanoms.where( (trefhtanoms >= binhigh_min) & (trefhtanoms < binhigh_max), nan).mean('time', skipna=True)
del(lhflxanoms)
##End LHFLX 


##FLNS
fpath = "/project/mojave/observations/ERA5_daily/STR/*.nc"
flns = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
flns = flns.sel(time=~((flns.time.dt.month==2) & (flns.time.dt.day == 29)))
flnsanoms = calcdeseas(flns.STR)
del(flns)
flnsanoms = flnsanoms.stack(time=("year","day"))

minflns = flnsanoms.where( (trefhtanoms >= binlow_min) & (trefhtanoms < binlow_max), nan).mean('time', skipna=True)
maxflns = flnsanoms.where( (trefhtanoms >= binhigh_min) & (trefhtanoms < binhigh_max), nan).mean('time', skipna=True)
del(flnsanoms)
##End FLNS


##FSNS
fpath = "/project/mojave/observations/ERA5_daily/SSR/*.nc"
fsns = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
fsns = fsns.sel(time=~((fsns.time.dt.month==2) & (fsns.time.dt.day == 29)))
fsnsanoms = calcdeseas(fsns.SSR)
del(fsns)
fsnsanoms = fsnsanoms.stack(time=("year","day"))

minfsns = fsnsanoms.where( (trefhtanoms >= binlow_min) & (trefhtanoms < binlow_max), nan).mean('time', skipna=True)
maxfsns = fsnsanoms.where( (trefhtanoms >= binhigh_min) & (trefhtanoms < binhigh_max), nan).mean('time', skipna=True)
del(fsnsanoms)
##End FSNS

##T850
fpath = "/project/haggis/ERA5/day/T850/t850*.nc"
t850 = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
t850 = t850.sel(time=~((t850.time.dt.month==2) & (t850.time.dt.day == 29)))
t850anoms = calcdeseas(t850.t850)
del(t850)
t850anoms = t850anoms.stack(time=("year","day"))

mint850 = t850anoms.where( (trefhtanoms >= binlow_min) & (trefhtanoms < binlow_max), nan).mean('time', skipna=True)
maxt850 = t850anoms.where( (trefhtanoms >= binhigh_min) & (trefhtanoms < binhigh_max), nan).mean('time', skipna=True)
del(t850anoms)
##End FSNS

##increments
fpath="/project/haggis/ERA5/forincrements/increments/t2m/*.nc"
increment = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
increment = increment.increment
increment = increment.sel(time=~((increment.time.dt.month==2) & (increment.time.dt.day == 29)))
incrementanoms = calcdeseas(increment)
del(increment)
incrementanoms = incrementanoms.stack(time=("year","day"))

minincrement = incrementanoms.where( (trefhtanoms >= binlow_min) & (trefhtanoms < binlow_max), nan).mean('time', skipna=True)
maxincrement = incrementanoms.where( (trefhtanoms >= binhigh_min) & (trefhtanoms < binhigh_max), nan).mean('time', skipna=True)
del(incrementanoms)
##End increment 

mintrefht = mintrefht.rename('mintrefht')
maxtrefht = maxtrefht.rename('maxtrefht')
minshflx = minshflx.rename('minshflx')
maxshflx = maxshflx.rename('maxshflx')
minlhflx = minlhflx.rename('minlhflx')
maxlhflx = maxlhflx.rename('maxlhflx')
minflns = minflns.rename('minflns')
maxflns = maxflns.rename('maxflns')
minfsns = minfsns.rename('minfsns')
maxfsns = maxfsns.rename('maxfsns')
mint850 = mint850.rename('mint850')
maxt850 = maxt850.rename('maxt850')
minincrement = minincrement.rename('minincrement')
maxincrement = maxincrement.rename('maxincrement')

mintrefht.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc")
maxtrefht.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
minshflx.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
maxshflx.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
minlhflx.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
maxlhflx.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
minflns.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
maxflns.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
minfsns.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
maxfsns.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
mint850.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
maxt850.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
minincrement.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")
maxincrement.to_netcdf(path=outpath+"ERA5_minmax_trefhtptilecomposite.nc", mode="a")










































