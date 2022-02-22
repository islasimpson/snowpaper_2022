import importlib
import xarray as xr
import numpy as np
import pandas as pd

from CASutils import filter_utils as filt
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal
from CASutils import linfit_utils as linfit

importlib.reload(filt)
importlib.reload(read)
importlib.reload(cal)

import sys

def calcdeseas(da):
    datseas = da.groupby('time.dayofyear').mean('time', skipna=True)
    dat4harm = filt.calc_season_nharm(datseas, 4, dimtime=0)
    anoms = da.groupby('time.dayofyear') - dat4harm
    datdeseas = cal.group_season_daily(anoms, 'DJF')
    seasmean = datdeseas.mean('day', skipna=True)
    datdeseas = datdeseas - seasmean
    #datdeseas = np.array(datdeseas).flatten()
    return datdeseas

cityname=['Saskatoon','Toronto','Siderovsk']
citylon=[253.330, 280.617, 82.3139]
citylat=[52.1579, 43.6532, 66.5973]

basedir="/project/cas02/islas/CLM5_CLM4/raw/"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/CAM/"
 
# calculate regression coefficients for SHFLX = a + b(TBOT-TS) for CAM6_CLM5

fpath=basedir+"Isla_CAM6_CLM5_002/day/SHFLX/*.nc"
dat = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
shflx = xr.DataArray(np.zeros([dat.time.size,3]), coords=[dat.time, cityname],
                    dims=['time','city'], name='shflx')
for icity in np.arange(0,len(cityname),1):
    shflx[:,icity] = dat.SHFLX.sel(lon=citylon[icity], lat=citylat[icity], method='nearest')

fpath=basedir+"Isla_CAM6_CLM5_002/day/TS/*.nc"
dat = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
ts = xr.DataArray(np.zeros([dat.time.size,3]), coords=[dat.time, cityname],
                    dims=['time','city'], name='ts')
for icity in np.arange(0,len(cityname),1):
    ts[:,icity] = dat.TS.sel(lon=citylon[icity], lat=citylat[icity], method='nearest')

fpath=basedir+"Isla_CAM6_CLM5_002/day/TBOT/*.nc"
dat = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
tbot = xr.DataArray(np.zeros([dat.time.size,3]), coords=[dat.time, cityname],
                    dims=['time','city'], name='tbot')
for icity in np.arange(0,len(cityname),1):
    tbot[:,icity] = dat.TBOT.sel(lon=citylon[icity], lat=citylat[icity], method='nearest')

shflx_deseas = calcdeseas(shflx)
ts_deseas = calcdeseas(ts)
tbot_deseas = calcdeseas(tbot)

a=np.zeros([len(cityname)]) ; b = np.zeros([len(cityname)])
for icity in np.arange(0,len(cityname),1):
    shflx = np.array(shflx_deseas.isel(city=icity)).flatten()
    ts = np.array(ts_deseas.isel(city=icity)).flatten()
    tbot = np.array(tbot_deseas.isel(city=icity)).flatten()
    a[icity], b[icity] = linfit.linfit_xy( tbot - ts, shflx)

expname=('Isla_CAM6_CLM5_002','CAM6_CLM5_snowdensity_002')

sys.exit()

for iexp in np.arange(0,len(expname),1):

    print(expname[iexp])
    
    fpath=basedir+expname[iexp]+"/day/SHFLX/*.nc"
    dat = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
    shflx = xr.DataArray(np.zeros([dat.time.size,3]), coords=[dat.time, cityname],
                        dims=['time','city'], name='shflx')
    for icity in np.arange(0,len(cityname),1):
        shflx[:,icity] = dat.SHFLX.sel(lon=citylon[icity], lat=citylat[icity], method='nearest')

    fpath=basedir+expname[iexp]+"/day/TS/*.nc"
    dat = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
    ts = xr.DataArray(np.zeros([dat.time.size,3]), coords=[dat.time, cityname],
                        dims=['time','city'], name='ts')
    for icity in np.arange(0,len(cityname),1):
        ts[:,icity] = dat.TS.sel(lon=citylon[icity], lat=citylat[icity], method='nearest')

    fpath=basedir+expname[iexp]+"/day/TBOT/*.nc"
    dat = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
    tbot = xr.DataArray(np.zeros([dat.time.size,3]), coords=[dat.time, cityname],
                        dims=['time','city'], name='tbot')
    for icity in np.arange(0,len(cityname),1):
        tbot[:,icity] = dat.TBOT.sel(lon=citylon[icity], lat=citylat[icity], method='nearest')

    shflx_deseas = calcdeseas(shflx)
    ts_deseas = calcdeseas(ts)
    tbot_deseas = calcdeseas(tbot)

    shflxconstruct = np.zeros([dat.time.size,3])
    for icity in np.arange(0,len(cityname),1):
        shflxconstruct[:,icity] = a[icity] + b[icity]*(tbot[:,icity] - ts[:,icity])

    shflxconstruct_xr = xr.DataArray(shflxconstruct, coords=shflx.coords, name='shflxconstruct')

    
    shflxconstruct_xr.to_netcdf(path=pathout+"SHFLXconstruct_"+expname[iexp]+".nc")
