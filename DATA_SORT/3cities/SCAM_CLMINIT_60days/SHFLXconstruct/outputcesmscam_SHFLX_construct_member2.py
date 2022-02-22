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

expname=("SCAM_CLM5_CLM5F_02","SCAM_SNOWD_CLM5F_02")

basedir="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM_CLMINIT_60days/"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM_CLMINIT_60days/SHFLXconstruct/"

shflx = xr.open_dataset(basedir+"SHFLX_SCAM_CLM5_CLM5F_02.nc")
shflx_deseas = calcdeseas(shflx.shflx)
ts = xr.open_dataset(basedir+"TS_SCAM_CLM5_CLM5F_02.nc")
ts_deseas = calcdeseas(ts.ts)
tbot = xr.open_dataset(basedir+"TBOT_SCAM_CLM5_CLM5F_02.nc")
tbot_deseas = calcdeseas(tbot.tbot)

cityname=['Saskatoon','Toronto','Siderovsk']
citylon=[253.330, 280.617, 82.3139]
citylat=[52.1579, 43.6532, 66.5973]

a = np.zeros([len(cityname)]) ; b = np.zeros([len(cityname)])
for icity in np.arange(0,len(cityname),1):
    shflx = np.array(shflx_deseas.isel(city=icity)).flatten()
    ts = np.array(ts_deseas.isel(city=icity)).flatten()
    tbot = np.array(tbot_deseas.isel(city=icity)).flatten()
    a[icity], b[icity] = linfit.linfit_xy(tbot -ts, shflx)


for iexp in np.arange(0,len(expname),1):
    print(expname[iexp])

    fpath=basedir+"SHFLX_"+expname[iexp]+".nc"
    dat = xr.open_dataset(fpath)
    shflx_deseas = calcdeseas(dat.shflx)
    
    fpath=basedir+"TS_"+expname[iexp]+".nc"
    dat = xr.open_dataset(fpath)
    ts_deseas = calcdeseas(dat.ts)

    fpath=basedir+"TBOT_"+expname[iexp]+".nc"
    dat = xr.open_dataset(fpath)
    tbot_deseas = calcdeseas(dat.tbot)

    nyears = shflx_deseas.year.size
    ndays = shflx_deseas.day.size

    shflxconstruct = np.zeros([nyears*ndays,3])
    for icity in np.arange(0,len(cityname),1):
        shflx = np.array(shflx_deseas.isel(city=icity)).flatten()
        tbot = np.array(tbot_deseas.isel(city=icity)).flatten()
        ts = np.array(ts_deseas.isel(city=icity)).flatten()
        shflxconstruct[:,icity] = a[icity] + b[icity]*(tbot[:] - ts[:])

    shflxconstruct = np.reshape(shflxconstruct,[nyears,ndays,3])

    shflxconstruct_xr = xr.DataArray(shflxconstruct, 
      coords=[np.arange(1979,1979+nyears,1),np.arange(0,ndays,1),cityname],
      dims=["year","day","city"],name='shflxconstruct')
    shflxconstruct_xr.to_netcdf(path=pathout+"SHFLXconstruct_"+expname[iexp]+".nc")
