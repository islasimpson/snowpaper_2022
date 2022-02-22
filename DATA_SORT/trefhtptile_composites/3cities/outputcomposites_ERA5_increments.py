import importlib
import xarray as xr
import numpy as np
import pandas as pd
import sys
import os

from CASutils import filter_utils as filt
from CASutils import calendar_utils as cal

importlib.reload(filt)
importlib.reload(cal)

def calcdeseas(da):
    datseas = da.groupby('time.dayofyear').mean('time', skipna=True)
    dat4harm = filt.calc_season_nharm(datseas, 4, dimtime=0)
    anoms = da.groupby('time.dayofyear') - dat4harm
    datdeseas = cal.group_season_daily(anoms, 'DJF')
    seasmean = datdeseas.mean('day', skipna=True)
    datdeseas = datdeseas - seasmean
    #datdeseas = np.array(datdeseas).flatten()
    return datdeseas

basepath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/CAM/"
trefht_clm5 = xr.open_dataset(basepath+"TREFHT_Isla_CAM6_CLM5_002.nc")
trefht_clm5_deseas = calcdeseas(trefht_clm5.trefht)

cities = trefht_clm5.city
ncities = trefht_clm5.city.size

for icity in range(0,ncities,1):
    trefht_clm5 = np.array(trefht_clm5_deseas[:,:,icity]).flatten()
    # calculate the ptile bin ranges
    nblocks = 10
    binmin  = np.empty([nblocks]) ; binmax = np.empty([nblocks])
    for iblock in np.arange(0,nblocks,1):
        binmin[iblock] = np.percentile(trefht_clm5,iblock*10)
        binmax[iblock] = np.percentile(trefht_clm5,iblock*10+10)
        if (iblock == 0):
            binmin[iblock] = np.percentile(trefht_clm5,1)
        if (iblock == (nblocks-1)):
            binmax[iblock] = np.percentile(trefht_clm5,99)

outpath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/trefhtptile_composites/3cities/"


basepath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/OBS/"
trefht = xr.open_dataset(basepath+"ERA5_TREFHT.nc")

basepath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/ERA5/"
dat = xr.open_dataset(basepath+"ERA5_increments.nc")
increments_deseas = calcdeseas(dat.increments)
forecast_deseas = calcdeseas(dat.forecast)
analysis_deseas = calcdeseas(dat.analysis)
trefht_deseas = calcdeseas(trefht.era5)

cities=dat.city
ncities = dat.city.size

for icity in range(0,ncities,1):
    trefht = np.array(trefht_deseas[:,:,icity]).flatten()
    increments = np.array(increments_deseas[:,:,icity]).flatten()
    forecast = np.array(forecast_deseas[:,:,icity]).flatten()
    analysis = np.array(analysis_deseas[:,:,icity]).flatten()

    if (icity == 0):
        incrementcomp = np.zeros([nblocks, ncities])
        forecastcomp = np.zeros([nblocks, ncities])
        analysiscomp = np.zeros([nblocks, ncities])

    for iblock in np.arange(0,nblocks,1):
        incrementcomp[iblock, icity] = \
          (increments[(analysis >= binmin[iblock]) & (analysis < binmax[iblock])]).mean()

        forecastcomp[iblock, icity] = \
          (forecast[(analysis >= binmin[iblock]) & (analysis < binmax[iblock])]).mean()

        analysiscomp[iblock, icity] = \
          (analysis[(analysis >= binmin[iblock]) & (analysis < binmax[iblock])]).mean()

increment_xr = xr.DataArray(incrementcomp,
        coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='increment')
forecast_xr = xr.DataArray(forecastcomp,
        coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='forecast')
analysis_xr = xr.DataArray(analysiscomp,
        coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='analysis')

increment_xr.to_netcdf(path=outpath+'trefhtptilecomposites_3cities_ERA5increments.nc')
forecast_xr.to_netcdf(path=outpath+'trefhtptilecomposites_3cities_ERA5increments.nc', mode='a')
analysis_xr.to_netcdf(path=outpath+'trefhtptilecomposites_3cities_ERA5increments.nc', mode='a')


 

