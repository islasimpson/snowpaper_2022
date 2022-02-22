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


#-----generate percentiles for compositing
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
shflx = xr.open_dataset(basepath+"ERA5_SHFLX.nc")
flns = xr.open_dataset(basepath+"ERA5_FLNS.nc")

cities = trefht.city
ncities = trefht.city.size

trefht_deseas = calcdeseas(trefht.era5)
shflx_deseas = calcdeseas(shflx.shflx)
flns_deseas = calcdeseas(flns.flns)

for icity in range(0,ncities,1):
    print(icity)
    trefht = np.array(trefht_deseas[:,:,icity]).flatten()
    shflx = np.array(shflx_deseas[:,:,icity]).flatten()
    flns = np.array(flns_deseas[:,:,icity]).flatten()

    if (icity == 0):
        trefhtcomp = np.zeros([nblocks, ncities])
        shflxcomp = np.zeros([nblocks, ncities])
        flnscomp = np.zeros([nblocks, ncities])

    for iblock in np.arange(0,nblocks,1):
        trefhtcomp[iblock,icity] = \
         (trefht[(trefht >= binmin[iblock]) & (trefht < binmax[iblock])]).mean()

        shflxcomp[iblock,icity] = \
         (shflx[(trefht >= binmin[iblock]) & (trefht < binmax[iblock])]).mean()

        flnscomp[iblock,icity] = \
         (flns[(trefht >= binmin[iblock]) & (trefht < binmax[iblock])]).mean()

trefht_xr = xr.DataArray(trefhtcomp,
    coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='trefht')
shflx_xr = xr.DataArray(shflxcomp,
    coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='shflx')
flns_xr = xr.DataArray(flnscomp,
    coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='flns')

try:
    os.remove(outpath+"trefhtptilecomposites_3cities_ERA5.nc")
except:
    pass

trefht_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_ERA5.nc")
shflx_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_ERA5.nc", mode="a")
flns_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_ERA5.nc", mode="a")
