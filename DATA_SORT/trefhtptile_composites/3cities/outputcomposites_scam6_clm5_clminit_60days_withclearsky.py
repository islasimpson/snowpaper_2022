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

basepath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM_CLMINIT_60days_withclearsky/"
outpath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/trefhtptile_composites/3cities/"

trefht_clm5_1 = xr.open_dataset(basepath+"TREFHT_SCAM_CLM5_CLM5F_01.nc")
trefht_snowd_1 = xr.open_dataset(basepath+"TREFHT_SCAM_SNOWD_CLM5F_01.nc")
flns_clm5_1 = xr.open_dataset(basepath+"FLNS_SCAM_CLM5_CLM5F_01.nc")
flns_snowd_1 = xr.open_dataset(basepath+"FLNS_SCAM_SNOWD_CLM5F_01.nc")
flnsc_clm5_1 = xr.open_dataset(basepath+"FLNSC_SCAM_CLM5_CLM5F_01.nc")
flnsc_snowd_1 = xr.open_dataset(basepath+"FLNSC_SCAM_SNOWD_CLM5F_01.nc")


cities = trefht_clm5_1.city
ncities = trefht_clm5_1.city.size


trefht_clm5_deseas = calcdeseas(trefht_clm5_1.trefht)
trefht_snowd_deseas = calcdeseas(trefht_snowd_1.trefht)

flns_clm5_deseas = calcdeseas(flns_clm5_1.flns)
flns_snowd_deseas = calcdeseas(flns_snowd_1.flns)

flnsc_clm5_deseas = calcdeseas(flnsc_clm5_1.flnsc)
flnsc_snowd_deseas = calcdeseas(flnsc_snowd_1.flnsc)


for icity in range(0,ncities,1):
    print(icity)
    trefht_clm5 = np.array(trefht_clm5_deseas[:,:,icity]).flatten()
    trefht_snowd = np.array(trefht_snowd_deseas[:,:,icity]).flatten()

    flns_clm5 = np.array(flns_clm5_deseas[:,:,icity]).flatten()
    flns_snowd = np.array(flns_snowd_deseas[:,:,icity]).flatten()

    flnsc_clm5 = np.array(flnsc_clm5_deseas[:,:,icity]).flatten()
    flnsc_snowd = np.array(flnsc_snowd_deseas[:,:,icity]).flatten()


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

    # composite
    if (icity == 0):
        trefhtcomp_clm5 = np.zeros([nblocks, ncities])
        trefhtcomp_snowd = np.zeros([nblocks, ncities])

        flnscomp_clm5 = np.zeros([nblocks, ncities])
        flnscomp_snowd = np.zeros([nblocks, ncities])

        flnsccomp_clm5 = np.zeros([nblocks, ncities])
        flnsccomp_snowd = np.zeros([nblocks, ncities])




    for iblock in np.arange(0,nblocks,1):
        trefhtcomp_clm5[iblock, icity] = \
        (trefht_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        trefhtcomp_snowd[iblock, icity] = \
        (trefht_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        flnscomp_clm5[iblock, icity] = \
        (flns_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        flnscomp_snowd[iblock, icity] = \
        (flns_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        flnsccomp_clm5[iblock, icity] = \
        (flnsc_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        flnsccomp_snowd[iblock, icity] = \
        (flnsc_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()




trefht_clm5_xr = xr.DataArray(trefhtcomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='trefht_clm5')
trefht_snowd_xr = xr.DataArray(trefhtcomp_snowd,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='trefht_snowd')

flns_clm5_xr = xr.DataArray(flnscomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='flns_clm5')
flns_snowd_xr = xr.DataArray(flnscomp_snowd,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='flns_snowd')

flnsc_clm5_xr = xr.DataArray(flnsccomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='flnsc_clm5')
flnsc_snowd_xr = xr.DataArray(flnsccomp_snowd,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='flnsc_snowd')



try:
    os.remove(outpath+"trefhtptilecomposites_3cities_scam_clminit_60days_clearsky.nc")
except:
    pass


trefht_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days_clearsky.nc")
trefht_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days_clearsky.nc", 
                  mode="a")
flns_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days_clearsky.nc", 
                  mode="a")
flns_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days_clearsky.nc", 
                  mode="a")
flnsc_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days_clearsky.nc", 
                  mode="a")
flnsc_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days_clearsky.nc", 
                  mode="a")
















