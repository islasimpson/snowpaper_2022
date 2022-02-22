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



basepath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/offlineland/"
outpath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/trefhtptile_composites/3cities/"

trefht_clm5 = xr.open_dataset(basepath+"TSA_offlineland_clm5.nc")
shflx_clm5 = xr.open_dataset(basepath+"FSH_offlineland_clm5.nc")
flns_clm5 = xr.open_dataset(basepath+"FLNS_offlineland_clm5.nc")
fgr_clm5 = xr.open_dataset(basepath+"FGR_offlineland_clm5.nc")

trefht_snowd = xr.open_dataset(basepath+"TSA_offlineland_snowrevert.nc")
shflx_snowd = xr.open_dataset(basepath+"FSH_offlineland_snowrevert.nc")
flns_snowd = xr.open_dataset(basepath+"FLNS_offlineland_snowrevert.nc")
fgr_snowd = xr.open_dataset(basepath+"FGR_offlineland_snowrevert.nc")

cities = trefht_clm5.city
ncities = trefht_clm5.city.size

trefht_clm5_deseas = calcdeseas(trefht_clm5.trefht)
shflx_clm5_deseas = calcdeseas(shflx_clm5.shflx)
flns_clm5_deseas = calcdeseas(flns_clm5.flns)
fgr_clm5_deseas = calcdeseas(fgr_clm5.fgr)

trefht_snowd_deseas = calcdeseas(trefht_snowd.trefht)
shflx_snowd_deseas = calcdeseas(shflx_snowd.shflx)
flns_snowd_deseas = calcdeseas(flns_snowd.flns)
fgr_snowd_deseas = calcdeseas(fgr_snowd.fgr)

for icity in range(0,ncities,1):
    print(icity)
    trefht_clm5 = np.array(trefht_clm5_deseas[:,:,icity]).flatten()
    trefht_snowd = np.array(trefht_snowd_deseas[:,:,icity]).flatten()
    shflx_clm5 = np.array(shflx_clm5_deseas[:,:,icity]).flatten()
    shflx_snowd = np.array(shflx_snowd_deseas[:,:,icity]).flatten()
    flns_clm5 = np.array(flns_clm5_deseas[:,:,icity]).flatten()
    flns_snowd = np.array(flns_snowd_deseas[:,:,icity]).flatten()
    fgr_clm5 = np.array(fgr_clm5_deseas[:,:,icity]).flatten()
    fgr_snowd = np.array(fgr_snowd_deseas[:,:,icity]).flatten()

#    # calculate the ptile bin ranges
#    nblocks = 10
#    binmin  = np.empty([nblocks]) ; binmax = np.empty([nblocks])
#    for iblock in np.arange(0,nblocks,1):
#        binmin[iblock] = np.percentile(trefht_clm5,iblock*10)
#        binmax[iblock] = np.percentile(trefht_clm5,iblock*10+10)
#        if (iblock == 0):
#            binmin[iblock] = np.percentile(trefht_clm5,1)
#        if (iblock == (nblocks-1)):
#            binmax[iblock] = np.percentile(trefht_clm5,99)


    # composite
    if (icity == 0):
        trefhtcomp_clm5 = np.zeros([nblocks, ncities])
        trefhtcomp_snowd = np.zeros([nblocks, ncities])
        shflxcomp_clm5 = np.zeros([nblocks, ncities])
        shflxcomp_snowd = np.zeros([nblocks, ncities])
        flnscomp_clm5 = np.zeros([nblocks, ncities])
        flnscomp_snowd = np.zeros([nblocks, ncities])
        fgrcomp_clm5 = np.zeros([nblocks, ncities])
        fgrcomp_snowd = np.zeros([nblocks, ncities])

    for iblock in np.arange(0,nblocks,1):
        trefhtcomp_clm5[iblock, icity] = \
         (trefht_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        trefhtcomp_snowd[iblock, icity] = \
         (trefht_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        shflxcomp_clm5[iblock, icity] = \
         (shflx_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        shflxcomp_snowd[iblock, icity] = \
         (shflx_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        flnscomp_clm5[iblock, icity] = \
         (flns_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        flnscomp_snowd[iblock, icity] = \
         (flns_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        fgrcomp_clm5[iblock, icity] = \
         (fgr_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        fgrcomp_snowd[iblock, icity] = \
         (fgr_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()


trefht_clm5_xr = xr.DataArray(trefhtcomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='trefht_clm5')
shflx_clm5_xr = xr.DataArray(shflxcomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='shflx_clm5')
flns_clm5_xr = xr.DataArray(flnscomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='flns_clm5')
fgr_clm5_xr = xr.DataArray(fgrcomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='fgr_clm5')

trefht_snowd_xr = xr.DataArray(trefhtcomp_snowd,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='trefht_snowd')
shflx_snowd_xr = xr.DataArray(shflxcomp_snowd,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='shflx_snowd')
flns_snowd_xr = xr.DataArray(flnscomp_snowd,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='flns_snowd')
fgr_snowd_xr = xr.DataArray(fgrcomp_snowd,
                    coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='fgr_snowd')

try:
    os.remove( outpath+"trefhtptilecomposites_3cities_offlineland.nc")
except:
    pass

trefht_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_offlineland.nc")
shflx_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_offlineland.nc", mode="a")
flns_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_offlineland.nc", mode="a")
fgr_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_offlineland.nc", mode="a")

trefht_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_offlineland.nc", mode="a")
shflx_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_offlineland.nc", mode="a")
flns_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_offlineland.nc", mode="a")
fgr_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_offlineland.nc", mode="a")
























