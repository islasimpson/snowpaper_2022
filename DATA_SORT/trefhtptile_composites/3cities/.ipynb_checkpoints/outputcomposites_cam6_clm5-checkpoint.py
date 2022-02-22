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

expname='Isla_CAM6_CLM5'

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
outpath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/trefhtptile_composites/3cities/"

trefht_clm5 = xr.open_dataset(basepath+"TREFHT_Isla_CAM6_CLM5_002.nc")
tbot_clm5 = xr.open_dataset(basepath+"TBOT_Isla_CAM6_CLM5_002.nc")
ts_clm5 = xr.open_dataset(basepath+"TS_Isla_CAM6_CLM5_002.nc")
shflx_clm5 = xr.open_dataset(basepath+"SHFLX_Isla_CAM6_CLM5_002.nc")
shflxconstruct_clm5 = xr.open_dataset(basepath+"SHFLXconstruct_Isla_CAM6_CLM5_002.nc")
lhflx_clm5 = xr.open_dataset(basepath+"LHFLX_Isla_CAM6_CLM5_002.nc")
flns_clm5 = xr.open_dataset(basepath+"FLNS_Isla_CAM6_CLM5_002.nc")
flds_clm5 = xr.open_dataset(basepath+"FLDS_Isla_CAM6_CLM5_002.nc")
fsns_clm5 = xr.open_dataset(basepath+"FSNS_Isla_CAM6_CLM5_002.nc")
fgr_clm5 = xr.open_dataset(basepath+"FGR_Isla_CAM6_CLM5_002.nc")
bulksnow_clm5 = xr.open_dataset(basepath+"BULKSNOW/BULKSNOW_Isla_CAM6_CLM5_002.nc")
bulksnow_condfix_clm5 = xr.open_dataset(basepath+"BULKSNOW/BULKSNOW_clm5conductance_Isla_CAM6_CLM5_002.nc")

trefht_snowd = xr.open_dataset(basepath+"TREFHT_CAM6_CLM5_snowdensity_002.nc")
tbot_snowd = xr.open_dataset(basepath+"TBOT_CAM6_CLM5_snowdensity_002.nc")
ts_snowd = xr.open_dataset(basepath+"TS_CAM6_CLM5_snowdensity_002.nc")
shflx_snowd = xr.open_dataset(basepath+"SHFLX_CAM6_CLM5_snowdensity_002.nc")
shflxconstruct_snowd = xr.open_dataset(basepath+"SHFLXconstruct_CAM6_CLM5_snowdensity_002.nc")
lhflx_snowd = xr.open_dataset(basepath+"LHFLX_CAM6_CLM5_snowdensity_002.nc")
flds_snowd = xr.open_dataset(basepath+"FLDS_CAM6_CLM5_snowdensity_002.nc")
flns_snowd = xr.open_dataset(basepath+"FLNS_CAM6_CLM5_snowdensity_002.nc")
fsns_snowd = xr.open_dataset(basepath+"FSNS_CAM6_CLM5_snowdensity_002.nc")
fgr_snowd = xr.open_dataset(basepath+"FGR_CAM6_CLM5_snowdensity_002.nc")
bulksnow_snowd = xr.open_dataset(basepath+"BULKSNOW/BULKSNOW_CAM6_CLM5_snowdensity_002.nc")
bulksnow_condfix_snowd = xr.open_dataset(basepath+"BULKSNOW/BULKSNOW_clm5conductance_CAM6_CLM5_snowdensity_002.nc")

cities = trefht_clm5.city
ncities = trefht_clm5.city.size

trefht_clm5_deseas = calcdeseas(trefht_clm5.trefht)
tbot_clm5_deseas = calcdeseas(tbot_clm5.tbot)
ts_clm5_deseas = calcdeseas(ts_clm5.ts)
shflx_clm5_deseas = calcdeseas(shflx_clm5.shflx)
shflxconstruct_clm5_deseas = calcdeseas(shflxconstruct_clm5.shflxconstruct)
lhflx_clm5_deseas = calcdeseas(lhflx_clm5.lhflx)
flns_clm5_deseas = calcdeseas(flns_clm5.flns)
flds_clm5_deseas = calcdeseas(flds_clm5.flds)
fsns_clm5_deseas = calcdeseas(fsns_clm5.fsns)
fgr_clm5_deseas = calcdeseas(fgr_clm5.fgr)
bulksnowflux_clm5_deseas = calcdeseas(bulksnow_clm5.snowflux)
bulksnowflux_condfix_clm5_deseas = calcdeseas(bulksnow_condfix_clm5.snowflux)

trefht_snowd_deseas = calcdeseas(trefht_snowd.trefht)
ts_snowd_deseas = calcdeseas(ts_snowd.ts)
tbot_snowd_deseas = calcdeseas(tbot_snowd.tbot)
shflx_snowd_deseas = calcdeseas(shflx_snowd.shflx)
shflxconstruct_snowd_deseas = calcdeseas(shflxconstruct_snowd.shflxconstruct)
lhflx_snowd_deseas = calcdeseas(lhflx_snowd.lhflx)
flns_snowd_deseas = calcdeseas(flns_snowd.flns)
flds_snowd_deseas = calcdeseas(flds_snowd.flds)
fsns_snowd_deseas = calcdeseas(fsns_snowd.fsns)
fgr_snowd_deseas = calcdeseas(fgr_snowd.fgr)
bulksnowflux_snowd_deseas = calcdeseas(bulksnow_snowd.snowflux)
bulksnowflux_condfix_snowd_deseas = calcdeseas(bulksnow_condfix_snowd.snowflux)



for icity in range(0,ncities,1):
    print(icity)
    trefht_clm5 = np.array(trefht_clm5_deseas[:,:,icity]).flatten()
    trefht_snowd = np.array(trefht_snowd_deseas[:,:,icity]).flatten()
    ts_clm5 = np.array(ts_clm5_deseas[:,:,icity]).flatten()
    ts_snowd = np.array(ts_snowd_deseas[:,:,icity]).flatten()
    tbot_clm5 = np.array(tbot_clm5_deseas[:,:,icity]).flatten()
    tbot_snowd = np.array(tbot_snowd_deseas[:,:,icity]).flatten()
    shflx_clm5 = np.array(shflx_clm5_deseas[:,:,icity]).flatten()
    shflx_snowd = np.array(shflx_snowd_deseas[:,:,icity]).flatten()
    shflxconstruct_clm5 = np.array(shflxconstruct_clm5_deseas[:,:,icity]).flatten()
    shflxconstruct_snowd = np.array(shflxconstruct_snowd_deseas[:,:,icity]).flatten()
    lhflx_clm5 = np.array(lhflx_clm5_deseas[:,:,icity]).flatten()
    lhflx_snowd = np.array(lhflx_snowd_deseas[:,:,icity]).flatten()
    flns_clm5 = np.array(flns_clm5_deseas[:,:,icity]).flatten()
    flns_snowd = np.array(flns_snowd_deseas[:,:,icity]).flatten()
    flds_clm5 = np.array(flds_clm5_deseas[:,:,icity]).flatten()
    flds_snowd = np.array(flds_snowd_deseas[:,:,icity]).flatten()
    fsns_clm5 = np.array(fsns_clm5_deseas[:,:,icity]).flatten()
    fsns_snowd = np.array(fsns_snowd_deseas[:,:,icity]).flatten()
    fgr_clm5 = np.array(fgr_clm5_deseas[:,:,icity]).flatten()
    fgr_snowd = np.array(fgr_snowd_deseas[:,:,icity]).flatten()
    bulksnowflux_clm5 = np.array(bulksnowflux_clm5_deseas[:,:,icity]).flatten()
    bulksnowflux_snowd = np.array(bulksnowflux_snowd_deseas[:,:,icity]).flatten()
    bulksnowflux_condfix_clm5 = np.array(bulksnowflux_condfix_clm5_deseas[:,:,icity]).flatten()
    bulksnowflux_condfix_snowd = np.array(bulksnowflux_condfix_snowd_deseas[:,:,icity]).flatten()

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
        tscomp_clm5 = np.zeros([nblocks, ncities])
        tscomp_snowd = np.zeros([nblocks, ncities])
        tbotcomp_clm5 = np.zeros([nblocks, ncities])
        tbotcomp_snowd = np.zeros([nblocks, ncities])
        shflxcomp_clm5 = np.zeros([nblocks, ncities])
        shflxcomp_snowd = np.zeros([nblocks, ncities])
        shflxconstructcomp_clm5 = np.zeros([nblocks, ncities])
        shflxconstructcomp_snowd = np.zeros([nblocks, ncities])
        lhflxcomp_clm5 = np.zeros([nblocks, ncities])
        lhflxcomp_snowd = np.zeros([nblocks, ncities])
        flnscomp_clm5 = np.zeros([nblocks, ncities])
        flnscomp_snowd = np.zeros([nblocks, ncities])
        fldscomp_clm5 = np.zeros([nblocks, ncities])
        fldscomp_snowd = np.zeros([nblocks, ncities])
        fsnscomp_clm5 = np.zeros([nblocks, ncities])
        fsnscomp_snowd = np.zeros([nblocks, ncities])
        fgrcomp_clm5 = np.zeros([nblocks, ncities])
        fgrcomp_snowd = np.zeros([nblocks, ncities])
        snowfluxcomp_clm5 = np.zeros([nblocks, ncities])
        snowfluxcomp_snowd = np.zeros([nblocks, ncities])
        snowfluxcomp_condfix_clm5 = np.zeros([nblocks, ncities])
        snowfluxcomp_condfix_snowd = np.zeros([nblocks, ncities])





    for iblock in np.arange(0,nblocks,1):
        trefhtcomp_clm5[iblock, icity] = \
         (trefht_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        trefhtcomp_snowd[iblock, icity] = \
         (trefht_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        tbotcomp_clm5[iblock, icity] = \
         (tbot_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        tbotcomp_snowd[iblock, icity] = \
         (tbot_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        tscomp_clm5[iblock, icity] = \
         (ts_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        tscomp_snowd[iblock, icity] = \
         (ts_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        shflxcomp_clm5[iblock, icity] = \
         (shflx_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        shflxcomp_snowd[iblock, icity] = \
         (shflx_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        shflxconstructcomp_clm5[iblock, icity] = \
         (shflxconstruct_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        shflxconstructcomp_snowd[iblock, icity] = \
         (shflxconstruct_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        lhflxcomp_clm5[iblock, icity] = \
         (lhflx_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        lhflxcomp_snowd[iblock, icity] = \
         (lhflx_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        flnscomp_clm5[iblock, icity] = \
         (flns_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        flnscomp_snowd[iblock, icity] = \
         (flns_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        fldscomp_clm5[iblock, icity] = \
         (flds_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        fldscomp_snowd[iblock, icity] = \
         (flds_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        fsnscomp_clm5[iblock, icity] = \
         (fsns_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        fsnscomp_snowd[iblock, icity] = \
         (fsns_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        fgrcomp_clm5[iblock, icity] = \
         (fgr_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        fgrcomp_snowd[iblock, icity] = \
         (fgr_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        snowfluxcomp_clm5[iblock, icity] = \
         np.nanmean(bulksnowflux_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])])

        snowfluxcomp_snowd[iblock, icity] = \
         np.nanmean(bulksnowflux_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])])

        snowfluxcomp_condfix_clm5[iblock, icity] = \
         np.nanmean(bulksnowflux_condfix_clm5[ 
           (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])])

        snowfluxcomp_condfix_snowd[iblock, icity] = \
         np.nanmean(bulksnowflux_condfix_snowd[ 
            (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])])




trefht_clm5_xr = xr.DataArray(trefhtcomp_clm5, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='trefht_clm5')
ts_clm5_xr = xr.DataArray(tscomp_clm5, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='ts_clm5')
tbot_clm5_xr = xr.DataArray(tbotcomp_clm5, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='tbot_clm5')
shflx_clm5_xr = xr.DataArray(shflxcomp_clm5, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='shflx_clm5')
shflxconstruct_clm5_xr = xr.DataArray(shflxconstructcomp_clm5, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='shflxconstruct_clm5')
lhflx_clm5_xr = xr.DataArray(lhflxcomp_clm5, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='lhflx_clm5')
flns_clm5_xr = xr.DataArray(flnscomp_clm5, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='flns_clm5')
flds_clm5_xr = xr.DataArray(fldscomp_clm5, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='flds_clm5')
fsns_clm5_xr = xr.DataArray(fsnscomp_clm5, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='fsns_clm5')
fgr_clm5_xr = xr.DataArray(fgrcomp_clm5, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='fgr_clm5')
snowflux_clm5_xr = xr.DataArray(snowfluxcomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='snowflux_clm5')
snowflux_condfix_clm5_xr = xr.DataArray(snowfluxcomp_condfix_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='snowflux_condfix_clm5')



trefht_snowd_xr = xr.DataArray(trefhtcomp_snowd, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='trefht_snowd')
ts_snowd_xr = xr.DataArray(tscomp_snowd, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='ts_snowd')
tbot_snowd_xr = xr.DataArray(tbotcomp_snowd, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='tbot_snowd')
shflx_snowd_xr = xr.DataArray(shflxcomp_snowd, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='shflx_snowd')
shflxconstruct_snowd_xr = xr.DataArray(shflxconstructcomp_snowd, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='shflxconstruct_snowd')
lhflx_snowd_xr = xr.DataArray(lhflxcomp_snowd, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='lhflx_snowd')
flns_snowd_xr = xr.DataArray(flnscomp_snowd, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='flns_snowd')
flds_snowd_xr = xr.DataArray(fldscomp_snowd, 
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='flds_snowd')
fsns_snowd_xr = xr.DataArray(fsnscomp_snowd,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='fsns_snowd')
fgr_snowd_xr = xr.DataArray(fgrcomp_snowd, 
                    coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='fgr_snowd')
snowflux_snowd_xr = xr.DataArray(snowfluxcomp_snowd,
                    coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='snowflux_snowd')
snowflux_condfix_snowd_xr = xr.DataArray(snowfluxcomp_condfix_snowd,
                    coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'], name='snowflux_condfix_snowd')


#!rm outpath"trefhtptilecomposites_3cities.nc"
os.remove( outpath+"trefhtptilecomposites_3cities.nc")


trefht_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc")
ts_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
tbot_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
shflx_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
shflxconstruct_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
lhflx_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
flns_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
flds_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
fsns_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
fgr_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
snowflux_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
snowflux_condfix_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")

trefht_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
ts_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
tbot_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
shflx_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
shflxconstruct_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
lhflx_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
flns_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
flds_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
fsns_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
fgr_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
snowflux_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")
snowflux_condfix_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities.nc", mode="a")




