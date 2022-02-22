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

basepath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM_CLMINIT_60days/"
outpath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/trefhtptile_composites/3cities/"

trefht_clm5_1 = xr.open_dataset(basepath+"TREFHT_SCAM_CLM5_CLM5F_01.nc")
trefht_clm5_2 = xr.open_dataset(basepath+"TREFHT_SCAM_CLM5_CLM5F_02.nc")
trefht_snowd_1 = xr.open_dataset(basepath+"TREFHT_SCAM_SNOWD_SNOWDF_01.nc")
trefht_snowd_2 = xr.open_dataset(basepath+"TREFHT_SCAM_SNOWD_SNOWDF_02.nc")
tbot_clm5_1 = xr.open_dataset(basepath+"TBOT_SCAM_CLM5_CLM5F_01.nc")
tbot_clm5_2 = xr.open_dataset(basepath+"TBOT_SCAM_CLM5_CLM5F_02.nc")
tbot_snowd_1 = xr.open_dataset(basepath+"TBOT_SCAM_SNOWD_SNOWDF_01.nc")
tbot_snowd_2 = xr.open_dataset(basepath+"TBOT_SCAM_SNOWD_SNOWDF_02.nc")
ts_clm5_1 = xr.open_dataset(basepath+"TS_SCAM_CLM5_CLM5F_01.nc")
ts_clm5_2 = xr.open_dataset(basepath+"TS_SCAM_CLM5_CLM5F_02.nc")
ts_snowd_1 = xr.open_dataset(basepath+"TS_SCAM_SNOWD_SNOWDF_01.nc")
ts_snowd_2 = xr.open_dataset(basepath+"TS_SCAM_SNOWD_SNOWDF_02.nc")
tsl_clm5_1 = xr.open_dataset(basepath+"TSL_SCAM_CLM5_CLM5F_01.nc")
tsl_clm5_2 = xr.open_dataset(basepath+"TSL_SCAM_CLM5_CLM5F_02.nc")
tsl_snowd_1 = xr.open_dataset(basepath+"TSL_SCAM_SNOWD_SNOWDF_01.nc")
tsl_snowd_2 = xr.open_dataset(basepath+"TSL_SCAM_SNOWD_SNOWDF_02.nc")
snot_clm5_1 = xr.open_dataset(basepath+"SNO_T_SCAM_CLM5_CLM5F_01.nc")
snot_clm5_2 = xr.open_dataset(basepath+"SNO_T_SCAM_CLM5_CLM5F_02.nc")
snot_snowd_1 = xr.open_dataset(basepath+"SNO_T_SCAM_SNOWD_SNOWDF_01.nc")
snot_snowd_2 = xr.open_dataset(basepath+"SNO_T_SCAM_SNOWD_SNOWDF_02.nc")
t_clm5_1 = xr.open_dataset(basepath+"T_SCAM_CLM5_CLM5F_01.nc")
t_clm5_2 = xr.open_dataset(basepath+"T_SCAM_CLM5_CLM5F_02.nc")
t_snowd_1 = xr.open_dataset(basepath+"T_SCAM_SNOWD_SNOWDF_01.nc")
t_snowd_2 = xr.open_dataset(basepath+"T_SCAM_SNOWD_SNOWDF_02.nc")





cities = trefht_clm5_1.city
ncities = trefht_clm5_1.city.size
levsno = snot_clm5_1.levsno
nlevsno = snot_clm5_1.levsno.size
lev = t_clm5_1.lev
nlev = t_clm5_1.lev.size


trefht_clm5_1_deseas = calcdeseas(trefht_clm5_1.trefht)
trefht_clm5_2_deseas = calcdeseas(trefht_clm5_2.trefht)
trefht_clm5_deseas = xr.concat([trefht_clm5_1_deseas,trefht_clm5_2_deseas], dim='year')

trefht_snowd_1_deseas = calcdeseas(trefht_snowd_1.trefht)
trefht_snowd_2_deseas = calcdeseas(trefht_snowd_2.trefht)
trefht_snowd_deseas = xr.concat([trefht_snowd_1_deseas,trefht_snowd_2_deseas], dim='year')

tbot_clm5_1_deseas = calcdeseas(tbot_clm5_1.tbot)
tbot_clm5_2_deseas = calcdeseas(tbot_clm5_2.tbot)
tbot_clm5_deseas = xr.concat([tbot_clm5_1_deseas,tbot_clm5_2_deseas], dim='year')

tbot_snowd_1_deseas = calcdeseas(tbot_snowd_1.tbot)
tbot_snowd_2_deseas = calcdeseas(tbot_snowd_2.tbot)
tbot_snowd_deseas = xr.concat([tbot_snowd_1_deseas,tbot_snowd_2_deseas], dim='year')

ts_clm5_1_deseas = calcdeseas(ts_clm5_1.ts)
ts_clm5_2_deseas = calcdeseas(ts_clm5_2.ts)
ts_clm5_deseas = xr.concat([ts_clm5_1_deseas,ts_clm5_2_deseas], dim='year')

ts_snowd_1_deseas = calcdeseas(ts_snowd_1.ts)
ts_snowd_2_deseas = calcdeseas(ts_snowd_2.ts)
ts_snowd_deseas = xr.concat([ts_snowd_1_deseas,ts_snowd_2_deseas], dim='year')

tsl_clm5_1_deseas = calcdeseas(tsl_clm5_1.tsl)
tsl_clm5_2_deseas = calcdeseas(tsl_clm5_2.tsl)
tsl_clm5_deseas = xr.concat([tsl_clm5_1_deseas,tsl_clm5_2_deseas], dim='year')

tsl_snowd_1_deseas = calcdeseas(tsl_snowd_1.tsl)
tsl_snowd_2_deseas = calcdeseas(tsl_snowd_2.tsl)
tsl_snowd_deseas = xr.concat([tsl_snowd_1_deseas,tsl_snowd_2_deseas], dim='year')

snot_clm5_1_deseas = calcdeseas(snot_clm5_1.snot)
snot_clm5_2_deseas = calcdeseas(snot_clm5_2.snot)
snot_clm5_deseas = xr.concat([snot_clm5_1_deseas,snot_clm5_2_deseas], dim='year')

snot_snowd_1_deseas = calcdeseas(snot_snowd_1.snot)
snot_snowd_2_deseas = calcdeseas(snot_snowd_2.snot)
snot_snowd_deseas = xr.concat([snot_snowd_1_deseas,snot_snowd_2_deseas], dim='year')

t_clm5_1_deseas = calcdeseas(t_clm5_1.T)
t_clm5_2_deseas = calcdeseas(t_clm5_2.T)
t_clm5_deseas = xr.concat([t_clm5_1_deseas,t_clm5_2_deseas], dim='year')

t_snowd_1_deseas = calcdeseas(t_snowd_1.T)
t_snowd_2_deseas = calcdeseas(t_snowd_2.T)
t_snowd_deseas = xr.concat([t_snowd_1_deseas,t_snowd_2_deseas], dim='year')

for icity in range(0,ncities,1):
    print(icity)
    trefht_clm5 = np.array(trefht_clm5_deseas[:,:,icity]).flatten()
    trefht_snowd = np.array(trefht_snowd_deseas[:,:,icity]).flatten()

    ts_clm5 = np.array(ts_clm5_deseas[:,:,icity]).flatten()
    ts_snowd = np.array(ts_snowd_deseas[:,:,icity]).flatten()

    tsl_clm5 = np.array(tsl_clm5_deseas[:,:,icity]).flatten()
    tsl_snowd = np.array(tsl_snowd_deseas[:,:,icity]).flatten()

    tbot_clm5 = np.array(tbot_clm5_deseas[:,:,icity]).flatten()
    tbot_snowd = np.array(tbot_snowd_deseas[:,:,icity]).flatten()

    snot_clm5 = np.array(snot_clm5_deseas[:,:,:,icity])
    snot_snowd = np.array(snot_snowd_deseas[:,:,:,icity])

    t_clm5 = np.array(t_clm5_deseas[:,:,:,icity])
    t_snowd = np.array(t_snowd_deseas[:,:,:,icity])

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
        tslcomp_clm5 = np.zeros([nblocks, ncities])
        tslcomp_snowd = np.zeros([nblocks, ncities])
        tbotcomp_clm5 = np.zeros([nblocks, ncities])
        tbotcomp_snowd = np.zeros([nblocks, ncities])

        snotcomp_clm5 = np.zeros([nblocks, nlevsno, ncities])
        snotcomp_snowd = np.zeros([nblocks, nlevsno, ncities])

        tcomp_clm5 = np.zeros([nblocks, nlev, ncities])
        tcomp_snowd = np.zeros([nblocks, nlev, ncities])


    for iblock in np.arange(0,nblocks,1):
        trefhtcomp_clm5[iblock, icity] = \
        (trefht_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        trefhtcomp_snowd[iblock, icity] = \
        (trefht_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        tscomp_clm5[iblock, icity] = \
        (ts_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        tscomp_snowd[iblock, icity] = \
        (ts_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        tslcomp_clm5[iblock, icity] = \
        (tsl_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        tslcomp_snowd[iblock, icity] = \
        (tsl_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()

        tbotcomp_clm5[iblock, icity] = \
        (tbot_clm5[ (trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])]).mean()

        tbotcomp_snowd[iblock, icity] = \
        (tbot_snowd[ (trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])]).mean()


        for ilev in np.arange(0,nlevsno,1):
            datclm5 = snot_clm5[:,:,ilev].flatten()
            snotcomp_clm5[iblock, ilev, icity] = np.nanmean(
             datclm5[(trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])])

            datsnowd = snot_snowd[:,:,ilev].flatten()
            snotcomp_snowd[iblock, ilev, icity] = np.nanmean(
             datsnowd[(trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])])


        for ilev in np.arange(0,nlev,1):
            datclm5 = t_clm5[:,:,ilev].flatten()
            tcomp_clm5[iblock, ilev, icity] = np.nanmean(
             datclm5[(trefht_clm5 >= binmin[iblock]) & (trefht_clm5 < binmax[iblock])])

            datsnowd = t_snowd[:,:,ilev].flatten()
            tcomp_snowd[iblock, ilev, icity] = np.nanmean(
             datsnowd[(trefht_snowd >= binmin[iblock]) & (trefht_snowd < binmax[iblock])])




trefht_clm5_xr = xr.DataArray(trefhtcomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='trefht_clm5')
trefht_snowd_xr = xr.DataArray(trefhtcomp_snowd,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='trefht_snowd')

tbot_clm5_xr = xr.DataArray(tbotcomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='tbot_clm5')
tbot_snowd_xr = xr.DataArray(tbotcomp_snowd,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='tbot_snowd')

ts_clm5_xr = xr.DataArray(tscomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='ts_clm5')
ts_snowd_xr = xr.DataArray(tscomp_snowd,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='ts_snowd')

tsl_clm5_xr = xr.DataArray(tslcomp_clm5,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='tsl_clm5')
tsl_snowd_xr = xr.DataArray(tslcomp_snowd,
                 coords=[np.arange(0,nblocks,1),cities], dims=['ptile','city'],
                 name='tsl_snowd')

snot_clm5_xr = xr.DataArray(snotcomp_clm5,
               coords=[np.arange(0,nblocks,1),levsno,cities], dims=['ptile','levsno','city'],
               name='snot_clm5')
snot_snowd_xr = xr.DataArray(snotcomp_snowd,
               coords=[np.arange(0,nblocks,1),levsno,cities], dims=['ptile','levsno','city'],
               name='snot_snowd')

t_clm5_xr = xr.DataArray(tcomp_clm5,
               coords=[np.arange(0,nblocks,1),lev,cities], dims=['ptile','lev','city'],
               name='t_clm5')
t_snowd_xr = xr.DataArray(tcomp_snowd,
               coords=[np.arange(0,nblocks,1),lev,cities], dims=['ptile','lev','city'],
               name='t_snowd')





try:
    os.remove(outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc")
except:
    pass


trefht_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc")
trefht_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc", mode="a")
ts_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc", mode="a")
ts_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc", mode="a")
tsl_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc", mode="a")
tsl_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc", mode="a")
tbot_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc", mode="a")
tbot_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc", mode="a")
snot_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc", mode="a")
snot_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc", mode="a")
t_clm5_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc", mode="a")
t_snowd_xr.to_netcdf(path=outpath+"trefhtptilecomposites_3cities_scam_clminit_60days.nc", mode="a")















