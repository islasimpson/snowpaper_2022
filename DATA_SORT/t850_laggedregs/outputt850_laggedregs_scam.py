import importlib
import xarray as xr
import numpy as np
import pandas as pd
import sys
import os
from math import nan
from CASutils import calendar_utils as cal
from CASutils import filter_utils as filt

importlib.reload(cal)
importlib.reload(filt)

pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/t850_laggedregs/"

def deseasonalize(dat):
    datseas = dat.groupby('time.dayofyear').mean('time')
    dat4harm = filt.calc_season_nharm(datseas,4,dimtime=1)
    datanoms = dat.groupby('time.dayofyear') - dat4harm
    datdjfanoms = cal.group_season_daily(datanoms,'DJF')
    datmean = datdjfanoms.mean('day')
    datdjfanoms = datdjfanoms - datmean
    return datdjfanoms

def readanddeseas(file1, file2, var):
    dat1 = xr.open_dataset(file1).sel(time=slice("1979-01-01","2014-12-31"))
    dat1 = dat1[var]
    dat1deseas = deseasonalize(dat1)

    dat2 = xr.open_dataset(file2).sel(time=slice("1979-01-01","2014-12-31"))
    dat2 = dat2[var]
    dat2deseas = deseasonalize(dat2)

    nyears = dat1deseas.year.size + dat2deseas.year.size
    datout = np.zeros([nyears, dat1deseas.day.size, dat1deseas.city.size])
    datout[0:dat1deseas.year.size,:,:] = dat1deseas
    datout[dat1deseas.year.size:nyears,:,:] = dat2deseas
 
    return datout

#def readanddeseas(file1,var):
#    dat = xr.open_dataset(file1).sel(time=slice("1979-01-01","2014-12-31"))
#    dat = dat[var]
#    datdeseas = deseasonalize(dat)
#    return datdeseas

cities=['Saskatoon','Toronto','Siderovsk']

filepath = "/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM_CLMINIT_60days/"

t850_clm5 = readanddeseas(filepath+"T850_SCAM_CLM5_CLM5F_01.nc",filepath+"T850_SCAM_CLM5_CLM5F_02.nc",
                          "t850")

t850_snowd = readanddeseas(filepath+"T850_SCAM_SNOWD_CLM5F_01.nc",filepath+"T850_SCAM_SNOWD_CLM5F_02.nc",
                           "t850")

trefht_clm5 = readanddeseas(filepath+"TREFHT_SCAM_CLM5_CLM5F_01.nc",
                            filepath+"TREFHT_SCAM_CLM5_CLM5F_02.nc",
                            "trefht")

trefht_snowd = readanddeseas(filepath+"TREFHT_SCAM_SNOWD_CLM5F_01.nc",
                             filepath+"TREFHT_SCAM_SNOWD_CLM5F_02.nc",
                             "trefht")

ts_clm5 = readanddeseas(filepath+"TS_SCAM_CLM5_CLM5F_01.nc",
                        filepath+"TS_SCAM_CLM5_CLM5F_02.nc",
                        "ts")

ts_snowd = readanddeseas(filepath+"TS_SCAM_SNOWD_CLM5F_01.nc",
                         filepath+"TS_SCAM_SNOWD_CLM5F_02.nc",
                         "ts")

tbot_clm5 = readanddeseas(filepath+"TBOT_SCAM_CLM5_CLM5F_01.nc",
                          filepath+"TBOT_SCAM_CLM5_CLM5F_02.nc",
                          "tbot")

tbot_snowd = readanddeseas(filepath+"TBOT_SCAM_SNOWD_CLM5F_01.nc",
                           filepath+"TBOT_SCAM_SNOWD_CLM5F_02.nc",
                           "tbot")

flns_clm5 = readanddeseas(filepath+"FLNS_SCAM_CLM5_CLM5F_01.nc",
                          filepath+"FLNS_SCAM_CLM5_CLM5F_02.nc",
                          "flns")

flns_snowd = readanddeseas(filepath+"FLNS_SCAM_SNOWD_CLM5F_01.nc",
                           filepath+"FLNS_SCAM_SNOWD_CLM5F_02.nc",
                           "flns")

flds_clm5 = readanddeseas(filepath+"FLDS_SCAM_CLM5_CLM5F_01.nc",
                          filepath+"FLDS_SCAM_CLM5_CLM5F_02.nc",
                          "flds")

flds_snowd = readanddeseas(filepath+"FLDS_SCAM_SNOWD_CLM5F_01.nc",
                           filepath+"FLDS_SCAM_SNOWD_CLM5F_02.nc",
                           "flds")

fsns_clm5 = readanddeseas(filepath+"FSNS_SCAM_CLM5_CLM5F_01.nc",
                          filepath+"FSNS_SCAM_CLM5_CLM5F_02.nc",
                          "fsns")

fsns_snowd = readanddeseas(filepath+"FSNS_SCAM_SNOWD_CLM5F_01.nc",
                           filepath+"FSNS_SCAM_SNOWD_CLM5F_02.nc",
                           "fsns")

lhflx_clm5 = readanddeseas(filepath+"LHFLX_SCAM_CLM5_CLM5F_01.nc",
                           filepath+"LHFLX_SCAM_CLM5_CLM5F_02.nc",
                           "lhflx")

lhflx_snowd = readanddeseas(filepath+"LHFLX_SCAM_SNOWD_CLM5F_01.nc",
                            filepath+"LHFLX_SCAM_SNOWD_CLM5F_02.nc",
                            "lhflx")

shflx_clm5 = readanddeseas(filepath+"SHFLX_SCAM_CLM5_CLM5F_01.nc",
                           filepath+"SHFLX_SCAM_CLM5_CLM5F_02.nc",
                           "shflx")

shflx_snowd = readanddeseas(filepath+"SHFLX_SCAM_SNOWD_CLM5F_01.nc",
                            filepath+"SHFLX_SCAM_SNOWD_CLM5F_02.nc",
                            "shflx")

fgr_clm5 = readanddeseas(filepath+"FGR_SCAM_CLM5_CLM5F_01.nc",
                           filepath+"FGR_SCAM_CLM5_CLM5F_02.nc",
                           "fgr")

fgr_snowd = readanddeseas(filepath+"FGR_SCAM_SNOWD_CLM5F_01.nc",
                            filepath+"FGR_SCAM_SNOWD_CLM5F_02.nc",
                            "fgr")

bulksnow_clm5 = readanddeseas(filepath+"BULKSNOW/BULKSNOW_SCAM_CLM5_CLM5F_01.nc",
                              filepath+"BULKSNOW/BULKSNOW_SCAM_CLM5_CLM5F_02.nc",
                              'snowflux')

bulksnow_snowd = readanddeseas(filepath+"BULKSNOW/BULKSNOW_SCAM_SNOWD_CLM5F_01.nc",
                              filepath+"BULKSNOW/BULKSNOW_SCAM_SNOWD_CLM5F_02.nc",
                              'snowflux')

bulksnow_condfix_clm5 = readanddeseas(filepath+"BULKSNOW/BULKSNOW_clm5conductance_SCAM_CLM5_CLM5F_01.nc",
                              filepath+"BULKSNOW/BULKSNOW_clm5conductance_SCAM_CLM5_CLM5F_02.nc",
                              'snowflux')

bulksnow_condfix_snowd = readanddeseas(filepath+"BULKSNOW/BULKSNOW_clm5conductance_SCAM_SNOWD_CLM5F_01.nc",
                              filepath+"BULKSNOW/BULKSNOW_clm5conductance_SCAM_SNOWD_CLM5F_02.nc",
                              'snowflux')

shflxconstructclm5_1 = xr.open_dataset(filepath+"SHFLXconstruct/SHFLXconstruct_SCAM_CLM5_CLM5F_01.nc")
shflxconstructclm5_2 = xr.open_dataset(filepath+"SHFLXconstruct/SHFLXconstruct_SCAM_CLM5_CLM5F_02.nc")
nyears1 = shflxconstructclm5_1.year.size ; nyears2 = shflxconstructclm5_2.year.size
nyears = nyears1+nyears2 ; ndays = shflxconstructclm5_1.day.size
ncities = shflxconstructclm5_1.city.size
shflxconstructclm5 = np.zeros([nyears, ndays, ncities])
shflxconstructclm5[0:nyears1,:,:] = shflxconstructclm5_1.shflxconstruct
shflxconstructclm5[nyears1:nyears,:,:] = shflxconstructclm5_2.shflxconstruct

shflxconstructsnowd_1 = xr.open_dataset(filepath+"SHFLXconstruct/SHFLXconstruct_SCAM_SNOWD_CLM5F_01.nc")
shflxconstructsnowd_2 = xr.open_dataset(filepath+"SHFLXconstruct/SHFLXconstruct_SCAM_SNOWD_CLM5F_02.nc")
nyears1 = shflxconstructsnowd_1.year.size ; nyears2 = shflxconstructsnowd_2.year.size
nyears = nyears1+nyears2 ; ndays = shflxconstructsnowd_1.day.size
ncities = shflxconstructsnowd_1.city.size
shflxconstructsnowd = np.zeros([nyears, ndays, ncities])
shflxconstructsnowd[0:nyears1,:,:] = shflxconstructsnowd_1.shflxconstruct
shflxconstructsnowd[nyears1:nyears,:,:] = shflxconstructsnowd_2.shflxconstruct


nyears = len(t850_clm5[:,0,0]) ; ndays = len(t850_clm5[0,:,0])


lag = np.arange(-10,11,1)

dat1clm5 = np.reshape(np.array(t850_clm5[:,10:ndays-10,:]),[nyears*(ndays-lag.size+1),3])
dat1snowd = np.reshape(np.array(t850_snowd[:,10:ndays-10,:]),[nyears*(ndays-lag.size+1),3])

t850regclm5 = np.zeros([lag.size,3]) ; t850regsnowd = np.zeros([lag.size,3])
trefhtregclm5 = np.zeros([lag.size,3]) ; trefhtregsnowd = np.zeros([lag.size,3])
tsregclm5 = np.zeros([lag.size,3]) ; tsregsnowd = np.zeros([lag.size,3])
tbotregclm5 = np.zeros([lag.size,3]) ; tbotregsnowd = np.zeros([lag.size,3])
flnsregclm5 = np.zeros([lag.size,3]) ; flnsregsnowd = np.zeros([lag.size,3])
fldsregclm5 = np.zeros([lag.size,3]) ; fldsregsnowd = np.zeros([lag.size,3])
fsnsregclm5 = np.zeros([lag.size,3]) ; fsnsregsnowd = np.zeros([lag.size,3])
shflxregclm5 = np.zeros([lag.size,3]) ; shflxregsnowd = np.zeros([lag.size,3])
lhflxregclm5 = np.zeros([lag.size,3]) ; lhflxregsnowd = np.zeros([lag.size,3])
fgrregclm5 = np.zeros([lag.size,3]) ; fgrregsnowd = np.zeros([lag.size,3])
bulksnowregclm5 = np.zeros([lag.size,3]) ; bulksnowregsnowd = np.zeros([lag.size,3])
bulksnowcondfixregclm5 = np.zeros([lag.size,3]) ; bulksnowcondfixregsnowd = np.zeros([lag.size,3])
shflxconstructregclm5 = np.zeros([lag.size,3]) ; shflxconstructregsnowd = np.zeros([lag.size,3])


for ilag in np.arange(0,lag.size,1):

    #------T850
    dat2clm5 = np.reshape(np.array(t850_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(t850_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        t850regclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        t850regsnowd[ilag,icity] = result[0]
    #------end T850

    #------TREFHT
    dat2clm5 = np.reshape(np.array(trefht_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(trefht_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        trefhtregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        trefhtregsnowd[ilag,icity] = result[0]
    #------end TREFHT

    #------TS
    dat2clm5 = np.reshape(np.array(ts_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(ts_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        tsregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        tsregsnowd[ilag,icity] = result[0]
    #------end TS

    #------TBOT
    dat2clm5 = np.reshape(np.array(tbot_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(tbot_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        tbotregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        tbotregsnowd[ilag,icity] = result[0]
    #------end TBOT

    #------SHFLX
    dat2clm5 = np.reshape(np.array(shflx_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(shflx_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        shflxregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        shflxregsnowd[ilag,icity] = result[0]
    #------end SHFLX

    #------LHFLX
    dat2clm5 = np.reshape(np.array(lhflx_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(lhflx_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        lhflxregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        lhflxregsnowd[ilag,icity] = result[0]
    #------end LHFLX

    #------FGR
    dat2clm5 = np.reshape(np.array(fgr_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(fgr_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        fgrregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        fgrregsnowd[ilag,icity] = result[0]
    #------end FGR

    #------FLNS
    dat2clm5 = np.reshape(np.array(flns_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(flns_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        flnsregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        flnsregsnowd[ilag,icity] = result[0]
    #------end FLNS

    #------FLDS
    dat2clm5 = np.reshape(np.array(flds_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(flds_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        fldsregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        fldsregsnowd[ilag,icity] = result[0]
    #------end FLDS

    #------FSNS
    dat2clm5 = np.reshape(np.array(fsns_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(fsns_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        fsnsregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        fsnsregsnowd[ilag,icity] = result[0]
    #------end FLNS

    #------SHFLX construct
    dat2clm5 = np.reshape(np.array(shflxconstructclm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(shflxconstructsnowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        shflxconstructregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        shflxconstructregsnowd[ilag,icity] = result[0]
    #------end SHFLX construct 



    #------BULKSNOW
    dat2clm5 = np.reshape(np.array(bulksnow_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(bulksnow_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        dat1city = dat1clm5[:,icity] ; dat2city = dat2clm5[:,icity]
        dat1cityuse = dat1city[~np.isnan(dat2city)]
        dat2cityuse = dat2city[~np.isnan(dat2city)]
        result = np.polyfit(dat1cityuse, dat2cityuse,1)
        bulksnowregclm5[ilag,icity] = result[0]

        dat1city = dat1snowd[:,icity] ; dat2city = dat2snowd[:,icity]
        dat1cityuse = dat1city[~np.isnan(dat2city)]
        dat2cityuse = dat2city[~np.isnan(dat2city)]
        result = np.polyfit(dat1cityuse, dat2cityuse,1)
        bulksnowregsnowd[ilag,icity] = result[0]
    #------end BULKSNOW 


    #------BULKSNOW conductance fixed
    dat2clm5 = np.reshape(np.array(bulksnow_condfix_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(bulksnow_condfix_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):

        dat1city = dat1clm5[:,icity] ; dat2city = dat2clm5[:,icity]
        dat1cityuse = dat1city[~np.isnan(dat2city)]
        dat2cityuse = dat2city[~np.isnan(dat2city)]
        result = np.polyfit(dat1cityuse, dat2cityuse,1)
        bulksnowcondfixregclm5[ilag,icity] = result[0]

        dat1city = dat1snowd[:,icity] ; dat2city = dat2snowd[:,icity]
        dat1cityuse = dat1city[~np.isnan(dat2city)]
        dat2cityuse = dat2city[~np.isnan(dat2city)]
        result = np.polyfit(dat1cityuse, dat2cityuse,1)
        bulksnowcondfixregsnowd[ilag,icity] = result[0]
    #------end BULKSNOW 


t850regclm5 = xr.DataArray(t850regclm5, name='t850regclm5',
                 coords=[lag,cities], dims=['lag','city'])
t850regsnowd = xr.DataArray(t850regsnowd, name='t850regsnowd',
                 coords=[lag,cities], dims=['lag','city'])

trefhtregclm5 = xr.DataArray(trefhtregclm5, name='trefhtregclm5',
                 coords=[lag,cities], dims=['lag','city'])
trefhtregsnowd = xr.DataArray(trefhtregsnowd, name='trefhtregsnowd',
                 coords=[lag,cities], dims=['lag','city'])

tbotregclm5 = xr.DataArray(tbotregclm5, name='tbotregclm5',
                 coords=[lag,cities], dims=['lag','city'])
tbotregsnowd = xr.DataArray(tbotregsnowd, name='tbotregsnowd',
                 coords=[lag,cities], dims=['lag','city'])

tsregclm5 = xr.DataArray(tsregclm5, name='tsregclm5',
                 coords=[lag,cities], dims=['lag','city'])
tsregsnowd = xr.DataArray(tsregsnowd, name='tsregsnowd',
                 coords=[lag,cities], dims=['lag','city'])

shflxregclm5 = xr.DataArray(shflxregclm5, name='shflxregclm5',
                 coords=[lag,cities], dims=['lag','city'])
shflxregsnowd = xr.DataArray(shflxregsnowd, name='shflxregsnowd',
                 coords=[lag,cities], dims=['lag','city'])

lhflxregclm5 = xr.DataArray(lhflxregclm5, name='lhflxregclm5',
                 coords=[lag,cities], dims=['lag','city'])
lhflxregsnowd = xr.DataArray(lhflxregsnowd, name='lhflxregsnowd',
                 coords=[lag,cities], dims=['lag','city'])

flnsregclm5 = xr.DataArray(flnsregclm5, name='flnsregclm5',
                 coords=[lag,cities], dims=['lag','city'])
flnsregsnowd = xr.DataArray(flnsregsnowd, name='flnsregsnowd',
                 coords=[lag,cities], dims=['lag','city'])

fldsregclm5 = xr.DataArray(fldsregclm5, name='fldsregclm5',
                 coords=[lag,cities], dims=['lag','city'])
fldsregsnowd = xr.DataArray(fldsregsnowd, name='fldsregsnowd',
                 coords=[lag,cities], dims=['lag','city'])

fsnsregclm5 = xr.DataArray(fsnsregclm5, name='fsnsregclm5',
                 coords=[lag,cities], dims=['lag','city'])
fsnsregsnowd = xr.DataArray(fsnsregsnowd, name='fsnsregsnowd',
                 coords=[lag,cities], dims=['lag','city'])

fgrregclm5 = xr.DataArray(fgrregclm5, name='fgrregclm5',
                 coords=[lag,cities], dims=['lag','city'])
fgrregsnowd = xr.DataArray(fgrregsnowd, name='fgrregsnowd',
                 coords=[lag,cities], dims=['lag','city'])

bulksnowregclm5 = xr.DataArray(bulksnowregclm5, name='bulksnowregclm5',
                 coords=[lag,cities], dims=['lag','city'])
bulksnowregsnowd = xr.DataArray(bulksnowregsnowd, name='bulksnowregsnowd',
                 coords=[lag,cities], dims=['lag','city'])

bulksnowcondfixregclm5 = xr.DataArray(bulksnowcondfixregclm5, name='bulksnowcondfixregclm5',
                 coords=[lag,cities], dims=['lag','city'])
bulksnowcondfixregsnowd = xr.DataArray(bulksnowcondfixregsnowd, name='bulksnowcondfixregsnowd',
                 coords=[lag,cities], dims=['lag','city'])


shflxconstructregclm5 = xr.DataArray(shflxconstructregclm5, name='shflxconstructregclm5',
                 coords=[lag,cities], dims=['lag','city'])
shflxconstructregsnowd = xr.DataArray(shflxconstructregsnowd, name='shflxconstructregsnowd',
                 coords=[lag,cities], dims=['lag','city'])


try:
    os.remove(pathout+"laggedreg_scam.nc")
except:
    pass

t850regclm5.to_netcdf(pathout+"laggedreg_scam.nc")
t850regsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
trefhtregclm5.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
trefhtregsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
tbotregclm5.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
tbotregsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
tsregclm5.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
tsregsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
shflxregclm5.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
shflxregsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
lhflxregclm5.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
lhflxregsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
flnsregclm5.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
flnsregsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
fldsregclm5.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
fldsregsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
fsnsregclm5.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
fsnsregsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
fgrregclm5.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
fgrregsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
bulksnowregclm5.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
bulksnowregsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
bulksnowcondfixregclm5.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
bulksnowcondfixregsnowd.to_netcdf(pathout+"laggedreg_scam.nc",mode="a")
shflxconstructregclm5.to_netcdf(pathout+"laggedreg_scam.nc", mode="a")
shflxconstructregsnowd.to_netcdf(pathout+"laggedreg_scam.nc", mode="a")




