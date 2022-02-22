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

def readanddeseas(file1, var):
    dat1 = xr.open_dataset(file1).sel(time=slice("1979-01-01","2014-12-31"))
    dat1 = dat1[var]
    dat1deseas = deseasonalize(dat1)

    return dat1deseas

#def readanddeseas(file1,var):
#    dat = xr.open_dataset(file1).sel(time=slice("1979-01-01","2014-12-31"))
#    dat = dat[var]
#    datdeseas = deseasonalize(dat)
#    return datdeseas

cities=['Saskatoon','Toronto','Siderovsk']

filepath = "/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM_CLMINIT_60days_withclearsky/"

t850_clm5 = readanddeseas(filepath+"T850_SCAM_CLM5_CLM5F_01.nc","t850")

t850_snowd = readanddeseas(filepath+"T850_SCAM_SNOWD_CLM5F_01.nc","t850")

trefht_clm5 = readanddeseas(filepath+"TREFHT_SCAM_CLM5_CLM5F_01.nc","t850")

trefht_snowd = readanddeseas(filepath+"TREFHT_SCAM_SNOWD_CLM5F_01.nc","t850")

flns_clm5 = readanddeseas(filepath+"FLNS_SCAM_CLM5_CLM5F_01.nc","flns")

flns_snowd = readanddeseas(filepath+"FLNS_SCAM_SNOWD_CLM5F_01.nc","flns")

flnsc_clm5 = readanddeseas(filepath+"FLNSC_SCAM_CLM5_CLM5F_01.nc","flnsc")

flnsc_snowd = readanddeseas(filepath+"FLNSC_SCAM_SNOWD_CLM5F_01.nc","flnsc")

nyears = len(t850_clm5[:,0,0]) ; ndays = len(t850_clm5[0,:,0])


lag = np.arange(-10,11,1)

dat1clm5 = np.reshape(np.array(t850_clm5[:,10:ndays-10,:]),[nyears*(ndays-lag.size+1),3])
dat1snowd = np.reshape(np.array(t850_snowd[:,10:ndays-10,:]),[nyears*(ndays-lag.size+1),3])

trefhtregclm5 = np.zeros([lag.size,3]) ; trefhtregsnowd = np.zeros([lag.size,3])
t850regclm5 = np.zeros([lag.size,3]) ; t850regsnowd = np.zeros([lag.size,3])
flnsregclm5 = np.zeros([lag.size,3]) ; flnsregsnowd = np.zeros([lag.size,3])
flnscregclm5 = np.zeros([lag.size,3]) ; flnscregsnowd = np.zeros([lag.size,3])

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

    #------FLNS
    dat2clm5 = np.reshape(np.array(flns_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(flns_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        flnsregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        flnsregsnowd[ilag,icity] = result[0]
    #------end FLNS


    #------FLNSC
    dat2clm5 = np.reshape(np.array(flnsc_clm5[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])
    dat2snowd = np.reshape(np.array(flnsc_snowd[:,10+lag[ilag]:ndays-10+lag[ilag],:]),[nyears*(ndays-lag.size+1),3])

    for icity in np.arange(0,3,1):
        result = np.polyfit(dat1clm5[:,icity], dat2clm5[:,icity],1)
        flnscregclm5[ilag,icity] = result[0]
        result = np.polyfit(dat1snowd[:,icity], dat2snowd[:,icity],1)
        flnscregsnowd[ilag,icity] = result[0]
    #------end FLNSC



t850regclm5 = xr.DataArray(t850regclm5, name='t850regclm5',
                 coords=[lag,cities], dims=['lag','city'])
t850regsnowd = xr.DataArray(t850regsnowd, name='t850regsnowd',
                 coords=[lag,cities], dims=['lag','city'])

flnsregclm5 = xr.DataArray(flnsregclm5, name='flnsregclm5',
                 coords=[lag,cities], dims=['lag','city'])
flnsregsnowd = xr.DataArray(flnsregsnowd, name='flnsregsnowd',
                 coords=[lag,cities], dims=['lag','city'])

flnscregclm5 = xr.DataArray(flnscregclm5, name='flnscregclm5',
                 coords=[lag,cities], dims=['lag','city'])
flnscregsnowd = xr.DataArray(flnscregsnowd, name='flnscregsnowd',
                 coords=[lag,cities], dims=['lag','city'])


try:
    os.remove(pathout+"laggedreg_scam_withclearsky.nc")
except:
    pass

t850regclm5.to_netcdf(pathout+"laggedreg_scam_withclearsky.nc")
t850regsnowd.to_netcdf(pathout+"laggedreg_scam_withclearsky.nc",mode="a")
flnsregclm5.to_netcdf(pathout+"laggedreg_scam_withclearsky.nc",mode="a")
flnsregsnowd.to_netcdf(pathout+"laggedreg_scam_withclearsky.nc",mode="a")
flnscregclm5.to_netcdf(pathout+"laggedreg_scam_withclearsky.nc",mode="a")
flnscregsnowd.to_netcdf(pathout+"laggedreg_scam_withclearsky.nc",mode="a")
