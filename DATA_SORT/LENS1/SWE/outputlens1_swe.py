import xarray as xr
import numpy as np
import importlib
import sys
import warnings

from CASutils import lensread_utils as lens
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

importlib.reload(lens)
importlib.reload(read)
importlib.reload(cal)
warnings.filterwarnings('ignore')

filepath="/project/mojave/cesm1/LENS/lnd/month_1/"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/LENS1/SWE/"
nmems = 40
memnames = lens.lens1memnamegen(nmems)

for imem in np.arange(0,len(memnames),1):
    print(imem)
    snowice_hist = xr.open_mfdataset(filepath+"SNOWICE/b.e11.B20TRC5CNBDRD.f09_g16."+memnames[imem]+"*.nc")
    snowice_rcp = xr.open_mfdataset(filepath+"SNOWICE/b.e11.BRCP85C5CNBDRD.f09_g16."+memnames[imem]+"*.nc")
    snowice = xr.concat([snowice_hist, snowice_rcp], dim='time')

    snowliq_hist = xr.open_mfdataset(filepath+"SNOWLIQ/b.e11.B20TRC5CNBDRD.f09_g16."+memnames[imem]+"*.nc")
    snowliq_rcp = xr.open_mfdataset(filepath+"SNOWLIQ/b.e11.BRCP85C5CNBDRD.f09_g16."+memnames[imem]+"*.nc")
    snowliq = xr.concat([snowliq_hist, snowliq_rcp], dim='time')

    snowice = read.fixcesmtime(snowice, timebndsvar='time_bounds')
    snowliq = read.fixcesmtime(snowliq, timebndsvar='time_bounds')

    snowice = snowice.sel(time=slice("1979-01","2014-12"))
    snowliq = snowliq.sel(time=slice("1979-01","2014-12"))

    swet = (( snowice.SNOWICE + snowliq.SNOWLIQ)/1000.)*1000.
    swe_djf = cal.season_ts(swet, "DJF")
   
    if (imem == 0):
        swe = xr.DataArray(
           np.zeros([len(memnames), swe_djf.time.size, swe_djf.lat.size, swe_djf.lon.size]),
           dims=['member','time','lat','lon'],
           coords=[memnames, swe_djf.time, swe_djf.lat, swe_djf.lon],
           name='swe')

    swe[imem,:,:,:] = np.array(swe_djf[:,:,:])

swe.to_netcdf(pathout+"swe_lens1_djf.nc")
