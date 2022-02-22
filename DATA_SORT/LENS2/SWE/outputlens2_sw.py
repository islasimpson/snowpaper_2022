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

filepath="/project/mojave/cesm2/LENS/lnd/month_1/"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/LENS2/SWE/"

nmems=100
memnames = lens.lens2memnamegen(nmems)
memnames1 = memnames[0:10]
memnames2 = memnames[20:100]
memnames = memnames1+memnames2
memnames.remove('1301.017')

for imem in np.arange(0,len(memnames),1):
    print(imem)
    snowice = xr.open_mfdataset(filepath+"SNOWICE/*"+memnames[imem]+"*SNOWICE*.nc")
    snowliq = xr.open_mfdataset(filepath+"SNOWLIQ/*"+memnames[imem]+"*SNOWLIQ*.nc")
    snowice = read.fixcesmtime(snowice, timebndsvar='time_bounds')
    snowliq = read.fixcesmtime(snowliq, timebndsvar='time_bounds')
    snowice = snowice.sel(time=slice("1979-01","2014-12"))
    snowliq = snowliq.sel(time=slice("1979-01","2014-12")) 
    swet = ((snowice.SNOWICE + snowliq.SNOWLIQ)/1000.)*1000.
    swe_djf = cal.season_ts(swet,"DJF")
    if (imem == 0):
        swe = xr.DataArray(np.zeros([len(memnames), swe_djf.time.size, swe_djf.lat.size, swe_djf.lon.size]),
                          dims=['member','time','lat','lon'],
                          coords=[memnames,swe_djf.time, swe_djf.lat, swe_djf.lon],
                          name='swe')

    swe[imem,:,:,:] = np.array(swe_djf[:,:,:])

swe.to_netcdf(pathout+"swe_lens2_djf.nc")

