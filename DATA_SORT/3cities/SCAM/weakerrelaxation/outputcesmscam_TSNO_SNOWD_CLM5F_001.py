import importlib
import xarray as xr
import numpy as np
import pandas as pd
import sys

from CASutils import filter_utils as filt
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

importlib.reload(filt)
importlib.reload(read)
importlib.reload(cal)


expname=['SASK_SNOWD_CLM5F_01.001.FSCAM.sask_1979_2014',
         'TOR_SNOWD_CLM5F_01.001.FSCAM.tor_1979_2014',
         'SID_SNOWD_CLM5F_01.001.FSCAM.sid_1979_2014']

outname='SCAM_SNOWD_CLM5F_01'

cityname=['Saskatoon','Toronto','Siderovsk']
citylon=[253.330, 280.617, 82.3139]
citylat=[52.1579, 43.6532, 66.5973]

for icity in np.arange(0,3,1):

    basedir="/project/cas02/islas/CLM5_CLM4/raw/SCAM_new_lowrelax/"
    pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM/new_lowrelax/"
    
    fpath=basedir+expname[icity]+"/lnd/hist/h2concat.nc"
    dat = xr.open_mfdataset(fpath, coords='minimal', join='override', decode_times=True)
    dat = dat.sel(time=slice("1979-01-01T00:00:00", "2014-12-31T23:50:00"))
    daystr = xr.DataArray(dat.indexes['time'].strftime('%Y-%m-%d'), coords = dat.time.coords, name='daystr')
    snot = dat.SNO_T
    snotdaily = snot.groupby(daystr).mean('time', skipna=True)
    time = dat.time.groupby(daystr).mean('time')
    snotdaily['daystr'] = time
    snotdaily = snotdaily.rename({'daystr':'time'})

    if (icity == 0): 
        snotout = xr.DataArray(np.zeros([snotdaily.time.size, snotdaily.levsno.size, 3]), 
                            coords=[snotdaily.time, snotdaily.levsno, cityname],
                                       dims=['time','levsno','city'], name='snot')
    
    snotout[:,:,icity] = snotdaily.isel(lon=0,lat=0)
    
    snotout.to_netcdf(path=pathout+"SNO_T_"+outname+".nc")
