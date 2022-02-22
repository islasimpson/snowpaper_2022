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


expname=['SASK_CLM5_CLM5F_02.001.FSCAM.sask2_1979_2014',
         'TOR_CLM5_CLM5F_02.001.FSCAM.tor2_1979_2014',
         'SID_CLM5_CLM5F_02.001.FSCAM.sid2_1979_2014']

outname='SCAM_CLM5_CLM5F_02'

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
    snoz = dat.SNO_Z
    snozdaily = snoz.groupby(daystr).mean('time', skipna=True)
    time = dat.time.groupby(daystr).mean('time')
    snozdaily['daystr'] = time
    snozdaily = snozdaily.rename({'daystr':'time'})

    if (icity == 0): 
        snozout = xr.DataArray(np.zeros([snozdaily.time.size, snozdaily.levsno.size, 3]), 
                            coords=[snozdaily.time, snozdaily.levsno, cityname],
                                       dims=['time','levsno','city'], name='snoz')
    
    snozout[:,:,icity] = snozdaily.isel(lon=0,lat=0)
    
    snozout.to_netcdf(path=pathout+"SNO_Z_"+outname+".nc")
