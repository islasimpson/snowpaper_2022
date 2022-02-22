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


expname=['SASK_SNOWD_SNOWDF_02.001.FSCAM.sasksnowd2',
         'TOR_SNOWD_SNOWDF_02.001.FSCAM.torsnowd2',
         'SID_SNOWD_SNOWDF_02.001.FSCAM.sidsnowd2']

outname='SCAM_SNOWD_SNOWDF_02'

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
    sno = dat.SNO_Z
    snodaily = sno.groupby(daystr).mean('time', skipna=True)
    time = dat.time.groupby(daystr).mean('time')
    snodaily['daystr'] = time
    snodaily = snodaily.rename({'daystr':'time'})

    if (icity == 0): 
        snoout = xr.DataArray(np.zeros([snodaily.time.size, snodaily.levsno.size, 3]), 
                            coords=[snodaily.time, snodaily.levsno, cityname],
                                       dims=['time','levsno','city'], name='snoz')
    
    snoout[:,:,icity] = snodaily.isel(lon=0,lat=0)
   
    snoout.to_netcdf(path=pathout+"SNO_Z_"+outname+".nc")
