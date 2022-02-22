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

    basedir="/project/cas02/islas/CLM5_CLM4/raw/SCAM_CLM_INIT_60days/"
    pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM_CLMINIT_60days/"
    
    fpath=basedir+expname[icity]+"/lnd/hist/h1concat.nc"
    print(fpath)
    dat = read.read_sfc_cesm(fpath,"1979-01-01T12:00:00","2014-12-31T12:00:00")
   
    if (icity == 0): 
        snottopl = xr.DataArray(np.zeros([dat.time.size, 3]), coords=[dat.time, cityname],
                                       dims=['time','city'], name='snottopl')
    
    snottopl[:,icity] = dat.SNOTTOPL.isel(lon=0,lat=0)
    
    snottopl.to_netcdf(path=pathout+"SNOTTOPL_"+outname+".nc")
