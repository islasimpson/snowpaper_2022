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

    basedir="/project/cas02/islas/CLM5_CLM4/raw/SCAM_CLM_INIT_60days/"
    pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM_CLMINIT_60days/"
    
    fpath=basedir+expname[icity]+"/atm/hist/h0concat.nc"
    print(fpath)
    dat = read.read_sfc_cesm(fpath,"1979-01-01T12:00:00","2014-12-31T12:00:00")
   
    if (icity == 0): 
        t850 = xr.DataArray(np.zeros([dat.time.size, 3]), coords=[dat.time, cityname],
                                       dims=['time','city'], name='t850')
    
    t850[:,icity] = dat.T850.isel(lon=0,lat=0)
    
    t850.to_netcdf(path=pathout+"T850_"+outname+".nc")
