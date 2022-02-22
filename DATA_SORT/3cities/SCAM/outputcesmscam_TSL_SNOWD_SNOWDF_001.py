import importlib
import xarray as xr
import numpy as np
import pandas as pd

from CASutils import filter_utils as filt
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

importlib.reload(filt)
importlib.reload(read)
importlib.reload(cal)


expname=['SASK_SNOWD_SNOWDF_01.001.FSCAM.sasksnowd1',
         'TOR_SNOWD_SNOWDF_01.001.FSCAM.torsnowd1',
         'SID_SNOWD_SNOWDF_01.001.FSCAM.sidsnowd1']

outname='SCAM_SNOWD_SNOWDF_001'

cityname=['Saskatoon','Toronto','Siderovsk']
citylon=[253.330, 280.617, 82.3139]
citylat=[52.1579, 43.6532, 66.5973]

for icity in np.arange(0,3,1):

    basedir="/project/cas02/islas/CLM5_CLM4/raw/SCAM_new_lowrelax/"
    pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/"
    
    fpath=basedir+expname[icity]+"/atm/hist/h0concat.nc"
    print(fpath)
    #dat = read.read_sfc_cesm(fpath,"1979-01-01T12:00:00","2014-12-31T12:00:00")
    dat = read.read_sfc_cesm_dailyavg(fpath,"1979-01-01","2014-12-31")  
 
    if (icity == 0): 
        tsl = xr.DataArray(np.zeros([dat.time.size, 3]), coords=[dat.time, cityname],
                                       dims=['time','city'], name='tsl')
    
    tsl[:,icity] = dat.TSL.isel(lon=0,lat=0)
    
    tsl.to_netcdf(path=pathout+"TSL_"+outname+".nc")