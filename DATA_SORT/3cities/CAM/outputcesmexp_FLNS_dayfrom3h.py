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

#expname=('Isla_CAM6_CLM5_002','Cecile_CAM6_CLM5','Cecile_CAM6_CLM4','CAM6_CLM5_snowdensity',
#'CAM6_CLM5_snowdensity_002')

#expname=['CAM6_CLM5_snowdensity_002']
expname=[ "CAM6_CLM5_snowdensity" ]
#expname=( 'Isla_CAM6_CLM5_002','CAM6_CLM5_snowdensity_002')

for iexp in np.arange(0,len(expname),1):

    print(expname[iexp])

    basedir="/project/cas02/islas/CLM5_CLM4/raw/"
    pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/"
    
    cityname=['Saskatoon','Toronto','Siderovsk']
    citylon=[253.330, 280.617, 82.3139]
    citylat=[52.1579, 43.6532, 66.5973]
    
    fpath=basedir+expname[iexp]+"/3h/FLNS/*.nc"
    print(fpath)
    dat3h = read.read_sfc_cesm_3hourly(fpath,"1979-01-01","2014-12-31")
    datestr = xr.DataArray( dat3h.indexes['time'].strftime('%Y-%m-%d'), coords = dat3h.time.coords, name='datestr')
   
    dat = dat3h.groupby(datestr).mean('time')
    timeout = dat3h.time.groupby(datestr).mean('time')
    dat = dat.rename({'datestr':'time'})
    dat['time'] = timeout
 
    flns = xr.DataArray(np.zeros([timeout.size, 3]), coords=[timeout, cityname],
                                       dims=['time','city'], name='flns')
    
    for icity in np.arange(0,len(cityname),1):
       flns[:,icity] = dat.FLNS.sel(lon=citylon[icity],lat=citylat[icity], method='nearest')
    
    flns.to_netcdf(path=pathout+"FLNS_"+expname[iexp]+".nc")
