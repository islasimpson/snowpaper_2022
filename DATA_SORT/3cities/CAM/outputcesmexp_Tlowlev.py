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

#expname=('Isla_CAM6_CLM5_002','Cecile_CAM6_CLM5','Cecile_CAM6_CLM4','CAM6_CLM5_snowdensity',
#'CAM6_CLM5_snowdensity_002')

expname=['Isla_CAM6_CLM5_002']

for iexp in np.arange(0,len(expname),1):

    print(expname[iexp])

    basedir="/project/cas02/islas/CLM5_CLM4/raw/"
    pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/CAM/"
    
    cityname=['Saskatoon','Toronto','Siderovsk']
    citylon=[253.330, 280.617, 82.3139]
    citylat=[52.1579, 43.6532, 66.5973]
    
    fpath=basedir+expname[iexp]+"/day/Tlowlev/*.nc"
    print(fpath)
    dat = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
    
    trefht = xr.DataArray(np.zeros([dat.time.size, dat.lev.size, 3]), 
                                       coords=[dat.time, dat.lev, cityname],
                                       dims=['time','lev','city'], name='T')
    
    for icity in np.arange(0,len(cityname),1):
       trefht[:,:,icity] = dat.T.sel(lon=citylon[icity],lat=citylat[icity], method='nearest')
    
    trefht.to_netcdf(path=pathout+"tlowlev_"+expname[iexp]+".nc")
