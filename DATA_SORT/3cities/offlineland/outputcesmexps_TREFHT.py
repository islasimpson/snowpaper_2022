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

expname=['offlineland_clm5','offlineland_snowrevert']

for iexp in np.arange(0,len(expname),1):
    print(expname[iexp])
    basedir="/project/cas02/islas/CLM5_CLM4/raw/offlineland/"
    pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/offlineland/"

    cityname=['Saskatoon','Toronto','Siderovsk']
    citylon=[253.330, 280.617, 82.3139]
    citylat=[52.1579, 43.6532, 66.5973]

    fpath=basedir+expname[iexp]+"/TSA/*.nc"
    print(fpath)
    dat = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
    trefht = xr.DataArray(np.zeros([dat.time.size, 3]), coords=[dat.time, cityname],
                       dims=['time','city'], name='trefht')

    for icity  in np.arange(0,len(cityname),1):
        trefht[:,icity] = dat.TSA.sel(lon=citylon[icity], lat=citylat[icity], method='nearest')

    trefht.to_netcdf(path=pathout+"TSA_"+expname[iexp]+".nc")
