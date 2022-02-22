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

filepath="/project/cas02/islas/CESM1LE/"

pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/LENS1/"

nmems = 40
memnames = lens.lens1memnamegen(nmems)

for imem in np.arange(0,len(memnames),1):
    print(imem)
    trefhtt=read.read_sfc_cesm(filepath+"TREFHT/TREFHT_"+memnames[imem]+"*.nc","1920-01","2100-12")
    trefht_djf = cal.season_ts(trefhtt.TREFHT,"DJF")
    if (imem == 0):
        trefht = xr.DataArray(
         np.zeros([len(memnames),trefht_djf.time.size, trefhtt.lat.size, trefhtt.lon.size]),
         dims=['member','time','lat','lon'], 
         coords=[memnames,trefht_djf.time, trefht_djf.lat, trefht_djf.lon], name='trefht')

    trefht[imem,:,:,:] = np.array(trefht_djf[:,:,:])
trefht.to_netcdf(pathout+"trefht_lens1_djf.nc")


