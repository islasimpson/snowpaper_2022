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
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/LENS1/BULKSNOW/"

nmems=40
memnames = lens.lens1memnamegen(nmems)

for imem in np.arange(0,len(memnames),1):
    print(imem)
    snowice = read.read_sfc_cesm(filepath+"SNOWICE/SNOWICE_"+memnames[imem]+"*.nc","1850-01","2100-12") 
    snowliq = read.read_sfc_cesm(filepath+"SNOWLIQ/SNOWLIQ_"+memnames[imem]+"*.nc","1850-01","2100-12")
    snowdp = read.read_sfc_cesm(filepath+"SNOWDP/SNOWDP_"+memnames[imem]+"*.nc","1850-01","2100-12")
    rhosnowt = (snowice.SNOWICE + snowliq.SNOWLIQ) / snowdp.SNOWDP
    rhosnow_djf = cal.season_ts(rhosnowt,"DJF")
    if (imem == 0):
        rhosnow = xr.DataArray(np.zeros([len(memnames),rhosnow_djf.time.size, rhosnow_djf.lat.size, rhosnow_djf.lon.size]),
                              dims=['member','time','lat','lon'], coords = [memnames,rhosnow_djf.time, rhosnow_djf.lat, rhosnow_djf.lon], name='rhosnow')
        
    rhosnow[imem,:,:,:] = rhosnow_djf[:,:,:]

rhosnow.to_netcdf(pathout+"rhosnow_lens1_djf.nc")
