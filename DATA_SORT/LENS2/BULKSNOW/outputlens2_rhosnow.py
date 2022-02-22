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

filepath="/project/cas02/islas/CESM2LE/mon/1850_2100/"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/LENS2/BULKSNOW/"

#nmems=50
#memnames = lens.lens2memnamegen_first50(nmems)

nmems = 100
memnames = lens.lens2memnamegen(nmems)
memnames1 = memnames[0:10]
memnames2 = memnames[20:100]
memnames = memnames1+memnames2
memnames.remove('1301.017')

for imem in np.arange(0,len(memnames),1):
    print(imem)
    snowice = read.read_sfc_cesm(filepath+"SNOWICE/SNOWICE_*"+memnames[imem]+".nc","1850-01","2100-12") 
    snowliq = read.read_sfc_cesm(filepath+"SNOWLIQ/SNOWLIQ_*"+memnames[imem]+".nc","1850-01","2100-12")
    snowdp = read.read_sfc_cesm(filepath+"SNOWDP/SNOWDP_*"+memnames[imem]+".nc","1850-01","2100-12")
    rhosnowt = (snowice.SNOWICE + snowliq.SNOWLIQ) / snowdp.SNOWDP
    rhosnow_djf = cal.season_ts(rhosnowt,"DJF")
    if (imem == 0):
        rhosnow = xr.DataArray(np.zeros([len(memnames),rhosnow_djf.time.size, rhosnow_djf.lat.size, rhosnow_djf.lon.size]),
                              dims=['member','time','lat','lon'], coords = [memnames,rhosnow_djf.time, rhosnow_djf.lat, rhosnow_djf.lon], name='rhosnow')
        
    rhosnow[imem,:,:,:] = np.array(rhosnow_djf[:,:,:])

rhosnow.to_netcdf(pathout+"rhosnow_lens2_djf.nc")
