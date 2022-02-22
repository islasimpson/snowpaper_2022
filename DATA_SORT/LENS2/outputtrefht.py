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

filepath="/project/mojave/cesm2/LENS/atm/month_1/TREFHT/"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/LENS2/"

#nmems=50
#memnames = lens.lens2memnamegen_first50(nmems)

nmems=100
memnames = lens.lens2memnamegen(nmems)

memnames1 = memnames[0:10]
memnames2 = memnames[20:100]
memnames = memnames1+memnames2
memnames.remove('1301.017')

#sys.exit()
#memnames.remove('1171.009')
#memnames.remove('1191.010')
#memnames.remove('1111.002')
#memnames.remove('1121.003')
#memnames.remove('1131.004')


for imem in np.arange(0,len(memnames),1):
    print(imem)
    trefhtt = read.read_sfc_cesm(filepath+"*-"+memnames[imem]+"*.nc","1920-01","2100-12")
    trefht_djf = cal.season_ts(trefhtt.TREFHT,"DJF")
    if (imem == 0):
        trefht = xr.DataArray(
         np.zeros([len(memnames), trefht_djf.time.size, trefht_djf.lat.size, trefht_djf.lon.size]),
         dims=['member','time','lat','lon'], 
         coords=[memnames, trefht_djf.time, trefht_djf.lat, trefht_djf.lon], name='trefht')
    trefht[imem,:,:,:] = trefht_djf[:,:,:]

trefht.to_netcdf(pathout+"trefht_lens2_djf.nc")

