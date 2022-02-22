import importlib
import xarray as xr
import numpy as np
import sys
from glob import glob
import os

from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal
from CASutils import filter_utils as filt

importlib.reload(read)
importlib.reload(cal)
importlib.reload(filt)

path="/project/cas02/islas/CESM1LE/day/TREFHT/1979_2014/"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/deseasonalized_tvar/LENS1_djf_var.nc"
#members = glob(path+"/*/")
#print(members)

members = [os.path.basename(path) for path in sorted(glob(path+'/*'))]

count=0
for imem in members:
    print(imem)
    fpath = path+imem+"/*.nc"
#    fout=pathout+"TVAR_"+iexp+".nc"   
    dat = read.read_sfc(fpath,"1979-01","2014-12") 
    datseason = dat.trefht.groupby('time.dayofyear').mean('time')
    trefht4harm = filt.calc_season_nharm(datseason, 4, dimtime=0)
    anoms = dat.trefht.groupby('time.dayofyear') - trefht4harm
    djfanoms = cal.group_season_daily(anoms,'DJF')
    djfmean = djfanoms.mean('day')
    djfanoms = djfanoms - djfmean

    djfvart = np.var(np.array(djfanoms), axis=(0,1))

    if (count == 0):
        djfvar = xr.DataArray(np.zeros([len(members),dat.lat.size,dat.lon.size]),
                 dims=['member','lat','lon'],
                 coords=[members,dat.lat,dat.lon],
                 name='djfvar')

    djfvar[count,:,:] = djfvart
    count=count+1

djfvar.to_netcdf(path=pathout)
