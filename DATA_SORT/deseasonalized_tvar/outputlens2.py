import importlib
import xarray as xr
import numpy as np
import sys
from glob import glob
import os

from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal
from CASutils import filter_utils as filt
from CASutils import lensread_utils as lens

importlib.reload(read)
importlib.reload(cal)
importlib.reload(filt)

path="/project/mojave/cesm2/LENS/atm/day_1/TREFHT/"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/deseasonalized_tvar/LENS2_djf_var.nc"
#members = glob(path+"/*/")
#print(members)

#members = [os.path.basename(path) for path in sorted(glob(path+'/*'))]
members = lens.lens2memnamegen(100)
members1 = members[0:10]
members2 = members[20:100]

members = members1+members2
members.remove('1301.017')


count=0
for imem in members:
    print(imem)
    fpath = path+"*-"+imem+"*.nc"
#    fout=pathout+"TVAR_"+iexp+".nc"   
    dat = read.read_sfc(fpath,"1979-01","2014-12") 
    datseason = dat.TREFHT.groupby('time.dayofyear').mean('time')
    trefht4harm = filt.calc_season_nharm(datseason, 4, dimtime=0)
    anoms = dat.TREFHT.groupby('time.dayofyear') - trefht4harm
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
