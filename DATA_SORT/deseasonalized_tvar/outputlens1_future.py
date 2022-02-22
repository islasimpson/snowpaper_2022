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

importlib.reload(lens)
importlib.reload(read)
importlib.reload(cal)
importlib.reload(filt)

path="/project/mojave/cesm1/LENS/atm/day_1/TREFHT/"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/deseasonalized_tvar/LENS1_djf_var_future.nc"

ystart=2060 ; yend = 2089

nmems=40
memnames = lens.lens1memnamegen(nmems)

count=0
for imem in np.arange(0,len(memnames),1):
    dat = read.read_sfc_cesm(path+"b.e11.BRCP85C5CNBDRD.f09_g16."+memnames[imem]+"*",
               str(ystart)+"-01",str(yend)+"-12")

    datseason = dat.TREFHT.groupby('time.dayofyear').mean('time')
    trefht4harm = filt.calc_season_nharm(datseason,4,dimtime=0)
    anoms = dat.TREFHT.groupby('time.dayofyear') - trefht4harm
    djfanoms = cal.group_season_daily(anoms,'DJF')
    djfmean = djfanoms.mean('day')
    djfanoms = djfanoms - djfmean

    djfvart = np.var(np.array(djfanoms), axis=(0,1))
    if (count == 0):
        djfvar = xr.DataArray(np.zeros([len(memnames),dat.lat.size,dat.lon.size]),
          dims=['member','lat','lon'],
          coords=[memnames, dat.lat, dat.lon],
          name='djfvar')

    djfvar[count,:,:] = djfvart
    count=count+1

djfvar.to_netcdf(path=pathout)

