import xarray as xr
import numpy as np
import sys
import importlib

from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal
from CASutils import filter_utils as filt

from scipy.stats import linregress

importlib.reload(read)
importlib.reload(cal)
importlib.reload(filt)

pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/flux_trefht_gradients/"

path="/project/cas02/islas/CLM5_CLM4/raw/offlineland/"
expnames=['offlineland_clm5','offlineland_snowrevert']

def deseasonalized(data):
    datseason = data.groupby('time.dayofyear').mean('time')
    dat4harm = filt.calc_season_nharm(datseason,4,dimtime=0)
    datanoms = data.groupby('time.dayofyear')-dat4harm
    datanoms = cal.group_season_daily(datanoms,'DJF')
    datmean = datanoms.mean('day')
    datanoms = datanoms - datmean
    return datanoms

for iexp in expnames:
    fpath=path+iexp+"/TSA/*.nc"
    trefht = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
    trefhtanoms = deseasonalized(trefht.TSA)

    fpath=path+iexp+"/FSH/*.nc"
    shflx = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
    shflxanoms = deseasonalized(shflx.FSH)

    fpath=path+iexp+"/FIRA/*.nc"
    flns = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
    flnsanoms = deseasonalized(flns.FIRA)

    trefhtstacked = trefhtanoms.stack(time=('year','day'))
    shflxstacked = shflxanoms.stack(time=('year','day'))
    flnsstacked = flnsanoms.stack(time=('year','day'))

    bshflx, ashflx, r, p, stderr = xr.apply_ufunc(
      linregress, trefhtstacked, shflxstacked, input_core_dims=[['time'],['time']],
      output_core_dims=[[],[],[],[],[]], vectorize=True)

    bflns, aflns, r, p, stderr = xr.apply_ufunc(
      linregress, trefhtstacked, flnsstacked, input_core_dims=[['time'],['time']],
      output_core_dims=[[],[],[],[],[]], vectorize=True)

    bshflx = bshflx.rename('bshflx')
    bflns = bflns.rename('bflns')

    bshflx.to_netcdf(path=pathout+"gradients_"+iexp+".nc")
    bflns.to_netcdf(path=pathout+"gradients_"+iexp+".nc", mode="a")



