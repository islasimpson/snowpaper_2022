import importlib
import xarray as xr
import numpy as np
import sys

from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal
from CASutils import filter_utils as filt

from scipy.stats import linregress

importlib.reload(read)
importlib.reload(cal)
importlib.reload(filt)

path="/project/cas02/islas/CLM5_CLM4/raw/"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/flux_trefht_gradients/"


#expnames=["Isla_CAM6_CLM5_002"]
expnames=["Isla_CAM6_CLM5_002","CAM6_CLM5_snowdensity_002"]

def deseasonalized(data):
    datseason = data.groupby('time.dayofyear').mean('time')
    dat4harm = filt.calc_season_nharm(datseason,4,dimtime=0)
    datanoms = data.groupby('time.dayofyear')-dat4harm
    datanoms = cal.group_season_daily(datanoms,'DJF')
    datmean = datanoms.mean('day')
    datanoms = datanoms - datmean
    return datanoms


for iexp in expnames:
   print(iexp)

   fpath=path+iexp+"/day/TREFHT/*.nc"
   trefht = read.read_sfc_cesm(fpath,"1979-01","2014-12")
   trefhtanoms = deseasonalized(trefht.TREFHT)
   del (trefht)
   print("after TREFHT")

   fpath=path+iexp+"/day/SHFLX/*.nc"
   shflx = read.read_sfc_cesm(fpath,"1979-01","2014-12")
   shflxanoms = deseasonalized(shflx.SHFLX)
   del (shflx)
   print("after SHFLX")

   fpath=path+iexp+"/day/LHFLX/*.nc"
   lhflx = read.read_sfc_cesm(fpath,"1979-01","2014-12")
   lhflxanoms = deseasonalized(lhflx.LHFLX)
   del (lhflx)
   print("after LHFLX")

   fpath=path+iexp+"/day/FLNS/*.nc"
   flns = read.read_sfc_cesm(fpath,"1979-01","2014-12")
   flnsanoms = deseasonalized(flns.FLNS)
   del (flns)
   print("after FLNS")

   fpath=path+iexp+"/day/FSNS/*.nc"
   fsns = read.read_sfc_cesm(fpath,"1979-01","2014-12")
   fsnsanoms = deseasonalized(fsns.FSNS)
   del (fsns)
   print("after FSNS")

#   fpath=path+iexp+"/day/TBOT/*.nc"
#   tbot = read.read_sfc_cesm(fpath,"1979-01","2014-12")
#   tbotanoms = deseasonalized(tbot.TBOT)
#   del (tbot)
#   print("after TBOT")
#
#   fpath=path+iexp+"/day/TS/*.nc"
#   ts = read.read_sfc_cesm(fpath,"1979-01","2014-12")
#   tsanoms = deseasonalized(ts.TS)
#   del (ts)
#   print("after TS")

   fpath=path+iexp+"/day/T850/*.nc"
   t850 = read.read_sfc_cesm(fpath,"1979-01","2014-12")
   t850anoms = deseasonalized(t850.T850)
   del (t850)
   print("after T850")

#   tdif = tsanoms - tbotanoms
   sumflux = -1.*fsnsanoms + flnsanoms + shflxanoms + lhflxanoms
   netrad = -1.*fsnsanoms + flnsanoms

#   trefhtstacked = trefhtanoms.stack(time=('year','day'), space=('lat','lon'))
#   shflxstacked = shflxanoms.stack(time=('year','day'), space=('lat','lon'))
#   lhflxstacked = lhflxanoms.stack(time=('year','day'), space=('lat','lon'))
#   flnsstacked = flnsanoms.stack(time=('year','day'), space=('lat','lon'))
#   fsnsstacked = fsnsanoms.stack(time=('year','day'), space=('lat','lon'))

   trefhtstacked = trefhtanoms.stack(time=('year','day'))
   shflxstacked = shflxanoms.stack(time=('year','day'))
   lhflxstacked = lhflxanoms.stack(time=('year','day'))
   flnsstacked = flnsanoms.stack(time=('year','day'))
   fsnsstacked = fsnsanoms.stack(time=('year','day'))
#   tdifstacked = tdif.stack(time=('year','day'))
   sumfluxstacked = sumflux.stack(time=('year','day'))
   netradstacked = netrad.stack(time=('year','day'))
   t850stacked = t850anoms.stack(time=('year','day'))

   bshflx, ashflx, r, p, stderr = xr.apply_ufunc(
      linregress, trefhtstacked, shflxstacked, input_core_dims=[['time'],['time']],
      output_core_dims=[[],[],[],[],[]], vectorize=True)

   blhflx, alhflx, r, p, stderr = xr.apply_ufunc(
      linregress, trefhtstacked, lhflxstacked, input_core_dims=[['time'],['time']],
      output_core_dims=[[],[],[],[],[]], vectorize=True)

   bflns, aflns, r, p, stderr = xr.apply_ufunc(
      linregress, trefhtstacked, flnsstacked, input_core_dims=[['time'],['time']],
      output_core_dims=[[],[],[],[],[]], vectorize=True)

   bfsns, afsns, r, p, stderr = xr.apply_ufunc(
      linregress, trefhtstacked, fsnsstacked, input_core_dims=[['time'],['time']],
      output_core_dims=[[],[],[],[],[]], vectorize=True)
  
#   btdif, atdif, r, p, stderr = xr.apply_ufunc(
#      linregress, trefhtstacked, tdifstacked, input_core_dims=[['time'],['time']],
#      output_core_dims=[[],[],[],[],[]], vectorize=True)

   bsumflux, asumflux, r, p, stderr = xr.apply_ufunc(
      linregress, trefhtstacked, sumfluxstacked, input_core_dims=[['time'],['time']],
      output_core_dims=[[],[],[],[],[]], vectorize=True)

   bnetrad, anetrad, r, p, stderr = xr.apply_ufunc(
      linregress, trefhtstacked, netradstacked, input_core_dims=[['time'],['time']],
      output_core_dims=[[],[],[],[],[]], vectorize=True)

   bt850, at850, r, p, stderr = xr.apply_ufunc(
      linregress, trefhtstacked, t850stacked, input_core_dims=[['time'],['time']],
      output_core_dims=[[],[],[],[],[]], vectorize=True)

#   bshflx = bshflx_stacked.unstack("space")
#   blhflx = blhflx_stacked.unstack("space")
#   bflns = bflns_stacked.unstack("space")
#   bfsns = bfsns_stacked.unstack("space")

   bshflx = bshflx.rename("bshflx")
   blhflx = blhflx.rename("blhflx")
   bflns = bflns.rename("bflns")
   bfsns = bfsns.rename("bfsns")
#   btdif = btdif.rename("btdif")
   bsumflux = bsumflux.rename("bsumflux")
   bnetrad = bnetrad.rename('bnetrad')
   bt850 = bt850.rename("bt850")

   bshflx.to_netcdf(path=pathout+"gradients_"+iexp+".nc")
   blhflx.to_netcdf(path=pathout+"gradients_"+iexp+".nc", mode="a")
   bflns.to_netcdf(path=pathout+"gradients_"+iexp+".nc",mode="a")
   bfsns.to_netcdf(path=pathout+"gradients_"+iexp+".nc",mode="a")
#   btdif.to_netcdf(path=pathout+"gradients_"+iexp+".nc",mode="a")
   bsumflux.to_netcdf(path=pathout+"gradients_"+iexp+".nc", mode="a")
   bnetrad.to_netcdf(path=pathout+'gradients_'+iexp+'.nc', mode='a')
   bt850.to_netcdf(path=pathout+"gradients_"+iexp+".nc", mode="a")
