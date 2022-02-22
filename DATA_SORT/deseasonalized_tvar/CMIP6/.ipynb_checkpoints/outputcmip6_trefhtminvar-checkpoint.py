import importlib
import pandas as pd
import xarray as xr
import numpy as np
from numpy import nan
import sys
import warnings
import xesmf as xe
from glob import glob

from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal
from CASutils import filter_utils as filt

importlib.reload(read)
importlib.reload(cal)

warnings.filterwarnings('ignore')

cesmdat = xr.open_dataset(
    '/project/cas/islas/cesmle/mon/clim/1979_2005/zg/zg_hist_001_1979_2005avg.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat)}, {'lon': (['lon'], cesmdat.lon)})


histpath="/project/cmip6/historical/day/"
var="tasmin"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/deseasonalized_tvar/CMIP6/"

cmip6models = pd.read_csv('cmip6csvinfo.csv')

ybegp = 1979 ; yendp = 2014 ; nyears=yendp-ybegp+1 # dates for Past period


models=cmip6models['Model']
nmodels=models.size

wgtfile = pathout+"wgtfile.nc"

def deseas(dat):
    daystr = xr.DataArray(dat.indexes['time'].strftime('%m-%d'), coords = dat.time.coords, name="daystr")
    datseason = dat.groupby(daystr).mean('time')
    dat4harm = filt.calc_season_nharm(datseason,4,dimtime=0)
    anoms = dat.groupby(daystr)-dat4harm
    djfanoms = cal.group_season_daily(anoms,'DJF')
    djfmean = djfanoms.mean('day')
    djfanoms = djfanoms - djfmean
    return djfanoms

for index, modname in models.iteritems():
    print(modname)
    memstr="r1i1p1f1" 
    histdir = glob(histpath+var+"/"+modname+"/"+memstr+"/*/")
    histdir = histdir[0]
    dat = read.read_sfc(histdir+"*.nc",str(ybegp)+'-01-01',str(yendp)+'-12-31')
    dat = dat.sel(time=~((dat.time.dt.month==2) & (dat.time.dt.day == 29)))
    if (dat.time.size != (nyears*365)):
        print("Something's wrong, ndays="+str(dat.time.size)+", expected"+str(nyears*365))

    djfanoms = deseas(dat[var])
    datvar = np.var(djfanoms, axis=(0,1))

    regridder = xe.Regridder(datvar, grid_out, 'bilinear', periodic=True, reuse_weights=False,
      filename=wgtfile)
    datvar_rg = regridder(datvar)
    datvar_rg = datvar_rg.rename('tasminvar')
    datvar_rg.to_netcdf(path=pathout+"tasminvar_"+modname+".nc")
