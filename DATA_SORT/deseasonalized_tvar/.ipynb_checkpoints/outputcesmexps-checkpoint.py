import importlib
import xarray as xr
import numpy as np

from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal
from CASutils import filter_utils as filt

importlib.reload(read)
importlib.reload(cal)
importlib.reload(filt)

path="/project/cas02/islas/CLM5_CLM4/raw/"
pathout="/project/cas/islas/python_savs/snowpaper/DATA_SORT/deseasonalized_tvar/"

#expnames=['Cecile_CAM6_CLM4','Cecile_CAM5_CLM4','Cecile_CAM5_CLM5','Isla_CAM6_CLM4p5','Cecile_CAM6_CLM4p5']
#expnames=['Cecile_CAM6_CLM5']
#expnames=['mam3.001','mg1.002',
# 'sb.002','tms.001','uw.003']
expnames=['CAM6_CLM5_snowdensity','CAM6_CLM5_snowdensity_002']


for iexp in expnames:
    print(iexp)
    #fpath=path+"f.e20.FHIST.f09_f09.cesm2_1_"+iexp+"/TREFHT/*.nc"
    fpath = path+iexp+"/day/TREFHT/*.nc"
    print(fpath)
    fout=pathout+"TVAR_"+iexp+".nc"
    dat = read.read_sfc_cesm(fpath,"1979-01","2014-12")
    datseason = dat.TREFHT.groupby('time.dayofyear').mean('time')
    trefht4harm = filt.calc_season_nharm(datseason, 4, dimtime=0)
    anoms = dat.TREFHT.groupby('time.dayofyear') - trefht4harm
    djfanoms = cal.group_season_daily(anoms,'DJF')
    djfmean = djfanoms.mean('day')
    djfanoms = djfanoms - djfmean

    djfvar = np.var(np.array(djfanoms), axis=(0,1))

    djfvar = xr.DataArray(djfvar, coords=[anoms.lat, anoms.lon], dims = ['lat','lon'], name='djfvar')
    djfvar.to_netcdf(path=fout)

print("******DONE**********")
