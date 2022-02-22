import xarray as xr
import sys
import pandas as pd
import numpy as np
import xesmf as xe
import warnings

cesmdat = xr.open_dataset(
    '/project/cas/islas/cesmle/mon/clim/1979_2005/zg/zg_hist_001_1979_2005avg.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat)}, {'lon': (['lon'], cesmdat.lon)})


basepath="/project/mojave/observations/FluxCom/"
outpath="/project/cas02/islas/FluxCom/SHFLX/"

ystart=2014
yend=2015

reusewgt = False
wgtfile=outpath+"wgtfile.nc"

for iyear in np.arange(ystart,yend,1):
    print(iyear)
    data = xr.open_dataset(basepath+"H.RS_METEO.EBC-ALL.MLM-ALL.METEO-CRUNCEP_v8.720_360.daily."\
            +str(iyear)+".nc", decode_times=False)
    time = pd.date_range(str(iyear)+'-01-01',str(iyear)+'-12-31')
    
    data['time'] = time

    shflx = data.H
    regridder = xe.Regridder(shflx, grid_out, 'bilinear', periodic=True, reuse_weights=reusewgt,
                filename=wgtfile)
    shflx_rg = regridder(shflx)

    shflx_rg = shflx_rg.rename("shflx")
    shflx_rg.to_netcdf(path=outpath+"SHFLX_"+str(iyear)+".nc")

    reusewgt = True  


