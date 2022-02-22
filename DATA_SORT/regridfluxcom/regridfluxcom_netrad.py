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
outpath="/project/cas02/islas/FluxCom/netrad/"

ystart=2014
yend=2015

reusewgt = False
wgtfile=outpath+"wgtfile.nc"

for iyear in np.arange(ystart,yend,1):
    print(iyear)
    data = xr.open_dataset(basepath+"Rn.RS_METEO.EBC-NONE.MLM-ALL.METEO-CRUNCEP_v8.720_360.daily."\
            +str(iyear)+".nc", decode_times=False)
    time = pd.date_range(str(iyear)+'-01-01',str(iyear)+'-12-31')
    
    data['time'] = time

    netrad = data.Rn
    regridder = xe.Regridder(netrad, grid_out, 'bilinear', periodic=True, reuse_weights=reusewgt,
                filename=wgtfile)
    netrad_rg = regridder(netrad)

    netrad_rg = netrad_rg.rename("netrad")
    netrad_rg.to_netcdf(path=outpath+"netrad_"+str(iyear)+".nc")

    reusewgt = True  


