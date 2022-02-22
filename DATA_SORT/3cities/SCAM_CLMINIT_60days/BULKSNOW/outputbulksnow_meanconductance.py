import xarray as xr
import sys
from math import nan
import numpy as np
import pandas as pd

# thermal conductivities (W/m)
lamice=2.29 ; lamair=0.023

expname=['SCAM_CLM5_CLM5F_02','SCAM_SNOWD_CLM5F_02']
#expname=['SCAM_CLM5_CLM5F_01','SCAM_CLM5_CLM5F_02','SCAM_SNOWD_CLM5F_01','SCAM_SNOWD_CLM5F_02']

datpath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM_CLMINIT_60days/"
outpath="/project/cas/islas/python_savs/snowpaper/DATA_SORT/3cities/SCAM_CLMINIT_60days/BULKSNOW/"

# calculate the mean conductance for the CLM5 simulations
snowice = xr.open_dataset(datpath+"SNOWICE_SCAM_CLM5_CLM5F_02.nc")
snowice = snowice.snowice
snowliq = xr.open_dataset(datpath+"SNOWLIQ_SCAM_CLM5_CLM5F_02.nc")
snowliq = snowliq.snowliq
snowdp = xr.open_dataset(datpath+"SNOWDP_SCAM_CLM5_CLM5F_02.nc")
snowdp = snowdp.snowdp
rhosnow = (snowice + snowliq)/snowdp
rhosnow = rhosnow.where(rhosnow != 0, nan)
rhosnow = rhosnow.rename('rhosnow')

# snow conductance
conductance = lamair + (7.75e-5*rhosnow + 1.105e-6*rhosnow**2.)*(lamice - lamair)
conductance = np.where(np.array(snowdp) > 0.05, conductance, nan)
conductance = xr.DataArray(conductance, coords=snowice.coords, name='conductance')

conductancedjf = conductance.where(conductance['time.season'] == "DJF").mean("time", skipna=True)



for iexp in expname:
    print(iexp)
    snowice = xr.open_dataset(datpath+"SNOWICE_"+iexp+".nc")
    snowice = snowice.snowice
    snowliq = xr.open_dataset(datpath+"SNOWLIQ_"+iexp+".nc")
    snowliq = snowliq.snowliq
    snowdp = xr.open_dataset(datpath+"SNOWDP_"+iexp+".nc")
    snowdp = snowdp.snowdp
    snottopl = xr.open_dataset(datpath+"SNOTTOPL_"+iexp+".nc")
    snottopl = snottopl.snottopl
    tsl = xr.open_dataset(datpath+"TSL_"+iexp+".nc")
    tsl = tsl.tsl

    # snow density 
    rhosnow = (snowice + snowliq)/snowdp
    rhosnow = rhosnow.where(rhosnow != 0,nan)
    rhosnow = rhosnow.rename('rhosnow')

    # snow conductance
#    conductance = lamair + (7.75e-5*rhosnow + 1.105e-6*rhosnow**2.)*(lamice - lamair)
#    conductance = conductance.rename('conductance')

    # temperature difference between top snow layer and soil
    tdif = snottopl - tsl
    tdif = tdif.rename('tdif')

    # diagnosed bulk flux across snow
#    snowflux = -1.*conductance*tdif/snowdp
#    snowflux = snowflux.rename('snowflux')
    snowflux = np.zeros([snowice.time.size,snowice.city.size])
    for icity in np.arange(0,snowice.city.size,1):
        snowflux[:,icity] = -1.*conductancedjf[icity]*tdif[:,icity]/snowdp[:,icity]


    # taper the winter season fluxes
    taper = np.zeros([365])
    minfull = 335-30 ; maxfull = 59+30 ; ntaper = 30
    taper[ (np.arange(0,365,1) >= minfull ) | (np.arange(0,365,1) <= maxfull)] = 1
    taper[minfull-ntaper:minfull] = 0.5*(1.-np.cos(np.pi*(np.arange(0,ntaper,1)+0.5)/float(ntaper)))
    taper[maxfull+1:maxfull+1+ntaper] = 1 - 0.5*(1.-np.cos(np.pi*(np.arange(0,ntaper,1)+0.5)/float(ntaper)))

    taper_3d = np.tile(taper,int(snowice.time.size/365)*3)
    taper_3d = np.reshape(taper_3d,[3,snowice.time.size])
    taper_3d = np.moveaxis(taper_3d,1,0)

    snowflux = snowflux*taper_3d
    snowflux = np.where(taper_3d == 0, taper_3d, snowflux)
    # only calculate the snow flux if there's more than 5cm of snow
    snowflux = np.where(np.array(snowdp) > 0.05, snowflux, 0)
    snowflux = xr.DataArray(snowflux, coords=rhosnow.coords, name='snowflux')

    rhosnow.to_netcdf(path=outpath+"BULKSNOW_clm5conductance_"+iexp+".nc")
    conductance.to_netcdf(path=outpath+"BULKSNOW_clm5conductance_"+iexp+".nc", mode="a")
    tdif.to_netcdf(path=outpath+"BULKSNOW_clm5conductance_"+iexp+".nc", mode="a")
    snowflux.to_netcdf(path=outpath+"BULKSNOW_clm5conductance_"+iexp+".nc", mode="a")



