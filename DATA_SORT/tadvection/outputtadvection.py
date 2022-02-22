import importlib
import xarray as xr
import numpy as np
import sys
import math

from CASutils import readdata_utils as read
from CASutils import filter_utils as filt

importlib.reload(read)

expnames=['Isla_CAM6_CLM5_002','CAM6_CLM5_snowdensity_002']

path="/project/cas02/islas/CLM5_CLM4/raw/"

for iexp in expnames:
    print(iexp)
    outdir=path+iexp+'/day/CLIM_ADVECTION/'
    fpath=path+iexp+'/day/TBOT/*.nc'
    tbot = read.read_sfc_cesm(fpath,"1979-01-01","2014-12-31")
    tbotseason = tbot.TBOT.groupby('time.dayofyear').mean('time')
    tbot4harm = filt.calc_season_nharm(tbotseason, 4, dimtime=0)

    del tbot

    lonrad = np.array((tbot4harm.lon/180.)*math.pi) ; latrad = np.array((tbot4harm.lat/180.)*math.pi)
    a = 6.371e6

    # calculate the climatological temperature gradient

    dtdx = np.zeros([365,latrad.size,lonrad.size])
    dtdy = np.zeros([365,latrad.size,lonrad.size])

    nlon = lonrad.size ; nlat = latrad.size

    for iday in np.arange(0,tbot4harm.dayofyear.size,1):
        for ilat in np.arange(0,latrad.size,1):
            tdat = np.array(tbot4harm.isel(dayofyear=iday, lat=ilat))
            tdatpad = np.tile(tdat,3)
            lonpad = np.tile(lonrad,3)
            lonpad[nlon:2*nlon] = lonpad[nlon:2*nlon]+2*math.pi
            lonpad[2*nlon:3*nlon] = lonpad[2*nlon:3*nlon] + 4*math.pi

            dtdx[iday,ilat,:] = (1./(a*np.cos(latrad[ilat])))* \
              (tdatpad[nlon+5:2*nlon+5] - tdatpad[nlon-5:2*nlon-5])/ \
              (lonpad[nlon+5:2*nlon+5] - lonpad[nlon-5:2*nlon-5])

        for ilon in np.arange(0,lonrad.size,1):
            tdat = np.array(tbot4harm.isel(dayofyear=iday, lon=ilon))
            dtdy[iday,5:nlat-5,ilon] = (1./a)* \
             ( (tdat[10:nlat] - tdat[0:nlat-10]) )/ \
             (latrad[10:nlat] - latrad[0:nlat-10])

    # calculate advection across climatological gradients
    ystart = 1979 ; yend = 2014 ; nyears=yend-ystart+1
    for iyear in np.arange(ystart,yend+1):
        print(iyear)
        ubotf = path+iexp+"/day/UBOT/UBOT_"+str(iyear)+".nc"
        vbotf = path+iexp+"/day/VBOT/VBOT_"+str(iyear)+".nc"
        ubot = read.read_sfc_cesm(ubotf,str(iyear)+"-01-01",str(iyear)+"-12-31")
        vbot = read.read_sfc_cesm(vbotf,str(iyear)+"-01-01",str(iyear)+"-12-31")

        advection = -1.*np.array(ubot.UBOT)*dtdx - 1.*np.array(vbot.VBOT)*dtdy

        advection_xr = xr.DataArray(advection, coords = [ubot.time, ubot.lat, ubot.lon], 
                       dims=['time','lat','lon'], name='advection')
        advection_xr.to_netcdf(path=outdir+'advection_'+str(iyear)+'.nc')















