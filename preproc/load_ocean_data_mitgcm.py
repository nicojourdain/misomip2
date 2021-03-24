# N. Jourdain, IGE-CNRS, MAR-2021

import numpy as np
import xarray as xr
import gsw
from pyproj import Proj
from .def_grids import grid_bounds

#====================================================================================
def load_oce_mod_mitgcm(files_in='MITgcm_output.nc',\
                        rho0=1026.0, teos10=False, region='Amundsen' ):
   """ Read MITgcm outputs on stereographic grid and define an xarray 
       dataset containing all variables required in MISOMIP2

       (XC,YC) -> gridT
       (XG,YC) -> gridU
       (XC,YG) -> gridV     

   """

   RT = 6.371e6 # Earth radius in meters

   ncTUV = xr.open_mfdataset(files_in,decode_coords=False)

   print(ncTUV)

   mtime = ncTUV.time.shape[0]

   # longitude & latitude on U, V, T grids
   if ( ( ncTUV.XC.min() < -180.1 ) | ( ncTUV.XC.max() > 360.1 ) ):
      print('!!! Assuming that (XC,YC) are stereographic coordinates (EPSG:3031) !!!')
      p = Proj('+init=EPSG:3031')
      XC2d, YC2d = np.meshgrid( ncTUV.XC.values, ncTUV.YC.values )
      XG2d, YG2d = np.meshgrid( ncTUV.XG.values, ncTUV.YG.values )
      lons, lats = p(XC2d, YC2d, inverse=True)
      lonT = xr.DataArray( lons, dims=['YC', 'XC'] )
      latT = xr.DataArray( lats, dims=['YC', 'XC'] )
      lons, lats = p(XG2d, YC2d, inverse=True)
      lonU = xr.DataArray( lons, dims=['YC', 'XG'] )
      latU = xr.DataArray( lats, dims=['YC', 'XG'] )
      lons, lats = p(XC2d, YG2d, inverse=True)
      lonV = xr.DataArray( lons, dims=['YG', 'XC'] )
      latV = xr.DataArray( lats, dims=['YG', 'XC'] )
   else:
      print('!!! Assuming that (XC,YC) are (longitude,latitude) !!!')
      lonT = ncTUV.XC
      latT = ncTUV.YC
      lonU = ncTUV.XG
      latU = ncTUV.XC
      lonV = ncTUV.XC
      latV = ncTUV.YG

   # save original domain boundaries:
   domain_minlat = latT.min().values
   domain_maxlat = latT.max().values
   domain_minlon = lonT.min().values
   domain_maxlon = lonT.max().values

   # grid mesh widths along x and y at C/T, U and V points [m]:
   dxT = xr.DataArray( 0.500000000 * (ncTUV.dxC.values+ncTUV.dxC.shift(XG=-1).values), dims=['YC', 'XC'] )
   dyT = xr.DataArray( 0.500000000 * (ncTUV.dyC.values+ncTUV.dyC.shift(YG=-1).values), dims=['YC', 'XC'] )
   dxU = ncTUV.dxC
   dxV = ncTUV.dxG

   dlatTdy = 0.500000000 * ( latT.shift(XC=1) - latT.shift(XC=-1) )
   dlatUdy = 0.500000000 * ( latU.shift(XG=1) - latU.shift(XG=-1) )
   dlatVdy = 0.500000000 * ( latV.shift(XC=1) - latV.shift(XC=-1) )

   # depth of U, V, C/T grids (neglecting the effects of partial steps in the interpolation) [m, positive in the ocean]
   depTUV=ncTUV.Z*(-1) 

   # local C/T, U, V grid rotation angle compared to the (zonal,meridional) direction [rad]
   thetaT = np.arcsin( RT*dlatTdy*np.pi/180. / dxT  )
   thetaU = np.arcsin( RT*dlatUdy*np.pi/180. / dxU  )
   thetaV = np.arcsin( RT*dlatVdy*np.pi/180. / dxV  )
   print('Minimum local grid angle in degrees w.r.t. (zonal,meridional):',thetaU.min().values*180./np.pi)
   print('Maximum local grid angle in degrees w.r.t. (zonal,meridional):',thetaU.max().values*180./np.pi)

   # Define useful masks on U, V, C/T grids 
   maskT = ncTUV.hFacC # 3d mask (=1 if ocean, =0 elsewhere)
   maskU = ncTUV.hFacW # 3d mask
   maskV = ncTUV.hFacS # 3d mask
   
   [mz, my, mx] = maskT.shape

   # Domain mask (ones with a halo of nans), used not to interpolate beyond:
   halonan = np.ones((my,mx))
   halonan[0,:] = np.nan ; halonan[-1,:] = np.nan
   halonan[:,0] = np.nan ; halonan[:,-1] = np.nan
   DOMMSKT = xr.DataArray( halonan, dims=['YC', 'XC'] )
   DOMMSKU = xr.DataArray( halonan, dims=['YC', 'XG'] )
   DOMMSKV = xr.DataArray( halonan, dims=['YG', 'XC'] )

   # Ocean fraction at each level:
   LEVOF = maskT*100.0

   # ice-shelf fraction:
   SFTFLI = LEVOF.max('Z').where( (LEVOF.max('Z')>1.0) & (LEVOF.isel(Z=0)==0.0), 0.0 )

   # Bathymetry (including under ice shelves) [m, positive]
   if ( "Depth" in ncTUV.data_vars ):
     DEPTHO = ncTUV.Depth
   else:
     print('@@@@@ WARNING @@@@@   No data found for DEPTHO  -->  filled with NaNs')
     DEPTHO = xr.DataArray( np.zeros((my,mx))*np.nan, dims=['YC', 'XC'] )
  
   # Depth of ice shelf draft [m]:
   if ( "isfdraft" in ncTUV.data_vars ):
     DEPFLI = ncTUV.isfdraft
   elif ( "ice_draft" in ncTUV.data_vars ):
     DEPFLI = ncTUV.ice_draft
   elif ( "ice_shelf_draft" in ncTUV.data_vars ):
     DEPFLI = ncTUV.ice_shelf_draft
   else:
     [IndIS_y,IndIS_x] = np.where( ( ncTUV.hFacC.isel(Z=0).values == 0. ) & ( ncTUV.hFacC.sum('Z').values > 1. ) )
     dz = ncTUV.Zl.values - ncTUV.Zu.values
     ISdraft = np.zeros((my,mx))
     for i in range(np.size(IndIS_x)):
       Ind_firstnonempty = np.where((ncTUV.hFacC.values[:,IndIS_y[i],IndIS_x[i]])>0)[0][0]
       ISdraft[IndIS_y[i],IndIS_x[i]] = ncTUV.Zl.values[Ind_firstnonempty]-ncTUV.hFacC.values[Ind_firstnonempty,IndIS_y[i],IndIS_x[i]]*dz[Ind_firstnonempty]
     DEPFLI = xr.DataArray( ISdraft, dims=['YC', 'XC'] )
 
   # ocean temperature [degC]
   if ( "toce" in ncTUV.data_vars ):
     TT = ncTUV.toce
   elif ( "thetao" in ncTUV.data_vars ):
     TT = ncTUV.thetao
   elif ( "THETA" in ncTUV.data_vars ):
     TT = ncTUV.THETA
   else:
     print('@@@@@ WARNING @@@@@   No data found for TT  -->  filled with NaNs')
     TT = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'Z', 'YC', 'XC'] )

   # ocean salinity [1.e-3]
   if ( "soce" in ncTUV.data_vars ):
     SS = ncTUV.soce
   elif ( "so" in ncTUV.data_vars ):
     SS = ncTUV.so
   elif ( "SALT" in ncTUV.data_vars ):
     SS = ncTUV.SALT
   else:
     print('@@@@@ WARNING @@@@@   No data found for SS  -->  filled with NaNs')
     SS = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'Z', 'YC', 'XC'] )

   # sea bottom ocean temperature [degC]
   if ( "sbt" in ncTUV.data_vars ):
     TTB = ncTUV.sbt
   elif ( "sosbt" in ncTUV.data_vars ):
     TTB = ncTUV.sosbt
   elif ( "tob" in ncTUV.data_vars ):
     TTB = ncTUV.tob
   else:
     print('@@@@@ WARNING @@@@@   No data found for TTB  -->  filled with NaNs')
     TTB = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )

   # sea bottom ocean salinity [1.e-3]
   if ( "sbs" in ncTUV.data_vars ):
     SSB = ncTUV.sbs
   elif ( "sosbs" in ncTUV.data_vars ):
     SSB = ncTUV.sosbs
   elif ( "sob" in ncTUV.data_vars ):
     SSB = ncTUV.sob
   else:
     print('@@@@@ WARNING @@@@@   No data found for SSB  -->  filled with NaNs')
     SSB = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )

   # ocean potential temperature and practical salinity :
   if (teos10):
     TOB = xr.apply_ufunc(gsw.pt_from_CT, SSB, TTB)
     SOB = xr.apply_ufunc(gsw.SP_from_SA, SSB, DEPTHO, lonT, latT)
     THETAO = xr.apply_ufunc(gsw.pt_from_CT, SS, TT)
     SO = xr.apply_ufunc(gsw.SP_from_SA, SS, DEPTHO, lonT, latT)
   else: 
     TOB = TTB
     SOB = SSB
     THETAO = TT
     SO = SS
   
   # ocean x-ward velocity [m s-1]
   if ( "uoce" in ncTUV.data_vars ):
     UX = ncTUV.uoce
   elif ( "UVEL" in ncTUV.data_vars ):
     UX = ncTUV.UVEL
   elif ( "uo" in ncTUV.data_vars ):
     UX = ncTUV.uo
   else:
     print('@@@@@ WARNING @@@@@   No data found for UX  -->  filled with NaNs')
     UX = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'Z', 'YC', 'XG'] )

   # ocean y-ward velocity [m s-1]
   if ( "voce" in ncTUV.data_vars ):
     VY = ncTUV.voce
   elif ( "VVEL" in ncTUV.data_vars ):
     VY = ncTUV.VVEL
   elif ( "vo" in ncTUV.data_vars ):
     VY = ncTUV.vo
   else:
     print('@@@@@ WARNING @@@@@   No data found for VY  -->  filled with NaNs')
     VY = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'Z', 'YG', 'XC'] )

   # surface stress received by the ocean along x [W m-1]
   if ( "utau" in ncTUV.data_vars ):
     TAUX = ncTUV.utau
   elif ( "TAUX" in ncTUV.data_vars ):
     TAUX = ncTUV.TAUX
   else:
     print('@@@@@ WARNING @@@@@   No data found for TAUX  -->  filled with NaNs')
     TAUX = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XG'] )

   # surface stress received by the ocean along x [W m-1]
   if ( "vtau" in ncTUV.data_vars ):
     TAUY = ncTUV.vtau
   elif ( "TAUY" in ncTUV.data_vars ):
     TAUY = ncTUV.TAUY
   else:
     print('@@@@@ WARNING @@@@@   No data found for TAUY  -->  filled with NaNs')
     TAUY = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YG', 'XC'] )

   # Sea surface height [m]
   if ( "zos" in ncTUV.data_vars ):
     ZOS = ncTUV.zos
   elif ( "SSH" in ncTUV.data_vars ):
     ZOS = ncTUV.SSH
   elif ( "ssh" in ncTUV.data_vars ):
     ZOS = ncTUV.ssh
   elif ( "ETAN" in ncTUV.data_vars ):
     ZOS = ncTUV.ETAN
   else:
     print('@@@@@ WARNING @@@@@   No data found for ZOS  -->  filled with NaNs')
     ZOS = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )

   # mass barotropic streamfunction
   # see Griffies et al. (2016, section H26): d(psi)/dy=-U (U: x-ward mass transport), d(psi)/dx=V (V: yward mass transport)
   if ( "sobarstf" in ncTUV.data_vars ):
     MSFTBAROT = ncTUV.sobarstf * rho0
   else:
     print('@@@@@ WARNING @@@@@   No data found for MSFTBAROT  -->  filled with NaNs')
     MSFTBAROT = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )

   # ice shelf melt [kg m-2 s-1, positive for actual melting] :
   if ( "fwfisf" in ncTUV.data_vars ):
     FICESHELF = ncTUV.fwfisf*(-1)
   elif ( "SHIfwFlx" in ncTUV.data_vars ):
     FICESHELF = ncTUV.SHIfwFlx*(-1)
   else:
     print('@@@@@ WARNING @@@@@   No data found for FICESHELF  -->  filled with NaNs')
     FICESHELF = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )

   # sea-ice concentration [0-100]
   if ( "siconc" in ncTUV.data_vars ):
     SICONC = ncTUV.siconc*100.0
     SICONC = SICONC.where( (~np.isnan(SICONC.values)) & (~np.isinf(SICONC.values)), 0.e0 )
   else:
     print('@@@@@ WARNING @@@@@   No data found for SICONC  -->  filled with NaNs')
     SICONC = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )   

   # sea-ice volume per area [m]
   if ( "sivolu" in ncTUV.data_vars ):
     SIVOL = ncTUV.sivolu
   elif ( "sivol" in ncTUV.data_vars ):
     SIVOL = ncTUV.sivol
   else:
     print('@@@@@ WARNING @@@@@   No data found for SIVOL  -->  filled with NaNs')
     SIVOL = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )

   # sea-ice x-ward velocity [m/s]
   if ( "sivelu" in ncTUV.data_vars ):
     SIUX = ncTUV.sivelu
   elif ("siu" in ncTUV.data_vars ):
     SIUX = ncTUV.siu
   else:
     print('@@@@@ WARNING @@@@@   No data found for SIUX  -->  filled with NaNs')
     SIUX = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )

   # sea-ice y-ward velocity [m/s]
   if ( "sivelv" in ncTUV.data_vars ):
     SIVY = ncTUV.sivelv
   elif ("siv" in ncTUV.data_vars ):
     SIVY = ncTUV.siv
   else:
     print('@@@@@ WARNING @@@@@   No data found for SIUY  -->  filled with NaNs')
     SIVY = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )

   # Total heat flux received by the ocean surface (including ice-shelf/ocean interface) [W m-2] 
   # see Griffies et al. (2016, section K4-K5) NB: here, including correction if any unlike Griffies (to avoid 2 variables)
   if ( "qt_oce" in ncTUV.data_vars ):
     HFDS = ncTUV.qt_oce
   else:
     print('@@@@@ WARNING @@@@@   No data found for HFDS  -->  filled with NaNs')
     HFDS = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )

   # Water flux entering the ocean due to sea-ice (melting-freezing) and surface correction (SSS restoring)
   # (= fsitherm + wfocorr in Griffies 2016 section K2) [kg m-2 s-1]
   if ( "wfocorr" in ncTUV.data_vars ):
     WFOCORR = - ncTUV.wfocorr
   else:
     WFOCORR = xr.DataArray( np.zeros((mtime,my,mx)), dims=['time', 'YC', 'XC'] )
   if ( "fsitherm" in ncTUV.data_vars ):
     WFOSICOR = WFOCORR - ncTUV.fsitherm
   else:
     print('@@@@@ WARNING @@@@@   No data found for WFOSICOR  -->  filled with NaNs')
     WFOSICOR = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )

   # Water flux entering the ocean due to rainfall, snowfall, condensation - evap, 
   # river runoff, iceberg and ice-shelf melt [kg m-2 s-1]  
   # (= pr+prs+evs+ficeberg+friver+ficeshelf in Griffies 2016, section K2)
   if ( "empmr" in ncTUV.data_vars ):
     WFOATRLI = - ncTUV.empmr + FICESHELF
   else:
     print('@@@@@ WARNING @@@@@   No data found for WFOATRLI  -->  filled with NaNs')
     WFOATRLI = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'YC', 'XC'] )
  
   #----------
   # Reduce the size of ocean dataset
 
   [lonmin,lonmax,latmin,latmax] = grid_bounds(region=region)
   lonmin=lonmin-1.1 # take a bit more for interpolation
   lonmax=lonmax+1.1
   latmin=latmin-1.1
   latmax=latmax+1.1

   condT2d = ( (latT >= latmin) & (latT <= latmax) & (lonT >= lonmin) & (lonT <= lonmax) )

   for ii in np.arange(latT.shape[1]):
      if ( np.sum(condT2d.isel(XC=ii).values) == 0 ):
        imin=ii
      else:
        imin=ii
        break
   for ii in np.arange(latT.shape[1]-1,0,-1):
      if ( np.sum(condT2d.isel(XC=ii).values) == 0 ):
        imax=ii
      else:
        imax=ii
        break
   for jj in np.arange(latT.shape[0]):
      if ( np.sum(condT2d.isel(YC=jj).values) == 0 ):
        jmin=jj
      else:
        jmin=jj
        break
   for jj in np.arange(latT.shape[0]-1,0,-1):
      if ( np.sum(condT2d.isel(YC=jj).values) == 0 ):
        jmax=jj
      else:
        jmax=jj
        break

   print([imin,imax,jmin,jmax])

   #----------
   # Create new xarray dataset including all useful variables:
   # reshaping (x,y) as 1-dimensional (sxy)

   nxy=(jmax-jmin+1)*(imax-imin+1)

   ds = xr.Dataset(
      {
       "SO":        (["time", "z", "sxy"], np.reshape( SO.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,mz,nxy)) ),
       "THETAO":    (["time", "z", "sxy"], np.reshape( THETAO.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,mz,nxy)) ),
       "UX":        (["time", "z", "sxy"], np.reshape( UX.isel(XG=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,mz,nxy)) ),
       "VY":        (["time", "z", "sxy"], np.reshape( VY.isel(XC=slice(imin,imax+1),YG=slice(jmin,jmax+1)).values, (mtime,mz,nxy)) ),
       "TAUX":      (["time", "sxy"], np.reshape( TAUX.isel(XG=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "TAUY":      (["time", "sxy"], np.reshape( TAUY.isel(XC=slice(imin,imax+1),YG=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "ZOS":       (["time", "sxy"], np.reshape( ZOS.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "TOB":       (["time", "sxy"], np.reshape( TOB.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SOB":       (["time", "sxy"], np.reshape( SOB.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "FICESHELF": (["time", "sxy"], np.reshape( FICESHELF.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "MSFTBAROT": (["time", "sxy"], np.reshape( MSFTBAROT.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "HFDS":      (["time", "sxy"], np.reshape( HFDS.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "WFOATRLI":  (["time", "sxy"], np.reshape( WFOATRLI.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "WFOSICOR":  (["time", "sxy"], np.reshape( WFOSICOR.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SICONC":    (["time", "sxy"], np.reshape( SICONC.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIVOL":     (["time", "sxy"], np.reshape( SIVOL.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIUX":      (["time", "sxy"], np.reshape( SIUX.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIVY":      (["time", "sxy"], np.reshape( SIVY.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "LEVOF":     (["z", "sxy"], np.reshape( LEVOF.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mz,nxy)) ),
       "maskT":     (["z", "sxy"], np.reshape( maskT.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mz,nxy)) ),
       "maskU":     (["z", "sxy"], np.reshape( maskU.isel(XG=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, (mz,nxy)) ),
       "maskV":     (["z", "sxy"], np.reshape( maskV.isel(XC=slice(imin,imax+1),YG=slice(jmin,jmax+1)).values, (mz,nxy)) ),
       "SFTFLI":    (["sxy"], np.reshape( SFTFLI.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "DEPFLI":    (["sxy"], np.reshape( DEPFLI.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "DEPTHO":    (["sxy"], np.reshape( DEPTHO.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKT":   (["sxy"], np.reshape( DOMMSKT.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKU":   (["sxy"], np.reshape( DOMMSKU.isel(XG=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKV":   (["sxy"], np.reshape( DOMMSKV.isel(XC=slice(imin,imax+1),YG=slice(jmin,jmax+1)).values, nxy) ),
       "lonT":      (["sxy"], np.reshape( lonT.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "lonU":      (["sxy"], np.reshape( lonU.isel(XG=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "lonV":      (["sxy"], np.reshape( lonV.isel(XC=slice(imin,imax+1),YG=slice(jmin,jmax+1)).values, nxy) ),
       "latT":      (["sxy"], np.reshape( latT.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "latU":      (["sxy"], np.reshape( latU.isel(XG=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "latV":      (["sxy"], np.reshape( latV.isel(XC=slice(imin,imax+1),YG=slice(jmin,jmax+1)).values, nxy) ),
       "dxT":       (["sxy"], np.reshape( dxT.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "dyT":       (["sxy"], np.reshape( dyT.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "thetaT":    (["sxy"], np.reshape( thetaT.isel(XC=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "thetaU":    (["sxy"], np.reshape( thetaU.isel(XG=slice(imin,imax+1),YC=slice(jmin,jmax+1)).values, nxy) ),
       "thetaV":    (["sxy"], np.reshape( thetaV.isel(XC=slice(imin,imax+1),YG=slice(jmin,jmax+1)).values, nxy) ),
       "depTUV":    (['z'], depTUV.values)
      },
      coords={
      "time": ncTUV.time.values
      },
      attrs={
      "original_minlat": domain_minlat,
      "original_maxlat": domain_maxlat,
      "original_minlon": domain_minlon,
      "original_maxlon": domain_maxlon
      },
   )

   return ds

 
