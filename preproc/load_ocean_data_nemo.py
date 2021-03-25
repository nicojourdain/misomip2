# N. Jourdain, IGE-CNRS, MAR-2021

import numpy as np
import xarray as xr
import gsw
from .def_grids import grid_bounds
from .def_attrs import rename_dimensions

#====================================================================================
def load_oce_mod_nemo(file_mesh_mask='mesh_mask.nc',\
                      file_bathy='bathy_meter.nc',\
                      files_gridT='nemo_grid_T.nc',\
                      files_gridU='nemo_grid_U.nc',\
                      files_gridV='nemo_grid_V.nc',\
                      files_SBC='nemo_flxT.nc',\
                      files_ice='nemo_icemod.nc',\
                      files_BSF='nemo_psi.nc',\
                      rho0=1026.0, teos10=False, region='Amundsen' ):
   """ Read NEMO outputs and define an xarray dataset containing 
       all variables required in MISOMIP2

       Adapted to NEMO's C-grid and standard mesh/mask variables.

       rho0 corresponds to rau0 value in NEMO's eosbn2.F90 i.e. reference volumic mass [kg m-3]

       teos10=False -> assumes the nemo outputs are in potential temperature & practical salinity (EOS80)
             =True  -> assumes the nemo outputs are in CT and AS and convert to PT and PS

       region = 'Amundsen' (default) or 'Weddell'

       NB: files_BSF contains the barotropic streamfunction calculated at U-points, e.g. using
           the cdfpsi function which is part of the cdftools (https://github.com/meom-group/CDFTOOLS) 

       Example1:
       ds = load_oce_mod_nemo()

       Example2:
       dir= 'datadir/model/'
       fT = [ dir+'nemo_y1990_grid_T.nc', dir+'nemo_y1991_grid_T.nc' ]
       fU = [ dir+'nemo_y1990_grid_U.nc', dir+'nemo_y1991_grid_U.nc' ]
       fV = [ dir+'nemo_y1990_grid_V.nc', dir+'nemo_y1991_grid_V.nc' ]
       fS = [ dir+'nemo_y1990_SBC.nc', dir+'nemo_y1991_SBC.nc' ]
       fI = [ dir+'nemo_y1990_icemod.nc', dir+'nemo_y1991_icemod.nc' ]
       fP = [ dir+'nemo_y1990_psi.nc', dir+'nemo_y1991_psi.nc' ]
       ds = load_oce_mod_nemo(files_gridT=fT,files_gridU=fU,files_gridV=fV,files_SBC=fS,files_ice=fI,files_BSF=fP,rho0=1028.0,teos10=True,region='Weddell')

   """

   RT = 6.371e6 # Earth radius in meters

   ncM = xr.open_dataset(file_mesh_mask,decode_coords=False) ; ncM=rename_dimensions(ncM); ncM=ncM.isel(time=0)
   ncB = xr.open_dataset(file_bathy,decode_coords=False) ; ncB=rename_dimensions(ncB)
   ncT = xr.open_mfdataset(files_gridT,decode_coords=False); ncT=rename_dimensions(ncT)
   ncU = xr.open_mfdataset(files_gridU,decode_coords=False); ncU=rename_dimensions(ncU)
   ncV = xr.open_mfdataset(files_gridV,decode_coords=False); ncV=rename_dimensions(ncV)
   ncS = xr.open_mfdataset(files_SBC,decode_coords=False); ncS=rename_dimensions(ncS)
   ncI = xr.open_mfdataset(files_ice,decode_coords=False); ncI=rename_dimensions(ncI)
   ncP = xr.open_mfdataset(files_BSF,decode_coords=False); ncP=rename_dimensions(ncP)

   mtime = ncT.time.shape[0]

   # Define useful masks on U, V, T grids 
   maskT = ncM.tmask # 3d mask (=1 if ocean, =0 elsewhere)
   maskU = ncM.umask # 3d mask
   maskV = ncM.vmask # 3d mask

   [mz, my, mx] = maskT.shape

   # Domain mask (ones with a halo of nans), used not to interpolate beyond:
   halonan = np.ones((my,mx))
   halonan[0,:] = np.nan ; halonan[-1,:] = np.nan
   halonan[:,0] = np.nan ; halonan[:,-1] = np.nan
   DOMMSKT = xr.DataArray( halonan, dims=['y', 'x'] )
   DOMMSKU = xr.DataArray( halonan, dims=['y', 'x'] )
   DOMMSKV = xr.DataArray( halonan, dims=['y', 'x'] )

   # longitude & latitude on U, V, T grids
   lonT = ncM.glamt ; latT = ncM.gphit
   lonU = ncM.glamu ; latU = ncM.gphiu
   lonV = ncM.glamv ; latV = ncM.gphiv

   # save original domain boundaries:
   domain_minlat = latT.min().values
   domain_maxlat = latT.max().values
   domain_minlon = lonT.min().values
   domain_maxlon = lonT.max().values

   dxT = ncM.e1t # grid mesh width along x-axis in meters for sea ice grid (T)
   dyT = ncM.e2t # grid mesh width along y-axis in meters for sea ice grid (T)
   dxU = ncM.e1u # grid mesh width along x-axis in meters for UX
   dxV = ncM.e1v # grid mesh width along x-axis in meters for VY
   dlatTdy = 0.500000000 * ( latT.shift(x=1) - latT.shift(x=-1) )
   dlatUdy = 0.500000000 * ( latU.shift(x=1) - latU.shift(x=-1) )
   dlatVdy = 0.500000000 * ( latV.shift(x=1) - latV.shift(x=-1) )

   # depth of U, V, T grids (neglecting the effects of partial steps in the interpolation)
   depTUV=ncM.gdept_1d

   # local T, U, V grid rotation angle compared to the (zonal,meridional) direction [rad]
   thetaT = np.arcsin( RT*dlatTdy*np.pi/180. / dxT  )
   thetaU = np.arcsin( RT*dlatUdy*np.pi/180. / dxU  )
   thetaV = np.arcsin( RT*dlatVdy*np.pi/180. / dxV  )
   print('Minimum local grid angle in degrees w.r.t. (zonal,meridional):',thetaU.min().values*180./np.pi)
   print('Maximum local grid angle in degrees w.r.t. (zonal,meridional):',thetaU.max().values*180./np.pi)

   # Ocean fraction at each level:
   LEVOF = maskT*100.0

   # 2d ice-shelf fractoin:
   SFTFLI = ncM.misf*1.e0
   SFTFLI = SFTFLI.where( (SFTFLI.values > 1.5), 0.e0 )
   SFTFLI = SFTFLI.where( (SFTFLI.values < 1.5), 100. )

   # Bathymetry (including under ice shelves) [m]
   # (if possible after NEMO's initialization, i.e. from mesh_mask) :
   if ( "bathy_meter" in ncM.data_vars ):
     DEPTHO = ncM.bathy_meter
   elif ( "bathy_metry" in ncM.data_vars ):
     DEPTHO = ncM.bathy_metry
   elif ( "Bathymetry_isf" in ncB.data_vars ):
     DEPTHO = ncB.Bathymetry_isf
   elif ( "Bathymetry" in ncB.data_vars ):
     DEPTHO = ncB.Bathymetry
   elif ( "bathy_meter" in ncB.data_vars ):
     DEPTHO = ncB.bathy_meter
   else:
     print('@@@@@ WARNING @@@@@   No data found for DEPTHO  -->  filled with NaNs')
     DEPTHO = xr.DataArray( np.zeros((my,mx))*np.nan, dims=['y', 'x'] )
  
   # Depth of ice shelf draft (if possible after NEMO's initialization, i.e. from mesh_mask) [m]:
   if ( "isfdraft" in ncM.data_vars ):
     DEPFLI = ncM.isfdraft
   elif ( "isfdraft" in ncB.data_vars ):
     DEPFLI = ncB.isfdraft
   else:
     print('@@@@@ WARNING @@@@@   No data found for DEPFLI  -->  filled with NaNs')
     DEPFLI = xr.DataArray( np.zeros((my,mx))*np.nan, dims=['y', 'x'] )
 
   # ocean temperature [degC]
   if ( "toce" in ncT.data_vars ):
     TT = ncT.toce
   elif ( "thetao" in ncT.data_vars ):
     TT = ncT.thetao
   elif ( "votemper" in ncT.data_vars ):
     TT = ncT.votemper
   else:
     print('@@@@@ WARNING @@@@@   No data found for TT  -->  filled with NaNs')
     TT = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'z', 'y', 'x'] )

   # ocean salinity [1.e-3]
   if ( "soce" in ncT.data_vars ):
     SS = ncT.soce
   elif ( "so" in ncT.data_vars ):
     SS = ncT.so
   elif ( "vosaline" in ncT.data_vars ):
     SS = ncT.vosaline
   else:
     print('@@@@@ WARNING @@@@@   No data found for SS  -->  filled with NaNs')
     SS = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'z', 'y', 'x'] )

   # sea bottom ocean temperature [degC]
   if ( "sbt" in ncT.data_vars ):
     TTB = ncT.sbt
   elif ( "sosbt" in ncT.data_vars ):
     TTB = ncT.sosbt
   elif ( "tob" in ncT.data_vars ):
     TTB = ncT.tob
   else:
     print('@@@@@ WARNING @@@@@   No data found for TTB  -->  filled with NaNs')
     TTB = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # sea bottom ocean salinity [1.e-3]
   if ( "sbs" in ncT.data_vars ):
     SSB = ncT.sbs
   elif ( "sosbs" in ncT.data_vars ):
     SSB = ncT.sosbs
   elif ( "sob" in ncT.data_vars ):
     SSB = ncT.sob
   else:
     print('@@@@@ WARNING @@@@@   No data found for SSB  -->  filled with NaNs')
     SSB = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

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
   if ( "uoce" in ncU.data_vars ):
     UX = ncU.uoce
   elif ( "vozocrtx" in ncU.data_vars ):
     UX = ncU.vozocrtx
   elif ( "uo" in ncU.data_vars ):
     UX = ncU.uo
   else:
     print('@@@@@ WARNING @@@@@   No data found for UX  -->  filled with NaNs')
     UX = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'z', 'y', 'x'] )

   # ocean y-ward velocity [m s-1]
   if ( "voce" in ncV.data_vars ):
     VY = ncV.voce
   elif ( "vomecrty" in ncV.data_vars ):
     VY = ncV.vomecrty
   elif ( "vo" in ncV.data_vars ):
     VY = ncV.vo
   else:
     print('@@@@@ WARNING @@@@@   No data found for VY  -->  filled with NaNs')
     VY = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'z', 'y', 'x'] )

   # surface stress received by the ocean along x [W m-1]
   if ( "utau" in ncU.data_vars ):
     TAUX = ncU.utau
   elif ( "sozotaux" in ncU.data_vars ):
     TAUX = ncU.sozotaux
   else:
     print('@@@@@ WARNING @@@@@   No data found for TAUX  -->  filled with NaNs')
     TAUX = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # surface stress received by the ocean along x [W m-1]
   if ( "vtau" in ncV.data_vars ):
     TAUY = ncV.vtau
   elif ( "sometauy" in ncV.data_vars ):
     TAUY = ncV.sometauy
   else:
     print('@@@@@ WARNING @@@@@   No data found for TAUY  -->  filled with NaNs')
     TAUY = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # Sea surface height [m]
   if ( "zos" in ncS.data_vars ):
     ZOS = ncS.zos
   elif ( "SSH" in ncS.data_vars ):
     ZOS = ncS.SSH
   elif ( "ssh" in ncS.data_vars ):
     ZOS = ncS.ssh
   elif ( "sossheig" in ncS.data_vars ):
     ZOS = ncS.sossheig
   elif ( "zos" in ncT.data_vars ):
     ZOS = ncT.zos
   elif ( "SSH" in ncT.data_vars ):
     ZOS = ncT.SSH
   elif ( "ssh" in ncT.data_vars ):
     ZOS = ncT.ssh
   elif ( "sossheig" in ncT.data_vars ):
     ZOS = ncT.sossheig
   else:
     print('@@@@@ WARNING @@@@@   No data found for ZOS  -->  filled with NaNs')
     ZOS = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # mass barotropic streamfunction
   # see Griffies et al. (2016, section H26): d(psi)/dy=-U (U: x-ward mass transport), d(psi)/dx=V (V: yward mass transport)
   if ( "sobarstf" in ncP.data_vars ):
     MSFTBAROT = ncP.sobarstf * rho0
   else:
     print('@@@@@ WARNING @@@@@   No data found for MSFTBAROT  -->  filled with NaNs')
     MSFTBAROT = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # ice shelf melt [kg m-2 s-1, positive for actual melting] :
   if ( "fwfisf" in ncS.data_vars ):
     FICESHELF = ncS.fwfisf*(-1)
   elif ( "sowflisf_cav" in ncS.data_vars ):
     FICESHELF = ncS.sowflisf_cav*(-1)
   else:
     print('@@@@@ WARNING @@@@@   No data found for FICESHELF  -->  filled with NaNs')
     FICESHELF = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # sea-ice concentration [0-100]
   if ( "siconc" in ncI.data_vars ):
     SICONC = ncI.siconc*100.0
     SICONC = SICONC.where( (~np.isnan(SICONC.values)) & (~np.isinf(SICONC.values)), 0.e0 )
   else:
     print('@@@@@ WARNING @@@@@   No data found for SICONC  -->  filled with NaNs')
     SICONC = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )   

   # sea-ice volume per area [m]
   if ( "sivolu" in ncI.data_vars ):
     SIVOL = ncI.sivolu
   elif ( "sivol" in ncI.data_vars ):
     SIVOL = ncI.sivol
   else:
     print('@@@@@ WARNING @@@@@   No data found for SIVOL  -->  filled with NaNs')
     SIVOL = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # sea-ice x-ward velocity [m/s]
   if ( "sivelu" in ncI.data_vars ):
     SIUX = ncI.sivelu
   elif ("siu" in ncI.data_vars ):
     SIUX = ncI.siu
   else:
     print('@@@@@ WARNING @@@@@   No data found for SIUX  -->  filled with NaNs')
     SIUX = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # sea-ice y-ward velocity [m/s]
   if ( "sivelv" in ncI.data_vars ):
     SIVY = ncI.sivelv
   elif ("siv" in ncI.data_vars ):
     SIVY = ncI.siv
   else:
     print('@@@@@ WARNING @@@@@   No data found for SIUY  -->  filled with NaNs')
     SIVY = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # Total heat flux received by the ocean surface (including ice-shelf/ocean interface) [W m-2] 
   # see Griffies et al. (2016, section K4-K5) NB: here, including correction if any unlike Griffies (to avoid 2 variables)
   if ( "qt_oce" in ncS.data_vars ):
     HFDS = ncS.qt_oce
   elif ( "sohefldo" in ncS.data_vars ):
     HFDS = ncS.sohefldo
   if ( "qisf" in ncS.data_vars ):
     HFDS = HFDS + ncS.qisf # not included in qt_oce in tested NEMO versions
   elif ( "qoceisf_cav" in ncS.data_vars ):
     HFDS = HFDS + ncS.qoceisf_cav # not included in sohefldo in tested NEMO versions

   # Water flux entering the ocean due to sea-ice (melting-freezing) and surface correction (SSS restoring)
   # (= fsitherm + wfocorr in Griffies 2016 section K2) [kg m-2 s-1]
   if ( "erp" in ncS.data_vars ):
     ERP = ncS.erp.where( (~np.isnan(ncS.erp.values)), 0.e0 ) # surface correction (SSS restoring)
   elif ( "sowafld" in ncS.data_vars ):
     ERP = ncS.sowafld.where( (~np.isnan(ncS.sowafld.values)), 0.e0 ) # surface correction (SSS restoring)
   else:
     print('@@@@@ WARNING @@@@@   No data found for ERP  -->  filled with NaNs')
     ERP = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   if ( "saltflx" in ncS.data_vars ): # NB: saltflx unit attribute is wrong in nico's output, it is actually in [1e-3 kg m-2 s-1]
     WFOSICOR = - ERP - ncS.saltflx / SS.isel(z=0)
   elif ( "sosfldow" in ncS.data_vars ):
     WFOSICOR = - ncS.sosfldow / SS.isel(z=0) # includes surface correction (ERP) in Pierre's version
   elif ( "sfx" in ncI.data_vars ):
     WFOSICOR = - ERP - ncI.sfx/86400.0 / SS.isel(z=0)
   else:
     print('@@@@@ WARNING @@@@@   No data found for WFOSICOR  -->  filled with NaNs')
     WFOSICOR = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # Water flux entering the ocean due to rainfall, snowfall, condensation - evap, 
   # river runoff, iceberg and ice-shelf melt [kg m-2 s-1]  
   # (= pr+prs+evs+ficeberg+friver+ficeshelf in Griffies 2016, section K2)
   if ( "empmr" in ncS.data_vars ):
     WFOATRLI = - ncS.empmr + FICESHELF
   elif ( "sowaflup" in ncS.data_vars ):
     WFOATRLI = - ncS.sowaflup - WFOSICOR + FICESHELF
   else:
     print('@@@@@ WARNING @@@@@   No data found for WFOATRLI  -->  filled with NaNs')
     WFOATRLI = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )
  
   #----------
   # Reduce the size of ocean dataset
 
   [lonmin,lonmax,latmin,latmax] = grid_bounds(region=region)
   lonmin=lonmin-1.1 # take a bit more for interpolation
   lonmax=lonmax+1.1
   latmin=latmin-1.1
   latmax=latmax+1.1

   condT2d = ( (latT >= latmin) & (latT <= latmax) & (lonT >= lonmin) & (lonT <= lonmax) )

   for ii in np.arange(latT.shape[1]):
      if ( np.sum(condT2d.isel(x=ii).values) == 0 ):
        imin=ii
      else:
        imin=ii
        break
   for ii in np.arange(latT.shape[1]-1,0,-1):
      if ( np.sum(condT2d.isel(x=ii).values) == 0 ):
        imax=ii
      else:
        imax=ii
        break
   for jj in np.arange(latT.shape[0]):
      if ( np.sum(condT2d.isel(y=jj).values) == 0 ):
        jmin=jj
      else:
        jmin=jj
        break
   for jj in np.arange(latT.shape[0]-1,0,-1):
      if ( np.sum(condT2d.isel(y=jj).values) == 0 ):
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
       "SO":        (["time", "z", "sxy"], np.reshape( SO.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,mz,nxy)) ),
       "THETAO":    (["time", "z", "sxy"], np.reshape( THETAO.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,mz,nxy)) ),
       "UX":        (["time", "z", "sxy"], np.reshape( UX.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,mz,nxy)) ),
       "VY":        (["time", "z", "sxy"], np.reshape( VY.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,mz,nxy)) ),
       "TAUX":      (["time", "sxy"], np.reshape( TAUX.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "TAUY":      (["time", "sxy"], np.reshape( TAUY.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "ZOS":       (["time", "sxy"], np.reshape( ZOS.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "TOB":       (["time", "sxy"], np.reshape( TOB.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SOB":       (["time", "sxy"], np.reshape( SOB.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "FICESHELF": (["time", "sxy"], np.reshape( FICESHELF.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "MSFTBAROT": (["time", "sxy"], np.reshape( MSFTBAROT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "HFDS":      (["time", "sxy"], np.reshape( HFDS.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "WFOATRLI":  (["time", "sxy"], np.reshape( WFOATRLI.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "WFOSICOR":  (["time", "sxy"], np.reshape( WFOSICOR.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SICONC":    (["time", "sxy"], np.reshape( SICONC.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIVOL":     (["time", "sxy"], np.reshape( SIVOL.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIUX":      (["time", "sxy"], np.reshape( SIUX.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIVY":      (["time", "sxy"], np.reshape( SIVY.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "LEVOF":     (["z", "sxy"], np.reshape( LEVOF.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mz,nxy)) ),
       "maskT":     (["z", "sxy"], np.reshape( maskT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mz,nxy)) ),
       "maskU":     (["z", "sxy"], np.reshape( maskU.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mz,nxy)) ),
       "maskV":     (["z", "sxy"], np.reshape( maskV.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mz,nxy)) ),
       "SFTFLI":    (["sxy"], np.reshape( SFTFLI.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "DEPFLI":    (["sxy"], np.reshape( DEPFLI.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "DEPTHO":    (["sxy"], np.reshape( DEPTHO.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKT":   (["sxy"], np.reshape( DOMMSKT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKU":   (["sxy"], np.reshape( DOMMSKU.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKV":   (["sxy"], np.reshape( DOMMSKV.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "lonT":      (["sxy"], np.reshape( lonT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "lonU":      (["sxy"], np.reshape( lonU.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "lonV":      (["sxy"], np.reshape( lonV.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "latT":      (["sxy"], np.reshape( latT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "latU":      (["sxy"], np.reshape( latU.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "latV":      (["sxy"], np.reshape( latV.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "dxT":       (["sxy"], np.reshape( dxT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "dyT":       (["sxy"], np.reshape( dyT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "thetaT":    (["sxy"], np.reshape( thetaT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "thetaU":    (["sxy"], np.reshape( thetaU.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "thetaV":    (["sxy"], np.reshape( thetaV.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "depTUV":    (['z'], depTUV.values)
      },
      coords={
      "time": ncT.time.values
      },
      attrs={
      "original_minlat": domain_minlat,
      "original_maxlat": domain_maxlat,
      "original_minlon": domain_minlon,
      "original_maxlon": domain_maxlon
      },
   )

   return ds

 