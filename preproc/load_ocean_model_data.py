import numpy as np
import xarray as xr
import gsw

#====================================================================================
def rename_dimensions(ds):
   """ rename dimensions of xarray dataset ds to standard (x,y,z,time) if needed

   """

   if ( "time_counter" in ds.dims ):
     ds=ds.rename({'time_counter':'time'})
   if ( "t" in ds.dims ):
     ds=ds.rename({'t':'time'})

   if ( "depth" in ds.dims ):
     ds=ds.rename({'depth':'z'})
   if ( "level" in ds.dims ):
     ds=ds.rename({'level':'z'})
   if ( "deptht" in ds.dims ):
     ds=ds.rename({'deptht':'z'})
   if ( "depthu" in ds.dims ):
     ds=ds.rename({'depthu':'z'})
   if ( "depthv" in ds.dims ):
     ds=ds.rename({'depthv':'z'})

   if ( "x_grid_T" in ds.dims ):
     ds=ds.rename({'x_grid_T':'x'})
   if ( "y_grid_T" in ds.dims ):
     ds=ds.rename({'y_grid_T':'y'})
   if ( "x_grid_U" in ds.dims ):
     ds=ds.rename({'x_grid_U':'x'})
   if ( "y_grid_U" in ds.dims ):
     ds=ds.rename({'y_grid_U':'y'})
   if ( "x_grid_V" in ds.dims ):
     ds=ds.rename({'x_grid_V':'x'})
   if ( "y_grid_V" in ds.dims ):
     ds=ds.rename({'y_grid_V':'y'})

#====================================================================================
def load_oce_mod_nemo(file_mesh_mask='mesh_mask.nc',\
                      file_bathy='bathy_meter.nc',\
                      files_gridT='nemo_grid_T.nc',\
                      files_gridU='nemo_grid_U.nc',\
                      files_gridV='nemo_grid_V.nc',\
                      files_SBC='nemo_flxT.nc',\
                      files_ice='nemo_icemod.nc',\
                      files_BSF='nemo_psi.nc',\
                      rho0=1026.0, teos10=False )
   """ Read NEMO outputs and define an xarray dataset containing 
       all variables required in MISOMIP2

       Adapted to NEMO's C-grid and standard mesh/mask variables.

       rho0 corresponds to rau0 value in NEMO's eosbn2.F90 i.e. reference volumic mass [kg m-3]

       teos10=False -> assumes the nemo outputs are in potential temperature & practical salinity (EOS80)
             =True  -> assumes the nemo outputs are in CT and AS and convert to PT and PS

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
       ds = load_oce_mod_nemo(files_gridT=fT,files_gridU=fU,files_gridV=fV,files_SBC=fS,files_ice=fI,files_BSF=fP)

   """

   ncM = xr.open_dataset(file_mesh,decode_coords=False) ; rename_dimensions(ncM)
   ncB = xr.open_dataset(file_bathy,decode_coords=False) ; rename_dimensions(ncB)
   ncT = xr.open_mfdataset(files_gridT,decode_coords=False); rename_dimensions(ncT)
   ncU = xr.open_mfdataset(files_gridU,decode_coords=False); rename_dimensions(ncU)
   ncV = xr.open_mfdataset(files_gridV,decode_coords=False); rename_dimensions(ncV)
   ncS = xr.open_mfdataset(files_SBC,decode_coords=False); rename_dimensions(ncS)
   ncI = xr.open_mfdataset(files_ice,decode_coords=False); rename_dimensions(ncI)
   ncP = xr.open_mfdataset(files_BSF,decode_coords=False); rename_dimensions(ncP)

   # Define useful masks on U, V, T grids 
   maskT = ncM.tmask.isel(t=0) # 3d mask (=1 if ocean, =0 elsewhere)
   maskU = ncM.umask.isel(t=0) # 3d mask
   maskV = ncM.vmask.isel(t=0) # 3d mask
   maskTnan = maskT.where( (maskT==1) ) # 3d mask (=1 if ocean, =nan elsewhere)
   maskUnan = maskU.where( (maskU==1) )
   maskVnan = maskV.where( (maskV==1) )
   maskTnan2d = maskT.max(dim='z').where( maskT.max(dim='z')==1 ) # 2d mask (=1 for open ocean and cavities, =nan elsewhere)
   maskUnan2d = maskU.max(dim='z').where( maskU.max(dim='z')==1 )
   maskVnan2d = maskV.max(dim='z').where( maskU.max(dim='z')==1 )

   # longitude & latitude on U, V, T grids
   lonT = ncM.glamt.isel(t=0) ; latT = ncM.gphit.isel(t=0)
   lonU = ncM.glamu.isel(t=0) ; latU = ncM.gphiu.isel(t=0)
   lonV = ncM.glamv.isel(t=0) ; latV = ncM.gphiv.isel(t=0)

   # Ocean fraction at each level:
   OCFRAC = maskT*100.0

   # 2d ice-shelf fractoin:
   ISFRAC = ncM.misf.isel(t=0)*1.e0
   ISFRAC = ISFRAC.where( (ISFRAC.values > 1.5), 0.e0 )
   ISFRAC = ISFRAC.where( (ISFRAC.values < 1.5), 100. )

   # Bathymetry (including under ice shelves) [m]
   # (if possible after NEMO's initialization, i.e. from mesh_mask) :
   if ( "bathy_meter" in ncM.data_vars ):
     OCEDEP = ncM.bathy_meter
   elif ( "bathy_metry" in ncM.data_vars ):
     OCEDEP = ncM.bathy_metry
   elif ( "Bathymetry_isf" in ncB.data_vars ):
     OCEDEP = ncB.Bathymetry_isf
   elif ( "Bathymetry" in ncB.data_vars ):
     OCEDEP = ncB.Bathymetry
   elif ( "bathy_meter" in ncB.data_vars ):
     OCEDEP = ncB.bathy_meter
  
   # Depth of ice shelf draft (if possible after NEMO's initialization, i.e. from mesh_mask) [m]:
   if ( "isfdraft" in ncM.data_vars ):
     ISFDEP = ncM.isfdraft.isel(t=0) 
   elif ( "isfdraft" in ncB.data_vars ):
     ISFDEP = ncB.isfdraft
 
   # ocean temperature [degC]
   if ( "toce" in ncT.data_vars ):
     TT = ncT.toce
   elif ( "thetao" in ncT.data_vars ):
     TT = ncT.thetao
   elif ( "votemper" in ncT.data_vars ):
     TT = ncT.votemper

   # ocean salinity [1.e-3]
   if ( "soce" in ncT.data_vars ):
     SS = ncT.soce
   if ( "so" in ncT.data_vars ):
     SS = ncT.so
   elif ( "vosaline" in ncT.data_vars ):
     SS = ncT.vosaline

   # sea bottom ocean temperature [degC]
   if ( "sbt" in ncT.data_vars ):
     TTB = ncT.sbt
   elif ( "sosbt" in ncT.data_vars ):
     TTB = ncT.sosbt
   elif ( "tob" in ncT.data_vars ):
     TTB = ncT.tob

   # sea bottom ocean salinity [1.e-3]
   if ( "sbs" in ncT.data_vars ):
     SSB = ncT.sbs
   elif ( "sosbs" in ncT.data_vars ):
     SSB = ncT.sosbs
   elif ( "sob" in ncT.data_vars ):
     SSB = ncT.sob

   # ocean potential temperature and practical salinity :
   if (teos10):
     TOB = xr.apply_ufunc(gsw.pt_from_CT, SSB, TTB)
     SOB = xr.apply_ufunc(gsw.SP_from_SA, SSB, OCEDEP, lonT, latT)
     THETAO = xr.apply_ufunc(gsw.pt_from_CT, SS, TT)
     SO = xr.apply_ufunc(gsw.SP_from_SA, SS, OCEDEP, lonT, latT)
   else: 
     TOB = TTB
     SOB = SSB
     THETAO = TT
     SO = SS
   
   # ocean x-ward velocity [m s-1]
   if ( "uoce" in ncU.data_vars ):
     UU = ncU.uoce
   elif ( "vozocrtx" in ncU.data_vars ):
     UU = ncU.vozocrtx
   elif ( "uo" in ncU.data_vars ):
     UU = ncU.uo

   # ocean y-ward velocity [m s-1]
   if ( "voce" in ncV.data_vars ):
     VV = ncV.voce
   elif ( "vomecrty" in ncV.data_vars ):
     VV = ncV.vomecrty
   elif ( "vo" in ncV.data_vars ):
     VV = ncV.vo

   # surface stress received by the ocean along x [W m-1]
   if ( "utau" in ncU.data_vars ):
     TAUUO = ncU.utau
   elif ( "sozotaux" in ncU.data_vars ):
     TAUUO = ncU.sozotaux

   # surface stress received by the ocean along x [W m-1]
   if ( "vtau" in ncV.data_vars ):
     TAUVO = ncV.vtau
   elif ( "sometauy" in ncV.data_vars ):
     TAUVO = ncV.sometauy

   # mass barotropic streamfunction
   # see Griffies et al. (2016, section H26): d(psi)/dy=-U (U: x-ward mass transport), d(psi)/dx=V (V: yward mass transport)
   MSFBA = ncP.sobarstf * rho0

   # ice shelf melt [kg m-2 s-1, positive for actual melting] :
   if ( "fwfisf" in ncS.data_vars ):
     ISFMLT = ncS.fwfisf*(-1)
   elif ( "sowflisf_cav" in ncS.data_vars ):
     ISFMLT = ncS.sowflisf_cav*(-1)

   # sea-ice concentration [0-100]
   SICONC = ncI.siconc*100.0
   SICONC = SICONC.where( (~np.isnan(SICONC.values)) & (~np.isinf(SICONC.values)), 0.e0 )

   # sea-ice volume per area [m]
   SIVOLU = ncI.sivolu

   # sea-ice x-ward velocity [m/s]
   if ( "sivelu" in ncI.data_vars ):
     SIVELU = ncI.sivelu
   elif ("siu" in ncI.data_vars ):
     SIVELU = ncI.siu

   # sea-ice y-ward velocity [m/s]
   if ( "sivelv" in ncI.data_vars ):
     SIVELV = ncI.sivelv
   elif ("siv" in ncI.data_vars ):
     SIVELV = ncI.siv


HFLUX = ncS.qt_oce + ncS.qisf # Total heat flux received by the ocean surface (including ice-shelf/ocean interface) [W m-2] 
                              # see Griffies et al. (2016, section K4-K5) NB: here, including correction if any unlike Griffies (to avoid 2 variables)
WFATRI = - ncS.empmr - ncS.fwfisf # Water flux entering the ocean due to rainfall + snowfall + evs (condensation - evap) + river runoff 
                                  # + iceberg melt + ice-shelf melt (= pr+prs+evs+ficeberg+friver+ficeshelf in Griffies 2016, section K2) [kg m-2 s-1]
ERP = ncS.erp.where( (~np.isnan(ncS.erp.values)), 0.e0 ) # surface correction (SSS restoring)
if ( "saltflx" in ncS.data_vars ): # NB: saltflx unit attribute is wrong in nico's output, it is actually in [1e-3 kg m-2 s-1]
  WFSICO = - ERP - ncS.saltflx / SS.isel(z=0) # Water flux entering the ocean due to sea-ice (melting-freezing) and surface correction (SSS restoring)
                                             # (= fsitherm + wfocorr in Griffies 2016 section K2) [kg m-2 s-1]
elif ( "sfx" in ncI.data_vars ):
  WFSICO = - ERP - ncI.sfx/86400.0 / SS.isel(z=0)

dxT = ncM.e1t.isel(t=0) # grid mesh width along x-axis in meters for sea ice grid (T)
dyT = ncM.e2t.isel(t=0) # grid mesh width along y-axis in meters for sea ice grid (T)
dxU = ncM.e1u.isel(t=0) # grid mesh width along x-axis in meters for UU
dxV = ncM.e1v.isel(t=0) # grid mesh width along x-axis in meters for VV
dlatTdy = 0.500000000 * ( latT.shift(x=1) - latT.shift(x=-1) )
dlatUdy = 0.500000000 * ( latU.shift(x=1) - latU.shift(x=-1) )
dlatVdy = 0.500000000 * ( latV.shift(x=1) - latV.shift(x=-1) )

depTUV=ncM.gdept_1d.isel(t=0) # depth of TT, SS, UU, VV (neglecting the effects of partial steps in the interpolation)

time0 = ncT.time  # central time of each month (in seconds since 1900-01-01 00:00:00)

RT = 6.371e6 # Earth radius in meters
thetaT = np.arcsin( RT*dlatTdy*np.pi/180. / dxT  ) # local T grid rotation angle compared to (zonal,meridional)
thetaU = np.arcsin( RT*dlatUdy*np.pi/180. / dxU  ) # local U grid rotation angle compared to (zonal,meridional)
thetaV = np.arcsin( RT*dlatVdy*np.pi/180. / dxV  ) # local V grid rotation angle compared to (zonal,meridional)
print('Minimum local grid angle in degrees w.r.t. (zonal,meridional):',thetaU.min().values*180./np.pi)
print('Maximum local grid angle in degrees w.r.t. (zonal,meridional):',thetaU.max().values*180./np.pi)

