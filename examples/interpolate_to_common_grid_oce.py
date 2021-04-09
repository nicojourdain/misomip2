###################################################################################
# Script to interpolate native ocean model outputs to the standard MISOMIP2 grids.
#
# Notes:
#     * interpolations are done on 1-dimensional vectors, so the script should be
#       easily adaptable to non-structured grids.
#     * needs some modificationds to work with sigma or HAL coordinates.
#     
# History: 
#     * 2021-03: initial code + tests for NEMO & MITgcm outputs [N. Jourdain, IGE-CNRS]
#
###################################################################################
import numpy as np
import xarray as xr
import sys
sys.path.append("/Users/jourdain/MY_SCRIPTS")
#sys.path.append("/Users/nakayama/Documents/GitHub/")
import misomip2.preproc as mp
import os
from datetime import datetime

startTime = datetime.now()

np.seterr(divide='ignore', invalid='ignore') # to avoid warning due to divide by zero

#--------------------------------------------------------------------------
# 0- General information:

# Official name in MISOMIP2:
model='NEMO_test'
#model='MITGCM_test'

reg='Amundsen' # 'Amundsen' or 'Weddell'

exp='A1' # MISOMIP2 experiment ('A1', 'W1', 'A2', ...)

data_dir='models/oce/'+model

missval=9.969209968386869e36 # missing value in created netcdf files

epsfr=1. # cut sea, sea-ice or ice-shelf fractions below epsfr (in %)
         # (useful to avoid extreme values that can occur with very small fractions)

#--------------------------------------------------------------------------
# 1- Files and variables

# loading an xarray dataset containing all useful variables with (x,y) reshaped 
# as 1 dimension in case of original structured grid :

if ( model[0:4] == 'NEMO' ):

   print('LOADING NEMO...')
   f_mesh   = data_dir+'/mesh_mask_test.nc'
   f_bathy  = data_dir+'/bathy_meter_test.nc'
   fs_gridT = [data_dir+'/'+model+'_m0'+month.astype('str')+'_grid_T.nc' for month in np.arange(1,3)]
   fs_gridU = [data_dir+'/'+model+'_m0'+month.astype('str')+'_grid_U.nc' for month in np.arange(1,3)]
   fs_gridV = [data_dir+'/'+model+'_m0'+month.astype('str')+'_grid_V.nc' for month in np.arange(1,3)]
   fs_SBC   = [data_dir+'/'+model+'_m0'+month.astype('str')+'_SBC.nc' for month in np.arange(1,3)]
   fs_ice   = [data_dir+'/'+model+'_m0'+month.astype('str')+'_icemod.nc' for month in np.arange(1,3)]
   # Barotropic Streamfunction calculated at U points using the cdfpsi function which is part of the cdftools (https://github.com/meom-group/CDFTOOLS):
   fs_BSF   = [data_dir+'/'+model+'_m0'+month.astype('str')+'_psi.nc' for month in np.arange(1,3)]

   oce = mp.load_oce_mod_nemo( file_mesh_mask=f_mesh, file_bathy=f_bathy, files_gridT=fs_gridT,\
                  files_gridU=fs_gridU, files_gridV=fs_gridV, files_SBC=fs_SBC, files_ice=fs_ice,\
                  files_BSF=fs_BSF, rho0=1026.0, teos10=True )

elif ( model[0:6] == 'MITGCM' ):

   print('LOADING MITGCM...')
   fT = data_dir+'/'+model+'_THETA.nc'
   fS = data_dir+'/'+model+'_SALT.nc'
   fU = data_dir+'/'+model+'_UVEL.nc'
   fV = data_dir+'/'+model+'_VVEL.nc'

   oce = mp.load_oce_mod_mitgcm( files_T=fT, files_S=fS, files_U=fU, files_V=fV, rho0=1026.0, teos10=False )

else:

   sys.exit('Unknown model ==> Write a function to load this model outputs')

#--------------------------------------------------------------------------
# 2- Global attributes of output netcdf :

def put_global_attrs(ds,experiment='TBD',avg_hor_res_73S=0.0,original_sim_name='None',\
                     original_min_lat=-90.0,original_max_lat=90.0,original_min_lon=-180.0,original_max_lon=180.0):
   """ Put global attributes to the ds xarray dataset

   """
   ds.attrs['project'] = 'MISOMIP2'
   ds.attrs['contact'] = 'Nicolas Jourdain <nicolas.jourdain@univ-grenoble-alpes.fr>' # name <email>
   ds.attrs['institute'] = 'CNRS-UGA-IGE'
   ds.attrs['computing_facility'] = 'occigen-CINES'                    # Computing center where the simulation was run
   ds.attrs['interpolation_method'] = 'linear triangular barycentric'  # in: 'linear triangular barycentric', 'bi-linear',
                                                                       #     'nearest-neighbor', 'conservative', ...       
   ds.attrs['ocean_model'] = 'NEMO3.6'                                 # Model name and version
   ds.attrs['reference'] = 'Jourdain et al. 2019 (doi:10.1016/j.ocemod.2018.11.001)'  # publication describing the simulation or a similar configuration
   ds.attrs['original_sim_name'] = original_sim_name                   # original simulation name
   ds.attrs['experiment'] = experiment                                 # in: 'A1', 'W1', 'A2', ...
   ds.attrs['bathymetry'] = 'BedMachine-v1.33'                         # Bathymetry dataset (specify exact version)
   ds.attrs['ice_shelf_draft'] = 'BedMachine-v1.33'                    # Ice draft depth dataset (specify exact version)
   ds.attrs['atmosphere'] = 'DFS5.2'                                   # in: 'ERA5', 'ERAint', 'CORE', 'MERRA2', 'JRA55do', ...
   ds.attrs['iceberg'] = 'Prescribed Freshwater'                       # in: 'Lagrangian Model', 'Prescribed Freshwater', 
                                                                       #     'Prescribed Freshwater and Heat', 'None'
   ds.attrs['sea_ice'] = 'Dynamics-Thermodynamics Model'               # in: 'Dynamics-Thermodynamics Model', 'Thermodynamics Model', 'Prescribed Freshwater and Heat'
   ds.attrs['ocean_lateral_bdy'] = 'Simulation'                        # in: 'None', 'Observational Data', 'Ocean Reanalysis', 'Simulation', 'Simulation with corrections' 
   ds.attrs['tides'] = 'Resolved (prescribed at bdy)'                  # in: 'None', 'Resolved (prescribed at bdy)', 'Resolved (tidal potential)',
                                                                       #     'Parameterized (uniform tidal velocity)', 'Parameterized (non-uniform tidal velocity),
                                                                       #     'Parameterized (other)'
   ds.attrs['vertical_coordinate'] = 'Stretched Geopotential (Zstar)'  # in: 'Geopotential (Z)', 'Stretched Geopotential (Zstar)', 'Pressure (P)', 
                                                                       #     'Stretched Pressure (P*)', 'Isopycnal', 'Terrain-Following (Sigma)', 
                                                                       #     'Arbitrary Lagrangian-Eulerian (ALE)'
   ds.attrs['is_melt_param'] = '3-equation (velocity-dependent gamma)' # in: '3-equation (constant gamma)', '3-equation (velocity-dependent gamma)',
                                                                       #     '3-equation (stability and velocity-dependent gamma)', ...
   ds.attrs['avg_hor_res_73S'] = avg_hor_res_73S                       # Average horizontal resolution (m) at 73degS in the MISOMIP2 domain (average of x and y resolution)
   ds.attrs['original_min_lat'] = original_min_lat                     # Minimum latitude of the original domain in [-180:180]
   ds.attrs['original_max_lat'] = original_max_lat                     # Minimum latitude of the original domain in [-180:180]
   ds.attrs['original_min_lon'] = original_min_lon                     # Minimum longitude of the original domain in [-90:90]
   ds.attrs['original_max_lon'] = original_max_lon                     # Minimum longitude of the original domain in [-90:90]

#--------------------------------------------------------------------------
# 3a- Interpolate to common 3d grid :


# Characteristics of MISOMIP 3d grid:
[lon_miso,lat_miso,dep_miso] = mp.generate_3d_grid_oce(region=reg)
mlon = np.size(lon_miso)
mlat = np.size(lat_miso)
mdep = np.size(dep_miso)
lon2d_miso, lat2d_miso = np.meshgrid( lon_miso, lat_miso )
lon_miso1d = np.reshape( lon2d_miso, mlon*mlat )
lat_miso1d = np.reshape( lat2d_miso, mlon*mlat )

# model coordinates:
lonT=oce.lonT ; lonU=oce.lonU ; lonV=oce.lonV
latT=oce.latT ; latU=oce.latU ; latV=oce.latV

# Some quantities needed to define the global attributes:
res_73S = 0.5 * (  oce.dxT.where(  (latT < -72.5) & (latT >= -73.5) & (lonT >= lon_miso.min()) & (lonT <= lon_miso.max()) \
                                 & (oce.LEVOFT.isel(z=0)>99.)   ).mean().values \
                 + oce.dyT.where(  (latT < -72.5) & (latT >= -73.5) & (lonT >= lon_miso.min()) & (lonT <= lon_miso.max()) \
                                 & (oce.LEVOFT.isel(z=0)>99.) ).mean().values )
print('Average horizontal resolution at 73S :',res_73S)

# 3d sea fraction
LEVOFT=oce.LEVOFT
LEVOFU=oce.LEVOFU
LEVOFV=oce.LEVOFV
# 2d sea fraction (including under-ice-shelf seas)
seafracT=LEVOFT.max('z',skipna=True) # 100 if any ocean point, including under ice shelves, =0 or =nan if no ocean point
seafracU=LEVOFU.max('z',skipna=True)
seafracV=LEVOFV.max('z',skipna=True)
# 2d sea fraction (NOT including under-ice-shelf seas)
openseafracT=LEVOFT.isel(z=0) # 100 only if open ocean point
openseafracU=LEVOFU.isel(z=0)
openseafracV=LEVOFV.isel(z=0)

# Lower and upper indices for vertical interpolation:
[kinf,ksup] = mp.vertical_interp(oce.depTUV.values,dep_miso)

mtime = np.shape(oce.SO)[0]

# mask showing the original domain on the MISOMIP grid (=nan where interpolation of any of T, U, V grid is nan, =1 elsewhere):
DOMMSK_miso =               mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.DOMMSKT )
DOMMSK_miso = DOMMSK_miso + mp.horizontal_interp( latU, lonU, mlat, mlon, lat_miso1d, lon_miso1d, oce.DOMMSKU )
DOMMSK_miso = DOMMSK_miso + mp.horizontal_interp( latV, lonV, mlat, mlon, lat_miso1d, lon_miso1d, oce.DOMMSKV )

# horizontal interpolation grid roation angles :
theT = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.thetaT )
theU = mp.horizontal_interp( latU, lonU, mlat, mlon, lat_miso1d, lon_miso1d, oce.thetaU )
theV = mp.horizontal_interp( latV, lonV, mlat, mlon, lat_miso1d, lon_miso1d, oce.thetaV )

# Ice shelf fraction on MISOMIP grid (in [0-100], =nan beyond model domain)
SFTFLI_miso = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.SFTFLI )
SFTFLI_miso[ np.isnan(DOMMSK_miso) ] = np.nan # will update to missval at the end of the calculation

# Depth of floating ice on MISOMIP grid (ice-shelf draft) (=nan where no ice shelf or beyond model domain)
DEPFLI_miso = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.DEPFLI, weight=oce.SFTFLI*oce.DOMMSKT, skipna=True, filnocvx=True, threshold=epsfr )
DEPFLI_miso[ np.isnan(DOMMSK_miso) | np.isnan(DEPFLI_miso) ] = missval

# Ocean depth on MISOMIP grid (=nan where land or grounded ice or beyond model domain)
DEPTHO_miso = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.DEPTHO, weight=seafracT, skipna=True, filnocvx=True, threshold=epsfr )
DEPTHO_miso[ np.isnan(DOMMSK_miso) | np.isnan(DEPTHO_miso) ] = missval 

# Sea area fraction at each vertical level on the MISOMIP grid: 
LEVOF_miso = np.zeros((mdep,mlat,mlon))
tmp_LEVT = np.zeros((mdep,LEVOFT.shape[1]))
tmp_LEVU = np.zeros((mdep,LEVOFU.shape[1]))
tmp_LEVV = np.zeros((mdep,LEVOFV.shape[1]))
for kk in np.arange(mdep):
  if ( kinf[kk] == ksup[kk] ):
    tmpaT = 1.e0
    tmpbT = 0.e0
  else:
    tmpaT = oce.depTUV.isel(z=ksup[kk]) - dep_miso[kk]
    tmpbT = dep_miso[kk] - oce.depTUV.isel(z=kinf[kk])
  tmpxT = tmpaT + tmpbT
  tmp_OC = ( LEVOFT.isel(z=kinf[kk]) * tmpaT + LEVOFT.isel(z=ksup[kk]) * tmpbT ) / tmpxT
  tmp_LEVT[kk,:] = tmp_OC # =LEVOF innterpolated vertically but not horizontally
  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, tmp_OC )
  tzz[ np.isnan(DOMMSK_miso) ] = np.nan # will update to missval at the end of the calculation
  tzz[ (DEPTHO_miso<dep_miso[kk]) ]=0.0 # criteria on DEPTHO_miso to account for partial step if needed.
  LEVOF_miso[kk,:,:] = tzz # saved to netcdf
  tmp_OC = ( LEVOFU.isel(z=kinf[kk]) * tmpaT + LEVOFU.isel(z=ksup[kk]) * tmpbT ) / tmpxT
  tmp_LEVU[kk,:] = tmp_OC
  tmp_OC = ( LEVOFV.isel(z=kinf[kk]) * tmpaT + LEVOFV.isel(z=ksup[kk]) * tmpbT ) / tmpxT
  tmp_LEVV[kk,:] = tmp_OC

#----- interpolation of time-varying fields:

SO_miso     = np.zeros((mtime,mdep,mlat,mlon)) + missval
THETAO_miso = np.zeros((mtime,mdep,mlat,mlon)) + missval
UO_miso     = np.zeros((mtime,mdep,mlat,mlon)) + missval
VO_miso     = np.zeros((mtime,mdep,mlat,mlon)) + missval
ZOS_miso    = np.zeros((mtime,mlat,mlon)) + missval
TOB_miso    = np.zeros((mtime,mlat,mlon)) + missval
SOB_miso    = np.zeros((mtime,mlat,mlon)) + missval
FICESHELF_miso = np.zeros((mtime,mlat,mlon)) + missval
DYDRFLI_miso = np.zeros((mtime,mlat,mlon)) + missval
THDRFLI_miso = np.zeros((mtime,mlat,mlon)) + missval
HADRFLI_miso = np.zeros((mtime,mlat,mlon)) + missval
MSFTBAROT_miso = np.zeros((mtime,mlat,mlon)) + missval
HFDS_miso = np.zeros((mtime,mlat,mlon)) + missval
WFOATRLI_miso = np.zeros((mtime,mlat,mlon)) + missval
WFOSICOR_miso = np.zeros((mtime,mlat,mlon)) + missval
SICONC_miso = np.zeros((mtime,mlat,mlon)) + missval
SIVOL_miso = np.zeros((mtime,mlat,mlon)) + missval
SIU_miso = np.zeros((mtime,mlat,mlon)) + missval
SIV_miso = np.zeros((mtime,mlat,mlon)) + missval
TAUUO_miso = np.zeros((mtime,mlat,mlon)) + missval
TAUVO_miso = np.zeros((mtime,mlat,mlon)) + missval

domcond = (np.isnan(DOMMSK_miso))

for ll in np.arange(mtime):

  ## fields interpolated from sea cells ( C/T grid ):

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.ZOS.isel(time=ll), weight=seafracT, skipna=True, filnocvx=True, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  ZOS_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.TOB.isel(time=ll), weight=seafracT, skipna=True, filnocvx=True, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  TOB_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.SOB.isel(time=ll), weight=seafracT, skipna=True, filnocvx=True, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  SOB_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.HFDS.isel(time=ll), weight=seafracT, skipna=True, filnocvx=True, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  HFDS_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.WFOATRLI.isel(time=ll), weight=seafracT, skipna=True, filnocvx=True, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  WFOATRLI_miso[ll,:,:] = tzz

  ## fields interpolated from sea cells ( U & V grids ):

  tzz = mp.horizontal_interp( latU, lonU, mlat, mlon, lat_miso1d, lon_miso1d, oce.MSFTBAROT.isel(time=ll), weight=seafracU, skipna=True, filnocvx=True, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  MSFTBAROT_miso[ll,:,:] = tzz

  TAUX_notrot = mp.horizontal_interp( latU, lonU, mlat, mlon, lat_miso1d, lon_miso1d, oce.TAUX.isel(time=ll), weight=seafracU, skipna=True, filnocvx=True, threshold=epsfr )
  TAUY_notrot = mp.horizontal_interp( latV, lonV, mlat, mlon, lat_miso1d, lon_miso1d, oce.TAUY.isel(time=ll), weight=seafracV, skipna=True, filnocvx=True, threshold=epsfr )
  TAUUO_miso[ll,:,:] = TAUX_notrot * np.cos(theU) + TAUY_notrot * np.sin(theV) # rotated to zonal
  TAUVO_miso[ll,:,:] = TAUY_notrot * np.cos(theV) - TAUX_notrot * np.sin(theU) # rotated to meridional
  TAUUO_miso[ ( np.isnan(TAUUO_miso) ) | domcond ] = missval
  TAUVO_miso[ ( np.isnan(TAUVO_miso) ) | domcond ] = missval

  ## fileds interpolated from open-sea cells :

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.WFOSICOR.isel(time=ll), weight=openseafracT, skipna=True, filnocvx=True, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  WFOSICOR_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.SICONC.isel(time=ll), weight=openseafracT, skipna=True, filnocvx=True, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = np.nan # will be updated to missval later on
  SICONC_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.SIVOL.isel(time=ll), weight=openseafracT, skipna=True, filnocvx=True, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  SIVOL_miso[ll,:,:] = tzz

  ## fileds interpolated from ice-shelf cells :

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.FICESHELF.isel(time=ll), weight=oce.SFTFLI, skipna=True, filnocvx=True, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  FICESHELF_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.DYDRFLI.isel(time=ll), weight=oce.SFTFLI, skipna=True, filnocvx=True, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  DYDRFLI_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.THDRFLI.isel(time=ll), weight=oce.SFTFLI, skipna=True, filnocvx=True, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  THDRFLI_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.HADRFLI.isel(time=ll), weight=oce.SFTFLI, skipna=True, filnocvx=True, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  HADRFLI_miso[ll,:,:] = tzz

  ## fileds interpolated from sea-ice cells :

  SIUX_notrot = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.SIUX.isel(time=ll), weight=oce.SICONC.isel(time=ll), skipna=True, filnocvx=True, threshold=epsfr )
  SIVY_notrot = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, oce.SIVY.isel(time=ll), weight=oce.SICONC.isel(time=ll), skipna=True, filnocvx=True, threshold=epsfr )
  SIU_miso[ll,:,:] = SIUX_notrot * np.cos(theT) + SIVY_notrot * np.sin(theT) # rotated to zonal
  SIV_miso[ll,:,:] = SIVY_notrot * np.cos(theT) - SIUX_notrot * np.sin(theT) # rotated to meridional
  SIU_miso [ np.isnan(SIU_miso) | domcond ] = missval
  SIV_miso [ np.isnan(SIV_miso) | domcond ] = missval

  for kk in np.arange(mdep):

    # vertical interpolation to common vertical grid :
    if ( kinf[kk] == ksup[kk] ):
      tmpaT = LEVOFT.isel(z=kinf[kk])
      tmpbT = tmpaT*0
    else:
      tmpaT = LEVOFT.isel(z=kinf[kk]) * (oce.depTUV.isel(z=ksup[kk])-dep_miso[kk])
      tmpbT = LEVOFT.isel(z=ksup[kk]) * (dep_miso[kk]-oce.depTUV.isel(z=kinf[kk]))
    tmpxT = tmpaT + tmpbT
    tmp_SS = ( oce.SO.isel(time=ll,z=kinf[kk]) * tmpaT + oce.SO.isel(time=ll,z=ksup[kk]) * tmpbT ) / tmpxT # Inf if no interpolable value
    tmp_TT = ( oce.THETAO.isel(time=ll,z=kinf[kk]) * tmpaT + oce.THETAO.isel(time=ll,z=ksup[kk]) * tmpbT ) / tmpxT
   
    if ( kinf[kk] == ksup[kk] ):
      tmpaU = LEVOFU.isel(z=kinf[kk])
      tmpbU = tmpaU*0
    else:
      tmpaU = LEVOFU.isel(z=kinf[kk]) * (oce.depTUV.isel(z=ksup[kk])-dep_miso[kk])
      tmpbU = LEVOFU.isel(z=ksup[kk]) * (dep_miso[kk]-oce.depTUV.isel(z=kinf[kk]))
    tmpxU = tmpaU + tmpbU
    tmp_UX = ( oce.UX.isel(time=ll,z=kinf[kk]) * tmpaU + oce.UX.isel(time=ll,z=ksup[kk]) * tmpbU ) / tmpxU # Inf if no interpolable value 
 
    if ( kinf[kk] == ksup[kk] ):
      tmpaV = LEVOFV.isel(z=kinf[kk])
      tmpbV = tmpaV*0
    else:
      tmpaV = LEVOFV.isel(z=kinf[kk]) * (oce.depTUV.isel(z=ksup[kk])-dep_miso[kk])
      tmpbV = LEVOFV.isel(z=ksup[kk]) * (dep_miso[kk]-oce.depTUV.isel(z=kinf[kk]))
    tmpxV = tmpaV + tmpbV
    tmp_VY = ( oce.VY.isel(time=ll,z=kinf[kk]) * tmpaV + oce.VY.isel(time=ll,z=ksup[kk]) * tmpbV ) / tmpxV # Inf if no interpolable value

    # horizontal interpolations of time-varying 3d fields to common horizontal grid : 

    condkk = ( (LEVOF_miso[kk,:,:] < epsfr) | domcond )

    tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, tmp_SS, weight=tmp_LEVT[kk,:], skipna=True, filnocvx=True, threshold=epsfr )
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    SO_miso[ll,kk,:,:] = tzz

    tzz = mp.horizontal_interp( latT, lonT, mlat, mlon, lat_miso1d, lon_miso1d, tmp_TT, weight=tmp_LEVT[kk,:], skipna=True, filnocvx=True, threshold=epsfr )
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    THETAO_miso[ll,kk,:,:] = tzz

    UX_notrot = mp.horizontal_interp( latU, lonU, mlat, mlon, lat_miso1d, lon_miso1d, tmp_UX, weight=tmp_LEVU[kk,:], skipna=True, filnocvx=True, threshold=epsfr )
    VY_notrot = mp.horizontal_interp( latV, lonV, mlat, mlon, lat_miso1d, lon_miso1d, tmp_VY, weight=tmp_LEVV[kk,:], skipna=True, filnocvx=True, threshold=epsfr )
    tzz = UX_notrot * np.cos(theU) + VY_notrot * np.sin(theV) # rotated to zonal
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    UO_miso[ll,kk,:,:] = tzz
    tzz = VY_notrot * np.cos(theV) - UX_notrot * np.sin(theU) # rotated to meridional
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    VO_miso[ll,kk,:,:] = tzz


LEVOF_miso[ np.isnan(LEVOF_miso) ] = missval
SFTFLI_miso[ np.isnan(SFTFLI_miso) ] = missval 
SICONC_miso[ np.isnan(SICONC_miso) ] = missval

#--------------------------------------------------------------------------
# 3b- Create new xarray dataset and save to netcdf

dsmiso3d = xr.Dataset(
    {
    "so":        (["time", "depth", "latitude", "longitude"], np.float32(SO_miso)),
    "thetao":    (["time", "depth", "latitude", "longitude"], np.float32(THETAO_miso)),
    "uo":        (["time", "depth", "latitude", "longitude"], np.float32(UO_miso)),
    "vo":        (["time", "depth", "latitude", "longitude"], np.float32(VO_miso)),
    "tauuo":     (["time", "latitude", "longitude"], np.float32(TAUUO_miso)),
    "tauvo":     (["time", "latitude", "longitude"], np.float32(TAUVO_miso)),
    "zos":       (["time", "latitude", "longitude"], np.float32(ZOS_miso)),
    "tob":       (["time", "latitude", "longitude"], np.float32(TOB_miso)),
    "sob":       (["time", "latitude", "longitude"], np.float32(SOB_miso)),
    "ficeshelf": (["time", "latitude", "longitude"], np.float32(FICESHELF_miso)),
    "dydrfli":   (["time", "latitude", "longitude"], np.float32(DYDRFLI_miso)),
    "thdrfli":   (["time", "latitude", "longitude"], np.float32(THDRFLI_miso)),
    "hadrfli":   (["time", "latitude", "longitude"], np.float32(HADRFLI_miso)),
    "msftbarot": (["time", "latitude", "longitude"], np.float32(MSFTBAROT_miso)),
    "hfds":      (["time", "latitude", "longitude"], np.float32(HFDS_miso)),
    "wfoatrli":  (["time", "latitude", "longitude"], np.float32(WFOATRLI_miso)),
    "wfosicor":  (["time", "latitude", "longitude"], np.float32(WFOSICOR_miso)),
    "siconc":    (["time", "latitude", "longitude"], np.float32(SICONC_miso)),
    "sivol":     (["time", "latitude", "longitude"], np.float32(SIVOL_miso)),
    "siu":       (["time", "latitude", "longitude"], np.float32(SIU_miso)),
    "siv":       (["time", "latitude", "longitude"], np.float32(SIV_miso)),
    "levof":     (["depth", "latitude", "longitude"], np.float32(LEVOF_miso)),
    "sftfli":    (["latitude", "longitude"], np.float32(SFTFLI_miso)),
    "depfli":    (["latitude", "longitude"], np.float32(DEPFLI_miso)),
    "deptho":    (["latitude", "longitude"], np.float32(DEPTHO_miso)),
    },
    coords={
    "longitude":np.float32(lon_miso),
    "latitude": np.float32(lat_miso),
    "depth": np.float32(dep_miso),
    "time": oce.time.values
    },
)

mp.add_standard_attributes_oce(dsmiso3d,miss=missval)
put_global_attrs(dsmiso3d,experiment=exp,avg_hor_res_73S=res_73S,original_sim_name='misomip2_test_case',\
                 original_min_lat=oce.attrs.get('original_minlat'),original_max_lat=oce.attrs.get('original_maxlat'),\
                 original_min_lon=oce.attrs.get('original_minlon'),original_max_lon=oce.attrs.get('original_maxlon') )

file_miso3d = 'Oce3d_'+model+'_'+exp+'.nc'

print('Creating ',file_miso3d)
dsmiso3d.to_netcdf(file_miso3d,unlimited_dims="time")

del dsmiso3d
del SO_miso, THETAO_miso, UO_miso, VO_miso, ZOS_miso, FICESHELF_miso, DYDRFLI_miso, THDRFLI_miso, HADRFLI_miso, MSFTBAROT_miso
del HFDS_miso, WFOATRLI_miso, WFOSICOR_miso, SICONC_miso, SIVOL_miso, SIU_miso, SIV_miso

print('  Execution time: ',datetime.now() - startTime)

#--------------------------------------------------------------------------
# 4a- Interpolate to common section :

# Characteristics of MISOMIP section grid:
[lon_sect1d,lat_sect1d,dep_sect] = mp.generate_section_grid_oce(region=reg)
mlonlatsec = np.size(lon_sect1d)
mdepsect = np.size(dep_sect)

# Reduce input data size
lonmin_sect = np.min( lon_sect1d ) - 1.1
lonmax_sect = np.max( lon_sect1d ) + 1.1
latmin_sect = np.min( lat_sect1d ) - 1.1
latmax_sect = np.max( lat_sect1d ) + 1.1

cond_secT=( (oce.latT>latmin_sect) & (oce.latT<latmax_sect) & (oce.lonT>lonmin_sect) & (oce.lonT<lonmax_sect) )

# model coordinates:
lonT=oce.lonT.where(cond_secT,drop=True)
latT=oce.latT.where(cond_secT,drop=True)

# 3d sea fraction
LEVOFT=oce.LEVOFT.where(cond_secT,drop=True)
# 2d sea fraction (including under-ice-shelf seas)
seafracT=LEVOFT.max('z',skipna=True) # 100 if any ocean point, including under ice shelves, =0 or =nan if no ocean point
# 2d sea fraction (NOT including under-ice-shelf seas)
openseafracT=LEVOFT.isel(z=0) # 100 only if open ocean point

# Lower and upper indices for vertical interpolation:
[kinf,ksup] = mp.vertical_interp(oce.depTUV.values,dep_sect)

# mask showing the original domain (nan where interpolation of any of T, U, V grid is nan):
DOMMSK_sect = np.squeeze(mp.horizontal_interp( latT, lonT, mlonlatsec, 1, lat_sect1d, lon_sect1d, oce.DOMMSKT.where(cond_secT,drop=True) ))

domcond = ( np.isnan(DOMMSK_sect) )

# Ocean depth on MISOMIP grid (=nan where land or grounded ice or beyond model domain)
DEPTHO_sect = np.squeeze( mp.horizontal_interp( latT, lonT, mlonlatsec, 1, lat_sect1d, lon_sect1d, oce.DEPTHO.where(cond_secT,drop=True), \
                                                weight=seafracT, skipna=True, filnocvx=True, threshold=epsfr ))
DEPTHO_sect[ domcond | np.isnan(DEPTHO_sect) ] = missval

# vertical then horizontal interpolation of constant 3d fields to common grid :
LEVOF_sect = np.zeros((mdepsect,mlonlatsec)) + missval
tmp_LEVT = np.zeros((mdepsect,LEVOFT.shape[1]))
for kk in np.arange(mdepsect):
  if ( kinf[kk] == ksup[kk] ):
    tmpaT = 1.e0
    tmpbT = 0.e0
  else:
    tmpaT = oce.depTUV.isel(z=ksup[kk]) - dep_sect[kk]
    tmpbT = dep_sect[kk] - oce.depTUV.isel(z=kinf[kk])
  tmpxT = tmpaT + tmpbT
  tmp_OC = ( LEVOFT.isel(z=kinf[kk]) * tmpaT + LEVOFT.isel(z=ksup[kk]) * tmpbT ) / tmpxT
  tmp_LEVT[kk,:] = tmp_OC # LEVOF innterpolated vertically but not horizontally
  tzz = np.squeeze(mp.horizontal_interp( latT, lonT, mlonlatsec, 1, lat_sect1d, lon_sect1d, tmp_OC ))
  tzz[ domcond ] = missval
  tzz[ (DEPTHO_sect>dep_sect[kk]) ]=0.0 # criteria on DEPTHO_sect to account for partial step if needed.
  LEVOF_sect[kk,:] = tzz

SO_sect     = np.zeros((mtime,mdepsect,mlonlatsec)) + missval
THETAO_sect = np.zeros((mtime,mdepsect,mlonlatsec)) + missval

for ll in np.arange(mtime):
  for kk in np.arange(mdepsect):

    condkk = ( (LEVOF_sect[kk,:] < epsfr) | domcond )

    if ( kinf[kk] == ksup[kk] ):
      tmpaT = LEVOFT.isel(z=kinf[kk])
      tmpbT = tmpaT*0
    else:
      tmpaT = LEVOFT.isel(z=kinf[kk]) * (oce.depTUV.isel(z=ksup[kk])-dep_sect[kk])
      tmpbT = LEVOFT.isel(z=ksup[kk]) * (dep_sect[kk]-oce.depTUV.isel(z=kinf[kk]))
    tmpxT = tmpaT + tmpbT
    tmp_SS = ( oce.SO.where(cond_secT,drop=True).isel(time=ll,z=kinf[kk]) * tmpaT + oce.SO.where(cond_secT,drop=True).isel(time=ll,z=ksup[kk]) * tmpbT ) / tmpxT # Inf if no interpolable value
    tmp_TT = ( oce.THETAO.where(cond_secT,drop=True).isel(time=ll,z=kinf[kk]) * tmpaT + oce.THETAO.where(cond_secT,drop=True).isel(time=ll,z=ksup[kk]) * tmpbT ) / tmpxT

    tzz = np.squeeze(mp.horizontal_interp( latT, lonT, mlonlatsec, 1, lat_sect1d, lon_sect1d, tmp_SS, weight=tmp_LEVT[kk,:], skipna=True, filnocvx=True, threshold=epsfr ) )
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    SO_sect[ll,kk,:] = tzz

    tzz = np.squeeze(mp.horizontal_interp( latT, lonT, mlonlatsec, 1, lat_sect1d, lon_sect1d, tmp_TT, weight=tmp_LEVT[kk,:], skipna=True, filnocvx=True, threshold=epsfr ) )
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    THETAO_sect[ll,kk,:] = tzz

#--------------------------------------------------------------------------
# 4b- Create new xarray dataset and save to netcdf

dssect = xr.Dataset(
    {
    "so":        (["time", "depth", "x"], np.float32(SO_sect)),
    "thetao":    (["time", "depth", "x"], np.float32(THETAO_sect)),
    "levof":     (["depth", "x"], np.float32(LEVOF_sect)),
    "longitude": (["x"], np.float32(lon_sect1d)),
    "latitude": (["x"], np.float32(lat_sect1d)),
    },
    coords={
    "depth": np.float32(dep_sect),
    "time": oce.time.values
    },
)

mp.add_standard_attributes_oce(dssect,miss=missval)
put_global_attrs(dssect,experiment=exp,avg_hor_res_73S=res_73S,original_sim_name='misomip2_test_case',\
                 original_min_lat=oce.attrs.get('original_minlat'),original_max_lat=oce.attrs.get('original_maxlat'),\
                 original_min_lon=oce.attrs.get('original_minlon'),original_max_lon=oce.attrs.get('original_maxlon') )

file_sect = 'OceSec_'+model+'_'+exp+'.nc'

print('Creating ',file_sect)
dssect.to_netcdf(file_sect,unlimited_dims="time")

del dssect
del SO_sect, THETAO_sect

print('  Execution time: ',datetime.now() - startTime)

#--------------------------------------------------------------------------
# 5a- Interpolate to common mooring location :


# Characteristics of MISOMIP mooring:
[lon_moor0d,lat_moor0d,dep_moor] = mp.generate_mooring_grid_oce(region=reg)
mdepmoor = np.size(dep_moor)

# Reduce input data size
lonmin_moor = np.min( lon_moor0d ) - 1.1
lonmax_moor = np.max( lon_moor0d ) + 1.1
latmin_moor = np.min( lat_moor0d ) - 1.1
latmax_moor = np.max( lat_moor0d ) + 1.1

cond_mooT=( (oce.latT>latmin_moor) & (oce.latT<latmax_moor) & (oce.lonT>lonmin_moor) & (oce.lonT<lonmax_moor) )

# model coordinates:
lonT=oce.lonT.where(cond_mooT,drop=True)
latT=oce.latT.where(cond_mooT,drop=True)

# 3d sea fraction
LEVOFT=oce.LEVOFT.where(cond_mooT,drop=True)
# 2d sea fraction (including under-ice-shelf seas)
seafracT=LEVOFT.max('z',skipna=True) # 100 if any ocean point, including under ice shelves, =0 or =nan if no ocean point
# 2d sea fraction (NOT including under-ice-shelf seas)
openseafracT=LEVOFT.isel(z=0) # 100 only if open ocean point

# Lower and upper indices for vertical interpolation:
[kinf,ksup] = mp.vertical_interp(oce.depTUV.values,dep_moor)

mtime = np.shape(oce.SO)[0]

# mask showing the original domain (nan where interpolation of any of T, U, V grid is nan):
DOMMSK_moor = np.squeeze(mp.horizontal_interp( latT, lonT, 1, 1, lat_moor0d, lon_moor0d, oce.DOMMSKT.where(cond_mooT,drop=True) ))

# Ocean depth on MISOMIP grid (=nan where land or grounded ice or beyond model domain)
DEPTHO_moor = np.squeeze(mp.horizontal_interp( latT, lonT, 1, 1, lat_moor0d, lon_moor0d, oce.DEPTHO.where(cond_mooT,drop=True), \
                                               weight=seafracT, skipna=True, filnocvx=True, threshold=epsfr ) )
DEPTHO_moor[ np.isnan(DOMMSK_moor) | np.isnan(DEPTHO_moor) ] = missval

# vertical then horizontal interpolation of constant 3d fields to common grid :
LEVOF_moor = np.zeros((mdepmoor)) + missval
tmp_LEVT = np.zeros((mdepmoor,LEVOFT.shape[1]))
for kk in np.arange(mdepmoor):
  if ( kinf[kk] == ksup[kk] ):
    tmpaT = 1.e0
    tmpbT = 0.e0
  else:
    tmpaT = oce.depTUV.isel(z=ksup[kk]) - dep_moor[kk]
    tmpbT = dep_moor[kk] - oce.depTUV.isel(z=kinf[kk])
  tmpxT = tmpaT + tmpbT
  tmp_OC = ( LEVOFT.isel(z=kinf[kk]) * tmpaT + LEVOFT.isel(z=ksup[kk]) * tmpbT ) / tmpxT
  tmp_LEVT[kk,:] = tmp_OC # LEVOF innterpolated vertically but not horizontally
  tzz = np.squeeze(mp.horizontal_interp( latT, lonT, 1, 1, lat_moor0d, lon_moor0d, tmp_OC ))
  if ( np.isnan(DOMMSK_moor) ):
     tzz = missval
  tzz[ (DEPTHO_moor<dep_moor[kk]) ]=0.0 # criteria on DEPTHO_moor to account for partial step if needed.
  LEVOF_moor[kk] = tzz
LEVOF_moor[ LEVOF_moor < epsfr ] = 0.e0
tmp_LEVT[ tmp_LEVT < epsfr ] = 0.e0

SO_moor     = np.zeros((mtime,mdepmoor)) + missval
THETAO_moor = np.zeros((mtime,mdepmoor)) + missval

for ll in np.arange(mtime):
  for kk in np.arange(mdepmoor):

    if ( kinf[kk] == ksup[kk] ):
      tmpaT = LEVOFT.isel(z=kinf[kk])
      tmpbT = tmpaT*0
    else:
      tmpaT = LEVOFT.isel(z=kinf[kk]) * (oce.depTUV.isel(z=ksup[kk])-dep_moor[kk])
      tmpbT = LEVOFT.isel(z=ksup[kk]) * (dep_moor[kk]-oce.depTUV.isel(z=kinf[kk]))
    tmpxT = tmpaT + tmpbT
    tmp_SS = ( oce.SO.where(cond_mooT,drop=True).isel(time=ll,z=kinf[kk]) * tmpaT + oce.SO.where(cond_mooT,drop=True).isel(time=ll,z=ksup[kk]) * tmpbT ) / tmpxT # Inf if no interpolable value
    tmp_TT = ( oce.THETAO.where(cond_mooT,drop=True).isel(time=ll,z=kinf[kk]) * tmpaT + oce.THETAO.where(cond_mooT,drop=True).isel(time=ll,z=ksup[kk]) * tmpbT ) / tmpxT

    condkk = ( (LEVOF_moor[kk] < epsfr) | (np.isnan(DOMMSK_moor)) )

    tzz = np.squeeze(mp.horizontal_interp( latT, lonT, 1, 1, lat_moor0d, lon_moor0d, tmp_SS, weight=tmp_LEVT[kk,:], skipna=True, filnocvx=True, threshold=epsfr ) )
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    SO_moor[ll,kk] = tzz

    tzz = np.squeeze(mp.horizontal_interp( latT, lonT, 1, 1, lat_moor0d, lon_moor0d, tmp_TT, weight=tmp_LEVT[kk,:], skipna=True, filnocvx=True, threshold=epsfr ) )
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    THETAO_moor[ll,kk] = tzz

#--------------------------------------------------------------------------
# 5b- Create new xarray dataset and save to netcdf

dsmoor = xr.Dataset(
    {
    "so":        (["time", "depth"], np.float32(SO_moor)),
    "thetao":    (["time", "depth"], np.float32(THETAO_moor)),
    "levof":     (["depth"], np.float32(LEVOF_moor)),
    },
    coords={
    "depth": np.float32(dep_moor),
    "time": oce.time.values
    },
)

mp.add_standard_attributes_oce(dsmoor,miss=missval)
put_global_attrs(dsmoor,experiment=exp,avg_hor_res_73S=res_73S,original_sim_name='misomip2_test_case',\
                 original_min_lat=oce.attrs.get('original_minlat'),original_max_lat=oce.attrs.get('original_maxlat'),\
                 original_min_lon=oce.attrs.get('original_minlon'),original_max_lon=oce.attrs.get('original_maxlon') )
dsmoor.attrs['mooring_longitude'] = np.float32(lon_moor0d)
dsmoor.attrs['mooring_latitude'] = np.float32(lat_moor0d)

file_moor = 'OceMoor_'+model+'_'+exp+'.nc'

print('Creating ',file_moor)
dsmoor.to_netcdf(file_moor,unlimited_dims="time")

print('  Execution time: ',datetime.now() - startTime)
