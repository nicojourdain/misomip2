###################################################################################
# Script to interpolate native ocean model outputs to the standard MISOMIP2 grid.
#
# Works for :
#          * A,B or C-grid (only tested with C-grid).
#          * non-regular grid such as NEMO's eORCA.
#
# Does not work for :
#          * unstructured grids (needs some modifications but most of it should work)
#
# History: 
#          * initial version and tests for NEMO outputs [N. Jourdain, March 2021]
#
###################################################################################
import numpy as np
import xarray as xr
import sys
import misomip2.preproc as misopre
from scipy import interpolate
import os

# Original NEMO config & case :
config='AMUXL025'
case='GNJ004_BM02'

# Official name in MISOMIP2:
model='NEMO3.6_IGE-CNRS-UGA_AMUXL025a'
#model='MITGCM_UNN_AMUNDSEN'

reg='Amundsen' # 'Amundsen' or 'Weddell'

exp='A1' # MISOMIP2 experiment ('A1', 'W1', 'A2', ...)

# Time period :
yeari=2010
yearf=2012

data_dir='/Users/jourdain/MISOMIP2/MODEL_OUTPUTS/'+model+'/ORIGINAL'

missval=9.969209968386869e36

#--------------------------------------------------------------------------
# 1- Files and variables

# loading an xarray dataset containing all useful variables with (x,y) reshaped 
# as 1 dimension in case of original structured grid :

if ( model[0:4] == 'NEMO' ):

   # Original grid coordinates and mask (mesh_mask) :
   f_mesh=data_dir+'/mesh_mask_AMUXL025_BedMachineAntarctica-2019-05-24.nc'
   # Bathymetry :
   f_bathy=data_dir+'/bathy_meter_AMUXL025_BedMachineAntarctica-2019-05-24.nc'
   # NEMO outputs :
   fs_gridT = [data_dir+'/'+config+'-'+case+'_1m_'+year.astype('str')+'0101_'+year.astype('str')+'1231_grid_T.nc' for year in np.arange(yeari,yearf+1)]
   fs_gridU = [data_dir+'/'+config+'-'+case+'_1m_'+year.astype('str')+'0101_'+year.astype('str')+'1231_grid_U.nc' for year in np.arange(yeari,yearf+1)]
   fs_gridV = [data_dir+'/'+config+'-'+case+'_1m_'+year.astype('str')+'0101_'+year.astype('str')+'1231_grid_V.nc' for year in np.arange(yeari,yearf+1)]
   fs_SBC   = [data_dir+'/'+config+'-'+case+'_1m_'+year.astype('str')+'0101_'+year.astype('str')+'1231_SBC.nc' for year in np.arange(yeari,yearf+1)]
   fs_ice   = [data_dir+'/'+config+'-'+case+'_1m_'+year.astype('str')+'0101_'+year.astype('str')+'1231_icemod.nc' for year in np.arange(yeari,yearf+1)]
   # Barotropic Streamfunction calculated at U points using the cdfpsi function which is part of the cdftools (https://github.com/meom-group/CDFTOOLS):
   fs_BSF   = [data_dir+'/'+config+'-'+case+'_1m_'+year.astype('str')+'0101_'+year.astype('str')+'1231_psi.nc' for year in np.arange(yeari,yearf+1)]

   print('LOADING NEMO...')
   oce = misopre.load_oce_mod_nemo( file_mesh_mask=f_mesh, file_bathy=f_bathy, files_gridT=fs_gridT,\
                  files_gridU=fs_gridU, files_gridV=fs_gridV, files_SBC=fs_SBC, files_ice=fs_ice,\
                  files_BSF=fs_BSF, rho0=1026.0, teos10=False )

elif ( model[0:10] == 'MITGCM_UNN' ):

   print('LOADING MITGCM...')
   fs = data_dir+'MITgcm_output_example.nc'
   oce = misopre.load_oce_mod_mitgcm( files_in=fs,\
                                      rho0=1026.0, teos10=False )

else:

   sys.exit('Unknown model ==> Write a function to load this model outputs')

print(oce)

#--------------------------------------------------------------------------
# 2- Local function to put the apropriate global attributes in output netcdf :

def put_global_attrs(ds,experiment='TBD',avg_hor_res_73S=0.0,original_sim_name='None',\
                     original_min_lat=-90.0,original_max_lat=90.0,original_min_lon=-180.0,original_max_lon=180.0):
   """ Put global attributes to the ds xarray dataset

   """
   ds.attrs['project'] = 'MISOMIP2'
   ds.attrs['contact'] = 'Nicolas Jourdain <nicolas.jourdain@univ-grenoble-alpes.fr>' # name <email>
   ds.attrs['institute'] = 'CNRS-UGA-IGE'
   ds.attrs['computing_facility'] = 'occigen-CINES'                    # Computing center where the simulation was run
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

maskTnan = oce.maskT.where( (oce.maskT==1) ) # 3d mask (=1 if ocean, =nan elsewhere)
maskUnan = oce.maskU.where( (oce.maskU==1) )
maskVnan = oce.maskV.where( (oce.maskV==1) )
maskTnan2d = oce.maskT.max(dim='z').where( oce.maskT.max(dim='z')==1 ) # 2d mask (=1 for open ocean and cavities, =nan elsewhere)
maskUnan2d = oce.maskU.max(dim='z').where( oce.maskU.max(dim='z')==1 )
maskVnan2d = oce.maskV.max(dim='z').where( oce.maskU.max(dim='z')==1 )

[lon_miso,lat_miso,dep_miso] = misopre.generate_3d_grid(region=reg)

mlon = np.size(lon_miso)
mlat = np.size(lat_miso)
mdep = np.size(dep_miso)

lon2d_miso, lat2d_miso = np.meshgrid( lon_miso, lat_miso )

[kinf,ksup] = misopre.vertical_interp(oce.depTUV.values,dep_miso)

mtime = np.shape(oce.SO)[0]

# useful quantities for interpolation and rotation of velocity vectors:
lonT1d = np.reshape( oce.lonT.values, np.size(oce.lonT) )
latT1d = np.reshape( oce.latT.values, np.size(oce.latT) )    
lonU1d = np.reshape( oce.lonU.values, np.size(oce.lonU) )
latU1d = np.reshape( oce.latU.values, np.size(oce.latU) )
lonV1d = np.reshape( oce.lonV.values, np.size(oce.lonV) )
latV1d = np.reshape( oce.latV.values, np.size(oce.latV) )
lon_miso1d = np.reshape( lon2d_miso, mlon*mlat )
lat_miso1d = np.reshape( lat2d_miso, mlon*mlat )

# mask of original domain (nan where interpolation of any of T, U, V grid is nan):
DOMMSK_miso =               misopre.horizontal_interp( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.DOMMSKT )
DOMMSK_miso = DOMMSK_miso + misopre.horizontal_interp( lonU1d, latU1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.DOMMSKU )
DOMMSK_miso = DOMMSK_miso + misopre.horizontal_interp( lonV1d, latV1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.DOMMSKV )

# vertical then horizontal interpolation of constant 3d fields to common grid :
LEVOF_miso = np.zeros((mdep,mlat,mlon)) + missval
for kk in np.arange(mdep):
  if ( kinf[kk] == ksup[kk] ):
    tmpaT = 1.e0
    tmpbT = 0.e0
  else:
    tmpaT = oce.depTUV.isel(z=ksup[kk]) - dep_miso[kk]
    tmpbT = dep_miso[kk] - oce.depTUV.isel(z=kinf[kk])
  tmpxT = tmpaT + tmpbT
  tmp_OC = ( oce.LEVOF.isel(z=kinf[kk]) * tmpaT + oce.LEVOF.isel(z=ksup[kk]) * tmpbT ) / tmpxT
  tzz = misopre.horizontal_interp( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, tmp_OC )
  tzz[ np.isnan(DOMMSK_miso) ] = missval
  LEVOF_miso[kk,:,:] = tzz

# horizontal interpolation of constant 2d fields to common horizontal grid :
theT = misopre.horizontal_interp( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.thetaT )
theU = misopre.horizontal_interp( lonU1d, latU1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.thetaU )
theV = misopre.horizontal_interp( lonV1d, latV1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.thetaV )

SFTFLI_miso = misopre.horizontal_interp( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.SFTFLI )
SFTFLI_miso[ np.isnan(DOMMSK_miso) ] = missval

DEPFLI = oce.DEPFLI.where( (oce.SFTFLI.values > 1.), np.nan ) # to avoid interpolation with 0 depth beyond the grounding line
DEPFLI_miso = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, DEPFLI )
DEPFLI_miso[ (SFTFLI_miso < 1.) | np.isnan(DEPFLI_miso) ] = 0.e0
DEPFLI_miso[ np.isnan(DOMMSK_miso) ] = missval

DEPTHO = oce.DEPTHO.where( (oce.SFTFLI.values > 1.) | (oce.LEVOF[0,:].values > 1.), np.nan )
DEPTHO_miso = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, DEPTHO )
DEPTHO_miso[ ((SFTFLI_miso < 1.) & (LEVOF_miso[0,:,:] < 1.)) | np.isnan(DEPTHO_miso) ] = 0.e0
DEPTHO_miso[ np.isnan(DOMMSK_miso) ] = missval

# interpolation of time-varying fields:

SO_miso     = np.zeros((mtime,mdep,mlat,mlon)) + missval
THETAO_miso = np.zeros((mtime,mdep,mlat,mlon)) + missval
UO_miso     = np.zeros((mtime,mdep,mlat,mlon)) + missval
VO_miso     = np.zeros((mtime,mdep,mlat,mlon)) + missval
ZOS_miso    = np.zeros((mtime,mlat,mlon)) + missval
TOB_miso    = np.zeros((mtime,mlat,mlon)) + missval
SOB_miso    = np.zeros((mtime,mlat,mlon)) + missval
FICESHELF_miso = np.zeros((mtime,mlat,mlon)) + missval
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

LEVOF_maxdep = np.amax(LEVOF_miso,axis=0)

# masking with nan to not consider points in interpolation 
# (instead of defining a mask for the bottom layer):
TOB=oce.TOB.where( (np.abs(oce.SOB.values) > 1.e-3) )
SOB=oce.SOB.where( (~np.isnan(TOB.values)) )

for ll in np.arange(2):
#for ll in np.arange(mtime):

  tzz = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.ZOS.isel(time=ll)*maskTnan.isel(z=0) )
  tzz[ (LEVOF_maxdep < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; ZOS_miso[ll,:,:] = tzz

  tzz = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, TOB.isel(time=ll) )
  tzz[ (LEVOF_maxdep < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; TOB_miso[ll,:,:] = tzz

  tzz = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, SOB.isel(time=ll) )
  tzz[ (LEVOF_maxdep < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; SOB_miso[ll,:,:] = tzz

  tzz = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.FICESHELF.isel(time=ll)*maskTnan.max('z',skipna=True) )
  tzz[ (SFTFLI_miso < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; FICESHELF_miso[ll,:,:] = tzz

  tzz = misopre.horizontal_interp_nonan( lonU1d, latU1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.MSFTBAROT.isel(time=ll)*maskUnan.isel(z=0) )
  tzz[ (LEVOF_maxdep < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; MSFTBAROT_miso[ll,:,:] = tzz

  tzz = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.HFDS.isel(time=ll)*maskTnan.isel(z=0) )
  tzz[ (LEVOF_maxdep < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; HFDS_miso[ll,:,:] = tzz

  tzz = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.WFOATRLI.isel(time=ll)*maskTnan.isel(z=0) )
  tzz[ (LEVOF_maxdep < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; WFOATRLI_miso[ll,:,:] = tzz

  tzz = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.WFOSICOR.isel(time=ll)*maskTnan.isel(z=0) )
  tzz[ (LEVOF_miso[0,:,:] < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; WFOSICOR_miso[ll,:,:] = tzz

  tzz = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.SICONC.isel(time=ll)*maskTnan.isel(z=0) )
  tzz[ (LEVOF_miso[0,:,:] < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; SICONC_miso[ll,:,:] = tzz

  tzz = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.SIVOL.isel(time=ll)*maskTnan.isel(z=0) )
  tzz[ (LEVOF_miso[0,:,:] < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; SIVOL_miso[ll,:,:] = tzz

  # sea-ice velocities: rotation and interpolation weighted by sea-ice concentration
  SIUX_notrot = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.SIUX.isel(time=ll)*oce.SICONC.isel(time=ll) )
  SIUX_notrot[ (LEVOF_miso[0,:,:] < 1.) ] = np.nan
  SIVY_notrot = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.SIVY.isel(time=ll)*oce.SICONC.isel(time=ll) )
  SIVY_notrot[ (LEVOF_miso[0,:,:] < 1.) ] = np.nan
  SIU_miso[ll,:,:] = ( SIUX_notrot * np.cos(theT) + SIVY_notrot * np.sin(theT) ) / SICONC_miso[ll,:,:] # rotated to zonal
  SIV_miso[ll,:,:] = ( SIVY_notrot * np.cos(theT) - SIUX_notrot * np.sin(theT) ) / SICONC_miso[ll,:,:] # rotated to meridional
  SIU_miso[ ( np.isnan(SIU_miso) ) | ( np.isinf(SIU_miso) ) | (SICONC_miso[ll,:,:] < 1.) | (np.isnan(DOMMSK_miso)) ] = missval
  SIV_miso[ ( np.isnan(SIV_miso) ) | ( np.isinf(SIV_miso) ) | (SICONC_miso[ll,:,:] < 1.) | (np.isnan(DOMMSK_miso)) ] = missval

  # wind stress: rotation and interpolation
  TAUX_notrot = misopre.horizontal_interp_nonan( lonU1d, latU1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.TAUX.isel(time=ll)*maskUnan2d )
  TAUX_notrot[ (LEVOF_miso[0,:,:] < 1.) & (SFTFLI_miso < 1.) ] = np.nan
  TAUY_notrot = misopre.horizontal_interp_nonan( lonV1d, latV1d, mlat, mlon, lon_miso1d, lat_miso1d, oce.TAUY.isel(time=ll)*maskVnan2d )
  TAUY_notrot[ (LEVOF_miso[0,:,:] < 1.) & (SFTFLI_miso < 1.) ] = np.nan
  TAUUO_miso[ll,:,:] = TAUX_notrot * np.cos(theU) + TAUY_notrot * np.sin(theV) # rotated to zonal
  TAUVO_miso[ll,:,:] = TAUY_notrot * np.cos(theV) - TAUX_notrot * np.sin(theU) # rotated to meridional
  TAUUO_miso[ ( np.isnan(TAUUO_miso) ) | ( np.isinf(TAUUO_miso) ) | (np.isnan(DOMMSK_miso)) ] = missval
  TAUVO_miso[ ( np.isnan(TAUVO_miso) ) | ( np.isinf(TAUVO_miso) ) | (np.isnan(DOMMSK_miso)) ] = missval

  for kk in np.arange(mdep):

    # vertical interpolation to common vertical grid :
    if ( kinf[kk] == ksup[kk] ):
      tmpaT = oce.maskT.isel(z=kinf[kk])
      tmpbT = tmpaT*0
    else:
      tmpaT = oce.maskT.isel(z=kinf[kk]) * (oce.depTUV.isel(z=ksup[kk])-dep_miso[kk])
      tmpbT = oce.maskT.isel(z=ksup[kk]) * (dep_miso[kk]-oce.depTUV.isel(z=kinf[kk]))
    tmpxT = tmpaT + tmpbT
    tmp_SS = ( oce.SO.isel(time=ll,z=kinf[kk]) * tmpaT + oce.SO.isel(time=ll,z=ksup[kk]) * tmpbT ) / tmpxT # Inf if no interpolable value
    tmp_TT = ( oce.THETAO.isel(time=ll,z=kinf[kk]) * tmpaT + oce.THETAO.isel(time=ll,z=ksup[kk]) * tmpbT ) / tmpxT

    if ( kinf[kk] == ksup[kk] ):
      tmpaU = oce.maskU.isel(z=kinf[kk])
      tmpbU = tmpaU*0
    else:
      tmpaU = oce.maskU.isel(z=kinf[kk]) * (oce.depTUV.isel(z=ksup[kk])-dep_miso[kk])
      tmpbU = oce.maskU.isel(z=ksup[kk]) * (dep_miso[kk]-oce.depTUV.isel(z=kinf[kk]))
    tmpxU = tmpaU + tmpbU
    tmp_UX = ( oce.UX.isel(time=ll,z=kinf[kk]) * tmpaU + oce.UX.isel(time=ll,z=ksup[kk]) * tmpbU ) / tmpxU # Inf if no interpolable value 
 
    if ( kinf[kk] == ksup[kk] ):
      tmpaV = oce.maskV.isel(z=kinf[kk])
      tmpbV = tmpaV*0
    else:
      tmpaV = oce.maskV.isel(z=kinf[kk]) * (oce.depTUV.isel(z=ksup[kk])-dep_miso[kk])
      tmpbV = oce.maskV.isel(z=ksup[kk]) * (dep_miso[kk]-oce.depTUV.isel(z=kinf[kk]))
    tmpxV = tmpaV + tmpbV
    tmp_VY = ( oce.VY.isel(time=ll,z=kinf[kk]) * tmpaV + oce.VY.isel(time=ll,z=ksup[kk]) * tmpbV ) / tmpxV # Inf if no interpolable value

    # horizontal interpolations of time-varying 3d fields to common horizontal grid : 
    tzz = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, tmp_SS )
    tzz[ (LEVOF_miso[kk,:,:] < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; SO_miso[ll,kk,:,:] = tzz

    tzz = misopre.horizontal_interp_nonan( lonT1d, latT1d, mlat, mlon, lon_miso1d, lat_miso1d, tmp_TT )
    tzz[ (LEVOF_miso[kk,:,:] < 1.) | (np.isnan(DOMMSK_miso)) ] = missval ; THETAO_miso[ll,kk,:,:] = tzz

    UX_notrot = misopre.horizontal_interp_nonan( lonU1d, latU1d, mlat, mlon, lon_miso1d, lat_miso1d, tmp_UX )
    VY_notrot = misopre.horizontal_interp_nonan( lonV1d, latV1d, mlat, mlon, lon_miso1d, lat_miso1d, tmp_VY )
    tzz = UX_notrot * np.cos(theU) + VY_notrot * np.sin(theV) # rotated to zonal
    tzz[ (LEVOF_miso[kk,:,:] < 1.) | (np.isnan(tzz)) | (np.isnan(DOMMSK_miso)) ] = missval
    UO_miso[ll,kk,:,:] = tzz
    tzz = VY_notrot * np.cos(theV) - UX_notrot * np.sin(theU) # rotated to meridional
    tzz[ (LEVOF_miso[kk,:,:] < 1.) | (np.isnan(tzz)) | (np.isnan(DOMMSK_miso)) ] = missval
    VO_miso[ll,kk,:,:] = tzz
 
# Some quantities needed to define the global attributes:
res_73S = 0.5 * (  oce.dxT.where(  (oce.latT < -72.5) & (oce.latT >= -73.5) & (oce.lonT >= lon_miso.min()) & (oce.lonT <= lon_miso.max()) \
                                 & (oce.maskT.isel(z=0)==1)   ).mean().values \
                 + oce.dyT.where(  (oce.latT < -72.5) & (oce.latT >= -73.5) & (oce.lonT >= lon_miso.min()) & (oce.lonT <= lon_miso.max()) \
                                 & (oce.maskT.isel(z=0)==1) ).mean().values )
print('Average horizontal resolution at 73S :',res_73S)


#--------------------------------------------------------------------------
# 3b- Create new xarray dataset and save to netcdf

dsmiso = xr.Dataset(
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

misopre.add_standard_attributes(dsmiso,miss=missval)
put_global_attrs(dsmiso,experiment=exp,avg_hor_res_73S=res_73S,original_sim_name='nemo_'+config+'_'+case,\
                 original_min_lat=oce.attrs.get('original_minlat'),original_max_lat=oce.attrs.get('original_maxlat'),\
                 original_min_lon=oce.attrs.get('original_minlon'),original_max_lon=oce.attrs.get('original_maxlon') )
print(dsmiso)
dsmiso.to_netcdf('test.nc',unlimited_dims="time")
