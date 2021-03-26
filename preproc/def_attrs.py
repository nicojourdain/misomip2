# 2021-03 : Initial code [N. Jourdain, IGE-CNRS]

#====================================================================================================
def add_standard_attributes(ds,miss=9.969209968386869e36):
  """Define standard netcdf attributes for variables already present in the ds xarray dataset.

     ds : xarray dataset
     miss : missing value (default=9.969209968386869e36)
  
  """

  ## coordinates :

  if ( "longitude" in ds.coords ):
    ds.longitude.encoding['_FillValue'] = False
    ds.longitude.attrs['units'] = 'degrees_east'
    ds.longitude.attrs['long_name'] = 'longitude'
    ds.longitude.attrs['standard_name'] = 'longitude'

  if ( "latitude" in ds.coords ):
    ds.latitude.encoding['_FillValue'] = False
    ds.latitude.attrs['units'] = 'degrees_north'
    ds.latitude.attrs['long_name'] = 'latitude'
    ds.latitude.attrs['standard_name'] = 'latitude'

  if ( "depth" in ds.coords ):
    ds.depth.encoding['_FillValue'] = False
    ds.depth.attrs['units'] = 'm'
    ds.depth.attrs['positive'] = 'down'
    ds.depth.attrs['long_name'] = 'depth'
    ds.depth.attrs['standard_name'] = 'depth'

  if ( "time" in ds.coords ):
    ds.time.encoding['units'] = 'days since 1900-01-01'
    ds.time.encoding['_FillValue'] = False
    ds.time.attrs['standard_name'] = 'time'


  ## variables :

  if ( "longitude" in ds.data_vars ):
    print('define attributes for variable longitude')
    ds.longitude.attrs['_FillValue'] = miss
    ds.longitude.attrs['units'] = 'degrees_east'
    ds.longitude.attrs['long_name'] = 'longitude'
    ds.longitude.attrs['standard_name'] = 'longitude'

  if ( "latitude" in ds.data_vars ):
    print('define attributes for variable latitude')
    ds.latitude.attrs['_FillValue'] = miss
    ds.latitude.attrs['units'] = 'degrees_north'
    ds.latitude.attrs['long_name'] = 'latitude'
    ds.latitude.attrs['standard_name'] = 'latitude'

  if ( "so" in ds.data_vars ):
     print('define attributes for variable so')
     ds.so.attrs['_FillValue'] = miss
     ds.so.attrs['units'] = '0.001'
     ds.so.attrs['long_name'] = 'Sea Water Salinity (practical salinity)'
     ds.so.attrs['standard_name'] = 'sea_water_salinity' 
     ds.so.attrs['cell_method'] = 'area: mean where sea'

  if ( "thetao" in ds.data_vars ):
     print('define attributes for variable thetao')
     ds.thetao.attrs['_FillValue'] = miss
     ds.thetao.attrs['units'] = 'degC'
     ds.thetao.attrs['long_name'] = 'Sea Water Potential Temperature'
     ds.thetao.attrs['standard_name'] = 'sea_water_potential_temperature'
     ds.thetao.attrs['cell_method'] = 'area: mean where sea'

  if ( "zos" in ds.data_vars ):
     print('define attributes for variable zos')
     ds.zos.attrs['_FillValue'] = miss
     ds.zos.attrs['units'] = 'm'
     ds.zos.attrs['long_name'] = 'Sea Surface Height Above Geoid'
     ds.zos.attrs['standard_name'] = 'sea_surface_height_above_geoid'
     ds.zos.attrs['cell_method'] = 'area: mean where sea'

  #if ( "bigthetao" in ds.data_vars ):
  #   print('define attributes for variable bigthetao')
  #   ds.bigthetao.attrs['_FillValue'] = miss
  #   ds.bigthetao.attrs['units'] = 'degC'
  #   ds.bigthetao.attrs['long_name'] = 'Sea Water Conservative Temperature'
  #   ds.bigthetao.attrs['standard_name'] = 'sea_water_conservative_temperature'
  #   ds.bigthetao.attrs['cell_method'] = 'area:'

  if ( "tob" in ds.data_vars ):
     print('define attributes for variable tob')
     ds.tob.attrs['_FillValue'] = miss
     ds.tob.attrs['units'] = 'degC'
     ds.tob.attrs['long_name'] = 'Sea Water Potential Temperature at Sea Floor'
     ds.tob.attrs['standard_name'] = 'sea_water_potential_temperature_at_sea_floor'
     ds.tob.attrs['cell_method'] = 'area: mean where sea'

  if ( "sob" in ds.data_vars ):
     print('define attributes for variable sob')
     ds.sob.attrs['_FillValue'] = miss
     ds.sob.attrs['units'] = '0.001'
     ds.sob.attrs['long_name'] = 'Sea Water Salinity at Sea Floor (Practical Salinity)'
     ds.sob.attrs['standard_name'] = 'sea_water_salinity_at_sea_floor'
     ds.sob.attrs['cell_method'] = 'area: mean where sea'

  if ( "uo" in ds.data_vars ):
     print('define attributes for variable uo')
     ds.uo.attrs['_FillValue'] = miss
     ds.uo.attrs['units'] = 'm s-1'
     ds.uo.attrs['long_name'] = 'Sea Water X Velocity (Zonal)'
     ds.uo.attrs['standard_name'] = 'sea_water_x_velocity'
     ds.uo.attrs['cell_method'] = 'area: mean where sea'

  if ( "vo" in ds.data_vars ):
     print('define attributes for variable vo')
     ds.vo.attrs['_FillValue'] = miss
     ds.vo.attrs['units'] = 'm s-1'
     ds.vo.attrs['long_name'] = 'Sea Water Y Velocity (Meridional)'
     ds.vo.attrs['standard_name'] = 'sea_water_y_velocity'
     ds.vo.attrs['cell_method'] = 'area: mean where sea'

  #if ( "wfo" in ds.data_vars ):
  #   print('define attributes for variable wfo')
  #   ds.wfo.attrs['_FillValue'] = miss
  #   ds.wfo.attrs['units'] = 'kg m-2 s-1'
  #   ds.wfo.attrs['long_name'] = 'Water Flux into Sea Water (precipitation, river runoff, sea ice, iceberg, ice shelf, restoring)'
  #   ds.wfo.attrs['standard_name'] = 'water_flux_into_sea_water'
  #   ds.wfo.attrs['cell_method'] = 'area: mean where sea'
  #   ds.wfo.attrs['positive'] = 'down'

  if ( "wfoatrli" in ds.data_vars ):
     print('define attributes for variable wfoatrli')
     ds.wfoatrli.attrs['_FillValue'] = miss
     ds.wfoatrli.attrs['units'] = 'kg m-2 s-1'
     ds.wfoatrli.attrs['long_name'] = 'Water Flux into Sea Water from Atmosphere, Rivers, and Land Ice'
     ds.wfoatrli.attrs['standard_name'] = 'water_flux_into_sea_water_from_atmosphere_rivers_land_ice'
     ds.wfoatrli.attrs['cell_method'] = 'area: mean where sea'
     ds.wfoatrli.attrs['positive'] = 'down'

  if ( "wfosicor" in ds.data_vars ):
     print('define attributes for variable wfosicor')
     ds.wfosicor.attrs['_FillValue'] = miss
     ds.wfosicor.attrs['units'] = 'kg m-2 s-1'
     ds.wfosicor.attrs['long_name'] = 'Water Flux into Sea Water Due to Sea Ice Thermodynamics and Correction (restoring)'
     ds.wfosicor.attrs['standard_name'] = 'water_flux_into_sea_water_due_to_sea_ice_thermodynamics_and_correction'
     ds.wfosicor.attrs['cell_method'] = 'area: mean where sea'
     ds.wfosicor.attrs['positive'] = 'down'

  if ( "hfds" in ds.data_vars ):
     print('define attributes for variable hfds')
     ds.hfds.attrs['_FillValue'] = miss
     ds.hfds.attrs['units'] = 'W m-2'
     ds.hfds.attrs['long_name'] = 'Downward Heat Flux at Sea Water Surface (net shortwave and longwave radiative fluxes, sensible and latent heat fluxes)'
     ds.hfds.attrs['standard_name'] = 'surface_downward_heat_flux_in_sea_water'
     ds.hfds.attrs['cell_method'] = 'area: mean where sea'
     ds.hfds.attrs['positive'] = 'down'

  if ( "tauuo" in ds.data_vars ):
     print('define attributes for variable tauuo')
     ds.tauuo.attrs['_FillValue'] = miss
     ds.tauuo.attrs['units'] = 'N m-2'
     ds.tauuo.attrs['long_name'] = 'Sea Water Surface Downward X Stress (stress on the liquid ocean from overlying atmosphere, sea ice, ice shelf, etc.)'
     ds.tauuo.attrs['standard_name'] = 'downward_x_stress_at_sea_water_surface'
     ds.tauuo.attrs['cell_method'] = 'area: mean where sea'
     ds.tauuo.attrs['positive'] = 'down'

  if ( "tauvo" in ds.data_vars ):
     print('define attributes for variable tauvo')
     ds.tauvo.attrs['_FillValue'] = miss
     ds.tauvo.attrs['units'] = 'N m-2'
     ds.tauvo.attrs['long_name'] = 'Sea Water Surface Downward Y Stress (stress on the liquid ocean from overlying atmosphere, sea ice, ice shelf, etc.)'
     ds.tauvo.attrs['standard_name'] = 'downward_y_stress_at_sea_water_surface'
     ds.tauvo.attrs['cell_method'] = 'area: mean where sea'
     ds.tauvo.attrs['positive'] = 'down'

  if ( "msftbarot" in ds.data_vars ):
     print('define attributes for variable vsftbarot')
     ds.msftbarot.attrs['_FillValue'] = miss
     ds.msftbarot.attrs['units'] = 'kg s-1'
     ds.msftbarot.attrs['long_name'] = 'Ocean Barotropic Mass Streamfunction'
     ds.msftbarot.attrs['standard_name'] = ''
     ds.msftbarot.attrs['cell_method'] = 'area: mean where sea'

  if ( "ficeshelf" in ds.data_vars ):
     print('define attributes for variable ficeshelf')
     ds.ficeshelf.attrs['_FillValue'] = miss
     ds.ficeshelf.attrs['units'] = 'kg m-2 s-1'
     ds.ficeshelf.attrs['long_name'] = 'Water Flux into Sea Water from Ice Shelf Basal Melting (negative for refreezing)'
     ds.ficeshelf.attrs['standard_name'] = 'water_flux_into_sea_water_from_iceshelf'
     ds.ficeshelf.attrs['cell_method'] = 'area: mean where sea'

  #if ( "ficeberg" in ds.data_vars ):
  #   print('define attributes for variable ficeberg')
  #   ds.ficeberg.attrs['_FillValue'] = miss
  #   ds.ficeberg.attrs['units'] = 'kg m-2 s-1'
  #   ds.ficeberg.attrs['long_name'] = 'Water Flux into Sea Water from Icebergs'
  #   ds.ficeberg.attrs['standard_name'] = 'water_flux_into_sea_water_from_icebergs'
  #   ds.ficeberg.attrs['cell_method'] = 'area: mean where sea'

  if ( "siconc" in ds.data_vars ):
     print('define attributes for variable siconc')
     ds.siconc.attrs['_FillValue'] = miss
     ds.siconc.attrs['units'] = '%'
     ds.siconc.attrs['long_name'] = 'Sea-Ice Area Percentage'
     ds.siconc.attrs['standard_name'] = 'sea_ice_area_fraction'
     ds.siconc.attrs['cell_method'] = 'area: mean where sea'

  if ( "sivol" in ds.data_vars ):
     print('define attributes for variable sivol')
     ds.sivol.attrs['_FillValue'] = miss
     ds.sivol.attrs['units'] = 'm'
     ds.sivol.attrs['long_name'] = 'Sea-Ice Volume per Area'
     ds.sivol.attrs['standard_name'] = 'sea_ice_thickness'
     ds.sivol.attrs['cell_method'] = 'area: mean where sea'

  if ( "siu" in ds.data_vars ):
     print('define attributes for variable siu')
     ds.siu.attrs['_FillValue'] = miss
     ds.siu.attrs['units'] = 'm s-1'
     ds.siu.attrs['long_name'] = 'X-Component of Sea-Ice Velocity (Zonal)'
     ds.siu.attrs['standard_name'] = 'sea_ice_x_velocity'
     ds.siu.attrs['cell_method'] = 'area: mean where sea-ice'

  if ( "siv" in ds.data_vars ):
     print('define attributes for variable siv')
     ds.siv.attrs['_FillValue'] = miss
     ds.siv.attrs['units'] = 'm s-1'
     ds.siv.attrs['long_name'] = 'Y-Component of Sea-Ice Velocity (Meridional)'
     ds.siv.attrs['standard_name'] = 'sea_ice_y_velocity'
     ds.siv.attrs['cell_method'] = 'area: mean where sea-ice'

  if ( "deptho" in ds.data_vars ):
     print('define attributes for variable deptho')
     ds.deptho.attrs['_FillValue'] = miss
     ds.deptho.attrs['units'] = 'm'
     ds.deptho.attrs['long_name'] = 'Sea Floor Depth Below Geoid' 
     ds.deptho.attrs['standard_name'] = 'sea_floor_depth_below_geoid'
     ds.deptho.attrs['cell_method'] = 'area: mean where sea'

  if ( "depfli" in ds.data_vars ):
     print('define attributes for variable depfli')
     ds.depfli.attrs['_FillValue'] = miss
     ds.depfli.attrs['units'] = 'm'
     ds.depfli.attrs['long_name'] = 'Depth of Floating Ice Base Below Geoid (Ice Shelf Draft)' 
     ds.depfli.attrs['standard_name'] = 'ice_shelf_base_depth_below_geoid'
     ds.depfli.attrs['cell_method'] = 'area: mean where floating ice'

  if ( "levof" in ds.data_vars ):
     print('define attributes for variable levof')
     ds.levof.attrs['_FillValue'] = miss
     ds.levof.attrs['units'] = '%'
     ds.levof.attrs['long_name'] = 'Sea Area Percentage (at each level)'
     ds.levof.attrs['standard_name'] = 'sea_area_fraction'
     ds.levof.attrs['cell_method'] = 'area: cell average'

  if ( "sftfli" in ds.data_vars ):
     print('define attributes for variable sftfli')
     ds.sftfli.attrs['_FillValue'] = miss
     ds.sftfli.attrs['units'] = '%'
     ds.sftfli.attrs['long_name'] = 'Surface Ice Shelf Percentage'
     ds.sftfli.attrs['standard_name'] = 'ice_shelf_fraction'
     ds.sftfli.attrs['cell_method'] = 'area: cell average'

#====================================================================================================
def rename_dimensions(ds):
   """ rename dimensions of xarray dataset ds to standard (x,y,z,time) if needed

       usage: 
              dds = rename_dimensions(dds)
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
   if ( "Z" in ds.dims ):
     ds=ds.rename({'Z':'z'})

   if ( "x_grid_T" in ds.dims ):
     ds=ds.rename({'x_grid_T':'x'})
   if ( "x_grid_U" in ds.dims ):
     ds=ds.rename({'x_grid_U':'x'})
   if ( "x_grid_V" in ds.dims ):
     ds=ds.rename({'x_grid_V':'x'})
   if ( "XC" in ds.dims ):
     ds=ds.rename({'XC':'x'})
   if ( "XG" in ds.dims ):
     ds=ds.rename({'XG':'x'})

   if ( "y_grid_T" in ds.dims ):
     ds=ds.rename({'y_grid_T':'y'})
   if ( "y_grid_U" in ds.dims ):
     ds=ds.rename({'y_grid_U':'y'})
   if ( "y_grid_V" in ds.dims ):
     ds=ds.rename({'y_grid_V':'y'})
   if ( "YC" in ds.dims ):
     ds=ds.rename({'YC':'y'}) 
   if ( "YG" in ds.dims ):
     ds=ds.rename({'YG':'y'})

   return ds
