# 2021-03 : Initial code [N. Jourdain, IGE-CNRS]
#============================================================================================
import numpy as np
from scipy import interpolate

#============================================================================================
def vertical_interp(original_depth,interpolated_depth):
   """ Find upper and lower bound indices for simple vertical 1d interpolation

   """
   if ( original_depth[1] < original_depth[0] ):
     ll_kupward = True
   else:
     ll_kupward = False
   nn = np.size(interpolated_depth)
   kinf=np.zeros(nn,dtype='int')
   ksup=np.zeros(nn,dtype='int')
   for k in np.arange(nn):
      knear = np.argmin( np.abs( original_depth - interpolated_depth[k] ) )
      if (original_depth[knear] > interpolated_depth[k]):
        ksup[k] = knear
        if ll_kupward:
          kinf[k] = np.min([ np.size(original_depth)-1, knear+1 ])
        else:
          kinf[k] = np.max([ 0, knear-1 ])
      else:
        kinf[k] = knear
        if ll_kupward:
          ksup[k] = np.max([ 0, knear-1 ])
        else:
          ksup[k] = np.min([ np.size(original_depth)-1, knear+1 ])
   return (kinf,ksup)

#============================================================================================
def vertical_interp_3d(original_depth,interpolated_depth):
   """ Find upper and lower bound indices for vertical 3d interpolation

       original_depth : 3d xarray data

   """
   if ( original_depth[1,0,0] < original_depth[0,0,0] ):
     ll_kupward = True
   else:
     ll_kupward = False
   nn = np.size(interpolated_depth)
   kinf=np.zeros(nn,dtype='int')
   ksup=np.zeros(nn,dtype='int')
   for k in np.arange(nn):
      diff = original_depth - interpolated_depth[k]
      knear = diff.argmin('Z')
      

      knear = np.argmin( np.abs( original_depth - interpolated_depth[k] ), axis=0 )
      

      if (original_depth[knear] > interpolated_depth[k]):
        ksup[k] = knear
        if ll_kupward:
          kinf[k] = np.min([ np.size(original_depth)-1, knear+1 ])
        else:
          kinf[k] = np.max([ 0, knear-1 ])
      else:
        kinf[k] = knear
        if ll_kupward:
          ksup[k] = np.max([ 0, knear-1 ])
        else:
          ksup[k] = np.min([ np.size(original_depth)-1, knear+1 ])
   return (kinf,ksup)

#============================================================================================
def horizontal_interp( lat_in_1d, lon_in_1d, mlat_misomip, mlon_misomip, lat_out_1d, lon_out_1d, \
                       var_in_1d, weight=[], threshold=1.e20, skipna=False, filnocvx=False ):
   """ Interpolates one-dimension data horizontally to a 2d numpy array reshaped to the misomip standard (lon,lat) format.

       Method: triangular linear barycentryc interpolation, using nans (i.e. gives nan if any nan in the triangle)

       lon\_in\_1d, lat\_in\_1d: 1d longitude and latitude of data to interpolate
 
       mlat\_misomip, mlon\_misomip: misomip grid size (nb points) alond latitude and longitude dimensions
 
       lon\_out\_1d, lat\_out\_1d: 1d longitude and latitude of the target misomip grid

       var\_in\_1d: 1d input data (same dimension as lon\_in\_1d and lat\_in\_1d)

       skipna = False to keep nans in interpolation, i.e. gives nan if any triangle node is nan [default]

              = True to find interpolation triangle nodes with non-nan values

       filnocvx = True to use nearest-neighbor to fill non-convex areas, i.e. for which no triangulation is possible [default]

                = False to fill non-convex areas with nans 

       weight = weights used for interpolation [optional]

       threshold = threshold below which weight value indicates a masked point [default=1.e20]

   """
   miss=-999999.99 # local variable, output missing values will be nan
   var1d_nonan = var_in_1d[ (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)) ]
   if ( np.size(var1d_nonan)==0 ):
     out = np.zeros((mlat_misomip, mlon_misomip))*np.nan
   else:
     if skipna:
       lon_in_1d_nonan = lon_in_1d[ (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)) ]
       lat_in_1d_nonan = lat_in_1d[ (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)) ]
       if ( np.size(weight) == 0 ):
         txxxx = interpolate.griddata( (lon_in_1d_nonan,lat_in_1d_nonan), var1d_nonan, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
       else:
         wgt1d_nonan = weight[ (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)) ]
         wgt1d_nonan[ np.isnan(wgt1d_nonan) | np.isinf(wgt1d_nonan) ] = 0.e0 # if nan in mask but not in input data
         txxxx = interpolate.griddata( (lon_in_1d_nonan,lat_in_1d_nonan), var1d_nonan*wgt1d_nonan, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
         wgt   = interpolate.griddata( (lon_in_1d_nonan,lat_in_1d_nonan), wgt1d_nonan, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
         txxxx = txxxx / wgt
         txxxx[ wgt < threshold ] = miss
       out = np.reshape( txxxx, (mlat_misomip, mlon_misomip) )
       if filnocvx:
         if ( np.size(weight) == 0 ):
           # fill non-convex areas with nearest point (whatever its value):
           tssss = interpolate.griddata( (lon_in_1d_nonan,lat_in_1d_nonan), var1d_nonan, (lon_out_1d,lat_out_1d), method='nearest' )
         else:
           # fill non-convex areas with nearest point having weight >= threshold :
           tmplon=lon_in_1d_nonan[ wgt1d_nonan >= threshold ]
           tmplat=lat_in_1d_nonan[ wgt1d_nonan >= threshold ]
           tmpvar=var1d_nonan[ wgt1d_nonan >= threshold ]
           if ( np.size(tmpvar) == 0 ):
             tssss = np.zeros((mlat_misomip*mlon_misomip)) * np.nan
           else:
             tssss = interpolate.griddata( (tmplon,tmplat), tmpvar, (lon_out_1d,lat_out_1d), method='nearest' )
         tmp   = np.reshape( tssss, (mlat_misomip, mlon_misomip) )
         out[ (np.isnan(out)) | (np.isinf(out)) ] = tmp[ (np.isnan(out)) | (np.isinf(out)) ] # points out of the convex area
         out[ out == miss ] = np.nan # points with weight below threshold
     else:
       if filnocvx:
         if ( np.size(weight) == 0 ):
           txxxx = interpolate.griddata( (lon_in_1d,lat_in_1d), var_in_1d, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
         else:
           txxxx = interpolate.griddata( (lon_in_1d,lat_in_1d), var_in_1d*weight, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
           wgt   = interpolate.griddata( (lon_in_1d,lat_in_1d), weight, (lon_out_1d,lat_out_1d), method='linear', fill_value=1.e0 )
           txxxx = txxxx / wgt
           txxxx[ wgt < threshold ] = miss
         out = np.reshape( txxxx, (mlat_misomip, mlon_misomip) )
         lon_in_1d_nonan = lon_in_1d[ (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)) ]
         lat_in_1d_nonan = lat_in_1d[ (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)) ]  
         if ( np.size(weight) == 0 ):
           # fill non-convex areas with nearest point:
           tssss = interpolate.griddata( (lon_in_1d_nonan,lat_in_1d_nonan), var1d_nonan, (lon_out_1d,lat_out_1d), method='nearest' )
         else:
           # fill non-convex areas with nearest point having weight >= threshold :
           wgt1d_nonan = weight[ (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)) ]
           wgt1d_nonan[ np.isnan(wgt1d_nonan) ] = 0.e0 # if nan in mask but not in input data
           tmplon=lon_in_1d_nonan[ wgt1d_nonan >= threshold ]
           tmplat=lat_in_1d_nonan[ wgt1d_nonan >= threshold ]
           tmpvar=var1d_nonan[ wgt1d_nonan >= threshold ] 
           tmpwgt=wgt1d_nonan[ wgt1d_nonan >= threshold ] 
           tssss = interpolate.griddata( (tmplon,tmplat), tmpvar, (lon_out_1d,lat_out_1d), method='nearest' )
           wgtss = interpolate.griddata( (tmplon,tmplat), tmpwgt, (lon_out_1d,lat_out_1d), method='nearest' )
           tssss[ wgtss < threshold ] = miss
         tmp = np.reshape( tssss, (mlat_misomip, mlon_misomip) )
         out[ (np.isnan(out)) | (np.isinf(out)) ] = tmp[ (np.isnan(out)) | (np.isinf(out)) ]
         out[ out == miss ] = np.nan
       else:
         # Simplest form of horizontal tirangular linear interpolation:
         if ( np.size(weight) == 0 ):
           txxxx = interpolate.griddata( (lon_in_1d,lat_in_1d), var_in_1d, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
         else:
           txxxx = interpolate.griddata( (lon_in_1d,lat_in_1d), var_in_1d*weight, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
           wgt   = interpolate.griddata( (lon_in_1d,lat_in_1d), weight, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
           txxxx = txxxx / wgt
           txxxx[ wgt < threshold ] = np.nan
         out   = np.reshape( txxxx, (mlat_misomip, mlon_misomip) )  
   return out
