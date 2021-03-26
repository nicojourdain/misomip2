# 2021-03 : Initial code [N. Jourdain, IGE-CNRS]
#============================================================================================
import numpy as np
from scipy import interpolate

#============================================================================================
def vertical_interp(original_depth,interpolated_depth):
   """ Find upper and lower bound indices for simple vertical interpolation

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
def horizontal_interp( lon_in_1d, lat_in_1d, mlat_misomip, mlon_misomip, lon_out_1d, lat_out_1d, var_in_1d ):
   """ Interpolates one-dimension data horizontally to a 2d numpy array reshaped to the misomip standard (lon,lat) format.

       Method: triangular linear barycentryc interpolation, using nans (i.e. gives nan if any nan in the triangle)

       lon\_in\_1d, lat\_in\_1d: 1d longitude and latitude of data to interpolate
 
       var\_in\_1d: 1d input data (same dimension as lon\_in\_1d and lat\_in\_1d)
 
       mlat\_misomip, mlon\_misomip: misomip grid size (nb points) alond latitude and longitude dimensions
 
       lon\_out\_1d, lat\_out\_1d: 1d longitude and latitude of the target misomip grid

   """
   txxxx = interpolate.griddata( (lon_in_1d,lat_in_1d), var_in_1d, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
   out   = np.reshape( txxxx, (mlat_misomip, mlon_misomip) )
   return out

#============================================================================================
def horizontal_interp_nonan( lon_in_1d, lat_in_1d, mlat_misomip, mlon_misomip, lon_out_1d, lat_out_1d, var_in_1d ):
   """ Interpolates one-dimension data horizontally to a 2d numpy array reshaped to the misomip standard (lon,lat) format.

       Method: triangular linear barycentryc interpolation, NOT using nans (i.e. find triangle with non-nan values)
               and nearest-neighbor interpolations for points not surrounded by 3 data points.

       lon\_in\_1d, lat\_in\_1d: 1d longitude and latitude of data to interpolate
 
       var\_in\_1d: 1d input data (same dimension as lon\_in\_1d and lat\_in\_1d)
 
       mlat\_misomip, mlon\_misomip: misomip grid size (nb points) alond latitude and longitude dimensions
 
       lon\_out\_1d, lat\_out\_1d: 1d longitude and latitude of the target misomip grid
   """
   var1d_nonan = var_in_1d[ (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)) ]
   if ( np.size(var1d_nonan)==0 ):
     out = np.zeros((mlat_misomip, mlon_misomip))*np.nan
   else:
     lon_in_1d_nonan = lon_in_1d[ (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)) ]
     lat_in_1d_nonan = lat_in_1d[ (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)) ]
     txxxx = interpolate.griddata( (lon_in_1d_nonan,lat_in_1d_nonan), var1d_nonan, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
     out   = np.reshape( txxxx, (mlat_misomip, mlon_misomip) )
     tssss = interpolate.griddata( (lon_in_1d_nonan,lat_in_1d_nonan), var1d_nonan, (lon_out_1d,lat_out_1d), method='nearest' )
     tmp   = np.reshape( tssss, (mlat_misomip, mlon_misomip) )
     out[ (np.isnan(out)) ] = tmp[ (np.isnan(out)) ]
   return out
