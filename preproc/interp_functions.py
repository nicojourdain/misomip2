# N. Jourdain, IGE-CNRS-UGA, MAR-2021
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
      print([k, interpolated_depth[k], kinf[k], original_depth[kinf[k]], ksup[k], original_depth[ksup[k]]]) 
   return (kinf,ksup)

#============================================================================================
def horizontal_interp( lon_in_1d, lat_in_1d, mlat_misomip, mlon_misomip, lon_out_1d, lat_out_1d, var_in_1d ):
   """ Horizontal interpolation of 1d array (NaNs used for interpolation)
       and reshape to misomip standard (lon,lat) format.

       Method: triangular linear barycentryc interpolation

   """
   txxxx = interpolate.griddata( (lon_in_1d,lat_in_1d), var_in_1d, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
   out   = np.reshape( txxxx, (mlat_misomip, mlon_misomip) )
   return out

#============================================================================================
def horizontal_interp_nonan( lon_in_1d, lat_in_1d, mlat_misomip, mlon_misomip, lon_out_1d, lat_out_1d, var_in_1d ):
   """ Horizontal interpolation of 1d array (NaNs in input data not used for interpolation)
       and reshape to misomip standard (lon,lat) format.     
 
       Method: triangular linear barycentryc interpolation
               and nearest-neighbour where no possible triangulation (on the edges and further).

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
