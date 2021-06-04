# 2021-03 : Initial code [N. Jourdain, IGE-CNRS]
#====================================================================================================
import numpy as np
import sys

#====================================================================================================
def grid_bounds_oce(region='Amundsen'):
   """ Gives minimum and maximum longitude and latitude for the common MISOMIP2 ocean grid

       region: 'Amundsen' (default), 'WeddellShelf', 'WeddellGyre'

       exemple: [lonmin,lonmax,latmin,latmax] = grid_bounds_oce(region='Amundsen')
   """
   if ( region == 'Amundsen' ):
     longitude_min = -140.0
     longitude_max =  -90.0
     latitude_min  =  -76.0
     latitude_max  =  -69.0
   elif ( region == 'WeddellShelf' ):
     longitude_min = -90.0
     longitude_max = -15.0
     latitude_min  = -85.0
     latitude_max  = -69.0
   elif ( region == 'WeddellGyre' ):
     longitude_min = -90.0
     longitude_max =  30.0
     latitude_min  = -85.0
     latitude_max  = -60.0
   else:
     sys.exit("~!@#$%^* error : region is not defined, choose either 'Amundsen' or 'Weddell'")

   return [longitude_min,longitude_max,latitude_min,latitude_max]

#====================================================================================================
def generate_3d_grid_oce(region='Amundsen'):
   """Generates (longitude, latitude, depth) of the common MISOMIP2 3d ocean grid

      region: 'Amundsen' (default), 'WeddellShelf', 'WeddellGyre'

      exemple: [lon,lat,depth]=generate_3d_grid_oce(region='Amundsen')
   """

   [lonmin,lonmax,latmin,latmax] = grid_bounds_oce(region=region)

   if ( region == 'Amundsen' ):
     longitude=np.arange(lonmin,lonmax+0.1,0.1)
     latitude=np.arange(latmin,latmax+1./30.,1./30.)
     depth=np.array([0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1500.])
   elif ( region == 'WeddellShelf' ):
     longitude=np.arange(lonmin,lonmax+0.15,0.15)
     latitude=np.arange(latmin,latmax+1./20.,1./20.)
     depth=np.array([0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1500.])
   elif ( region == 'WeddellGyre' ):
     longitude=np.arange(lonmin,lonmax+0.5,0.5)
     latitude=np.arange(latmin,latmax+1./6.,1./6.)
     depth=np.array([0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1500.])
   else:
     sys.exit("~!@#$%^* error : region is not defined, choose either 'Amundsen' or 'WeddellShelf' or 'WeddellGyre'")

   return [longitude,latitude,depth]

#====================================================================================================
def generate_section_grid_oce(region='Amundsen',section=1):
   """Generates (longitude, latitude, depth) of the common MISOMIP2 ocean section

      region: 'Amundsen' (default), 'Weddell'

      section: 1 (default) -> Pine Island Trough for Amundsen
                           -> xxxxxxxxxxxxxxxxxx for Weddell
               2           -> Dotson Trough for Amundsen

      exemple: [lon,lat,depth]=generate_section_grid_oce(region='Amundsen')
   """

   if ( region == 'Amundsen' & section == 1 ):
     longitude=np.array([  -97.643 , -97.994 , -98.344 , -98.694 , -99.045 , -99.395 , -99.746, -100.096,
                          -100.446, -100.797, -101.147, -101.497, -101.847, -102.056, -102.264, -102.473,
                          -102.681, -102.89 , -103.098, -103.306, -103.515, -103.723, -103.932, -104.069,
                          -104.150, -104.231, -104.311, -104.392, -104.472, -104.553, -104.633, -104.714,
                          -104.795, -104.875, -105.005, -105.147, -105.288, -105.430, -105.572, -105.713,
                          -105.855, -105.997, -106.138, -106.280, -106.349, -106.371, -106.393, -106.416,
                          -106.438, -106.46 , -106.483, -106.505, -106.528, -106.550, -106.572, -106.595,
                          -106.617, -106.639, -106.662, -106.684, -106.707, -106.716, -106.687, -106.659,
                          -106.630, -106.602, -106.573, -106.545, -106.516, -106.488, -106.460, -106.431,
                          -106.403, -106.374, -106.330, -106.230, -106.130, -106.030, -105.930, -105.830,
                          -105.730, -105.63 , -105.530, -105.430, -105.330, -105.230, -105.130, -105.030,
                          -104.942, -104.873, -104.803, -104.733, -104.663, -104.593, -104.523, -104.454,
                          -104.384, -104.314, -104.244, -104.174, -104.104, -104.017, -103.929, -103.841,
                          -103.753, -103.665, -103.578, -103.490, -103.402, -103.314, -103.226, -103.138,
                          -103.050, -103.003, -102.963, -102.923, -102.883, -102.843, -102.803, -102.763,
                          -102.724, -102.684, -102.644, -102.604, -102.563, -102.518, -102.472, -102.427,
                          -102.382, -102.338, -102.294, -102.251, -102.208, -102.164, -102.121, -102.104,
                          -102.093, -102.082, -102.071, -102.059, -102.048, -102.037, -102.026, -102.014,
                          -102.003, -101.992, -101.981, -101.969, -101.958, -101.947, -101.936, -101.942,
                          -101.951, -101.96 , -101.969, -101.978, -101.987, -101.996, -102.005, -102.015,
                          -102.024, -102.033, -102.042, -102.051, -102.060 ])
     latitude=np.arange(-75.5+1./30.,-70.0+1./30.,1./30.) 
     depth=np.arange(10.,1510.,10.)
   elif ( region == 'Amundsen' & section == 2 ):
     longitude=np.array([ -114.313, -114.127, -113.94 , -113.753, -113.567, -113.380, -113.193, -113.058,
                          -112.975, -112.892, -112.808, -112.725, -112.642, -112.575, -112.525, -112.475,
                          -112.425, -112.375, -112.325, -112.318, -112.353, -112.389, -112.424, -112.460,
                          -112.495, -112.538, -112.587, -112.635, -112.684, -112.733, -112.781, -112.830,
                          -112.878, -112.927, -112.975, -113.024, -113.079, -113.177, -113.275, -113.373,
                          -113.471, -113.569, -113.667, -113.765, -113.863, -113.961, -114.076, -114.208,
                          -114.340, -114.472, -114.604, -114.735, -114.867, -114.999, -115.123, -115.247,
                          -115.371, -115.495, -115.619, -115.743, -115.867, -115.991, -116.115, -116.239,
                          -116.363, -116.487, -116.580, -116.669, -116.758, -116.847, -116.936, -117.025,
                          -117.114, -117.203, -117.292, -117.381, -117.470, -117.559, -117.648, -117.730,
                          -117.785, -117.840, -117.896, -117.951, -118.006, -118.061, -118.117, -118.172,
                          -118.227, -118.282, -118.338, -118.393, -118.448 ])
     latitude=np.arange(-75.05+1./30.,-71.95+1./30.,1./30.)
     depth=np.arange(10.,1510.,10.)
   elif ( region == 'WeddellShelf' ):
     longitude=np.array([-45.,-46.]) # to update
     latitude=np.array([-80.,-70.]) # to update
     depth=np.arange(10.,1510.,10.) # to update
   elif ( region == 'WeddellGyre' ):
     longitude=np.array([-45.,-46.]) # to update
     latitude=np.array([-80.,-70.]) # to update
     depth=np.arange(10.,1510.,10.) # to update
   else:
     sys.exit("~!@#$%^* error : region is not defined, choose either 'Amundsen' or 'Weddell'")

   if ( np.size(latitude) != np.size(longitude) ):
     sys.exit("~!@#$%^* error : section must be defined with equal longitude and latitude values")

   return [longitude,latitude,depth]

#====================================================================================================
def generate_mooring_grid_oce(region='Amundsen'):
   """Generates (longitude, latitude, depth) of the common MISOMIP2 mooring

      region: 'Amundsen' (default), 'WeddellShelf', 'WeddellGyre'

      exemple: [lon,lat,depth]=generate_mooring_grid_oce(region='Amundsen')
   """

   if ( region == 'Amundsen' ):
     longitude=np.array([ -102. ])
     latitude=np.array([ -75. ])
     depth=np.arange(10.,1210.,10.)
   elif ( region == 'WeddellShelf' ):
     longitude=np.array([ -45. ]) # to update
     latitude=np.array([ -76. ]) # to update
     depth=np.arange(10.,1210.,10.) # to update
   elif ( region == 'WeddellGyre' ):
     longitude=np.array([ -45. ]) # to update
     latitude=np.array([ -76. ]) # to update
     depth=np.arange(10.,1210.,10.) # to update
   else:
     sys.exit("~!@#$%^* error : region is not defined, choose either 'Amundsen' or 'Weddell'")

   return [longitude,latitude,depth]

