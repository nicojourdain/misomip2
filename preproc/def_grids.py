# 2021-03 : Initial code [N. Jourdain, IGE-CNRS]
#====================================================================================================
import numpy as np
import sys

#====================================================================================================
def grid_bounds_oce(region='Amundsen'):
   """ Gives minimum and maximum longitude and latitude for the common MISOMIP2 ocean grid

       region: 'Amundsen' (default), 'Weddell'

       exemple: [lonmin,lonmax,latmin,latmax] = grid_bounds_oce(region='Amundsen')
   """
   if ( region == 'Amundsen' ):
     longitude_min = -140.0
     longitude_max = -90.0
     latitude_min = -76.0
     latitude_max = -69.0
   elif ( region == 'Weddell' ):
     longitude_min = -90.0
     longitude_max = 0.0
     latitude_min = -85.0
     latitude_max = -68.9
   else:
     sys.exit("~!@#$%^* error : region is not defined, choose either 'Amundsen' or 'Weddell'")

   return [longitude_min,longitude_max,latitude_min,latitude_max]

#====================================================================================================
def generate_3d_grid_oce(region='Amundsen'):
   """Generates (longitude, latitude, depth) of the common MISOMIP2 3d ocean grid

      region: 'Amundsen' (default), 'Weddell'

      exemple: [lon,lat,depth]=generate_3d_grid_oce(region='Amundsen')
   """

   [lonmin,lonmax,latmin,latmax] = grid_bounds_oce(region=region)

   if ( region == 'Amundsen' ):
     longitude=np.arange(lonmin,lonmax+0.1,0.1)
     latitude=np.arange(latmin,latmax+1./30.,1./30.)
     depth=np.array([0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1500.])
   elif ( region == 'Weddell' ):
     longitude=np.arange(lonmin,lonmax+0.1,0.1)
     latitude=np.arange(latmin,latmax+1./3.,1./3.)
     depth=np.array([0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1500.])
   else:
     sys.exit("~!@#$%^* error : region is not defined, choose either 'Amundsen' or 'Weddell'")

   return [longitude,latitude,depth]

#====================================================================================================
def generate_section_grid_oce(region='Amundsen'):
   """Generates (longitude, latitude, depth) of the common MISOMIP2 ocean section

      region: 'Amundsen' (default), 'Weddell'

      exemple: [lon,lat,depth]=generate_section_grid_oce(region='Amundsen')
   """

   if ( region == 'Amundsen' ):
     longitude=np.array([ -99.345382557121,  -99.553861837121,  -99.762341117121,  -99.970820397121, -100.179299677121, -100.387778957121,
                         -100.596258237121, -100.804737517121, -101.013216797121, -101.221696077121, -101.430175357121, -101.638654637121,
                         -101.847133917121, -102.055613197435, -102.264092477749, -102.472571758063, -102.681051038377, -102.889530318691, 
                         -103.098009599005, -103.306488879319, -103.514968159633, -103.723447439946, -103.931926720260, -104.069419090148, 
                         -104.149994027351, -104.230568964555, -104.311143901758, -104.391718838962, -104.472293776165, -104.552868713369, 
                         -104.633443650572, -104.714018587776, -104.794593524979, -104.875168462183, -105.004979380028, -105.146641141848, 
                         -105.288302903669, -105.429964665490, -105.571626427310, -105.713288189131, -105.854949950952, -105.996611712772, 
                         -106.138273474593, -106.279935236414, -106.348594068822, -106.370969998271, -106.393345927719, -106.415721857168,
                         -106.438097786616, -106.460473716064, -106.482849645513, -106.505225574961, -106.527601504409, -106.549977433858,
                         -106.572353363306, -106.594729292755, -106.617105222203, -106.639481151651, -106.661857081100, -106.684233010548,
                         -106.706608939997, -106.715616260163, -106.687160975610, -106.658705691057, -106.630250406504, -106.601795121951,
                         -106.573339837398, -106.544884552846, -106.516429268293, -106.487973983740, -106.459518699187, -106.431063414634,
                         -106.402608130081, -106.374152845528, -106.329601630080, -106.229608874882, -106.129616119684, -106.029623364486,
                         -105.929630609288, -105.829637854090, -105.729645098892, -105.629652343693, -105.529659588495, -105.429666833297,
                         -105.329674078099, -105.229681322901, -105.129688567703, -105.029695812505, -104.942426822603, -104.872585170002,
                         -104.802743517402, -104.732901864802, -104.663060212202, -104.593218559601, -104.523376907001, -104.453535254401,
                         -104.383693601800, -104.313851949200, -104.244010296600, -104.174168644000, -104.104326991399, -104.016812200543,
                         -103.928955116490, -103.841098032436, -103.753240948382, -103.665383864329, -103.577526780275, -103.489669696221,
                         -103.401812612168, -103.313955528114, -103.226098444060, -103.138241360007, -103.050384275953, -103.002564445588,
                         -102.962713660921, -102.922862876254, -102.883012091587, -102.843161306920, -102.803310522254, -102.763459737587,
                         -102.723608952920, -102.683758168253, -102.643907383586, -102.604056598919, -102.563241914847, -102.517819457313,
                         -102.472396999779, -102.426974542246, -102.381552084712, -102.337500808811, -102.294220822578, -102.250940836345,
                         -102.207660850112, -102.164380863879, -102.121100877646, -102.104378099355, -102.093133345862, -102.081888592370,
                         -102.070643838877, -102.059399085385, -102.048154331893, -102.036909578400, -102.025664824908, -102.014420071415,
                         -102.003175317923, -101.991930564430, -101.980685810938, -101.969441057445, -101.958196303953, -101.946951550460,
                         -101.935706796968, -101.942049707602, -101.951114035088, -101.960178362573, -101.969242690058, -101.978307017544,
                         -101.987371345029, -101.996435672515, -102.005500000000, -102.014564327485, -102.023628654971, -102.032692982456,
                         -102.041757309942, -102.050821637427, -102.059885964912 ])
     latitude=np.arange(-75.5+1./30.,-70.0+1./30.,1./30.) 
     depth=np.arange(10.,1510.,10.)
   elif ( region == 'Weddell' ):
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

      region: 'Amundsen' (default), 'Weddell'

      exemple: [lon,lat,depth]=generate_mooring_grid_oce(region='Amundsen')
   """

   if ( region == 'Amundsen' ):
     longitude=np.array([ -102. ])
     latitude=np.array([ -75. ])
     depth=np.arange(10.,1210.,10.)
   elif ( region == 'Weddell' ):
     longitude=np.array([ -45. ]) # to update
     latitude=np.array([ -76. ]) # to update
     depth=np.arange(10.,1210.,10.) # to update
   else:
     sys.exit("~!@#$%^* error : region is not defined, choose either 'Amundsen' or 'Weddell'")

   return [longitude,latitude,depth]

