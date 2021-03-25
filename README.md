# misomip2
A python package to postprocess model outputs to standard MISOMIP2 format and to analyse MISOMIP2 multi-model outputs.

### Contributors
Nicolas C. Jourdain (IGE, CNRS-UGA, Grenoble, France)

-----

## Preprocessing
Contains scripts that facilitate interpolation and formatting to the MISOMIP2 standards.

To use the preprocessing tools, start by specifying:
```bash
import misomip2.preproc as mp
```

To generate the standard MISOMIP2 [lon,lat,depth] grids, use one of these fucntions:

### misomip2.preproc.generate\_3d\_grid(region='Amundsen'):
> Generates (longitude, latitude, depth) of the common MISOMIP2 3d grid
> 
>    region: 'Amundsen' (default), 'Weddell'
>
_Exemple_: 
```bash 
[lon,lat,depth]=mp.generate_3d_grid(region='Weddell')
```

### misomip2.preproc.generate\_section\_grid(region='Amundsen'):
> Generates (longitude, latitude, depth) of the common MISOMIP2 section
> 
>    region: 'Amundsen' (default), 'Weddell'
> 
_Exemple_: 
```bash
[lon,lat,depth]=mp.generate_section_grid(region='Weddell')
```

### misomip2.preproc.generate\_mooring\_grid(region='Amundsen'):
> Generates (longitude, latitude, depth) of the common MISOMIP2 mooring
>
>    region: 'Amundsen' (default), 'Weddell'
> 
_Exemple_:
```bash
[lon,lat,depth]=mp.generate_mooring_grid(region='Weddell')
```

To put the MISOMIP2 standard attributes to the xarray dataset that will be saves as netcdf:

### add\_standard\_attributes(ds,miss=9.969209968386869e36):
> Define standard netcdf attributes for variables that are already present in the ds xarray dataset.
> (these variables must have the MISOMIP2 standard variable names)
> 
>    ds: xarray dataset
>
>    miss: missing value (default=9.969209968386869e36)
>
_Example_:
```bash
mp.add_standard_attributes(dsmiso,miss=1.e20)
```

-----

## Multi-model Analysis 
Contains scripts to quickly plot multi-model diagnostics.

To use the analysis tools, start by specifying:
```bash
import misomip2.analysis as ma
```

## Examples
