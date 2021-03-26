# misomip2
A python package to postprocess model outputs to standard MISOMIP2 format and to analyse MISOMIP2 multi-model outputs.

### Contributors
Nicolas C. Jourdain (IGE, CNRS-UGA, Grenoble, France)

### Install
This may be moved to anaconda, but for now, here is the way to proceed:
```bash
export MYPACK=/User/wmunk/MY_PACKAGES # to be adapted
cd $MYPACK

git clone https://github.com/nicojourdain/misomip2.git
# or git clone git@github.com:nicojourdain/misomip2.git

cat << EOF >> ~/.bashrc # or .bash_profile or .profile or equivalent
export PYTHONPATH="${MYPACK}:\$PYTHONPATH"
EOF
```
Then, the misomip2 fucntions can be imported from anywhere.

----------
----------

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

### misomip2.preproc.add\_standard\_attributes(ds,miss=9.969209968386869e36):
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

To load model outputs, either use an existing function or create a similar one if your model is not covered yet:

### misomip2.preproc.load\_oce\_mod\_mitgcm(files\_in='MITgcm\_output.nc', rho0=1026.0, teos10=False, region='Amundsen' ):
> Read MITgcm outputs and define an xarray dataset containing 
> all variables required in MISOMIP2. It automatically detects
> whether coordinates are stereographic or lon-lat.
>
>    files\_in: file or list of files containing all the variables
>
>    rho0: volumic mass of seawater used in ocean model
>
>    teos10=False -> assumes the nemo outputs are in potential temperature & practical salinity (EOS80)
>
>          =True  -> assumes the nemo outputs are in CT and AS and convert to PT and PS
>
_Example_:
```bash
dir= 'datadir/model/'
ff = [ dir+'MITgcm_y2009.nc', dir+'MITgcm_y2010.nc', dir+'MITgcm_y2011.nc' ]
ds = load_oce_mod_mitgcm(files_in=ff, rho0=1028.0, region='Weddell')
```

### misomip2.preproc.load\_oce\_mod\_nemo(file\_mesh\_mask='mesh\_mask.nc', file\_bathy='bathy\_meter.nc', files\_gridT='nemo\_grid\_T.nc', files\_gridU='nemo\_grid\_U.nc', files\_gridV='nemo\_grid\_V.nc', files\_SBC='nemo\_flxT.nc', files\_ice='nemo\_icemod.nc', files\_BSF='nemo\_psi.nc', rho0=1026.0, teos10=False, region='Amundsen' ):
> Read NEMO outputs and define an xarray dataset containing 
> all variables required in MISOMIP2
>
> Adapted to NEMO's C-grid and standard mesh/mask variables.
>
>    rho0 corresponds to rau0 value in NEMO's eosbn2.F90 i.e. reference volumic mass [kg m-3]
>
>    teos10=False -> assumes the nemo outputs are in potential temperature & practical salinity (EOS80)
>
>          =True  -> assumes the nemo outputs are in CT and AS and convert to PT and PS
>
>    region = 'Amundsen' (default) or 'Weddell'
>
>    NB: files\_BSF contains the barotropic streamfunction calculated at U-points, e.g. using
>        the cdfpsi function which is part of the [CDFTOOLS](https://github.com/meom-group/CDFTOOLS)
>
_Example_:
```bash
dir= 'datadir/model/'
fM = dir+'mesh_mask.nc'
fB = dir+'bathy_meter.nc'
fT = [ dir+'nemo_y2009_grid_T.nc', dir+'nemo_y2010_grid_T.nc' ]
fU = [ dir+'nemo_y2009_grid_U.nc', dir+'nemo_y2010_grid_U.nc' ]
fV = [ dir+'nemo_y2009_grid_V.nc', dir+'nemo_y2010_grid_V.nc' ]
fS = [ dir+'nemo_y2009_SBC.nc', dir+'nemo_y2010_SBC.nc' ]
fI = [ dir+'nemo_y2009_icemod.nc', dir+'nemo_y2010_icemod.nc' ]
fP = [ dir+'nemo_y2009_psi.nc', dir+'nemo_y2010_psi.nc' ]
ds = load_oce_mod_nemo(files_gridT=fT,files_gridU=fU,files_gridV=fV,files_SBC=fS,files_ice=fI,files_BSF=fP,rho0=1028.0,teos10=True,region='Weddell')
```


----------
----------

## Multi-model Analysis 
Contains scripts to quickly plot multi-model diagnostics.

To use the analysis tools, start by specifying:
```bash
import misomip2.analysis as ma
```

## Examples
