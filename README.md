# misomip2
A python package to postprocess model outputs to standard MISOMIP2 format and to analyse MISOMIP2 multi-model outputs.

This package is largely based on [xarray](http://xarray.pydata.org) and [scipy.interpolate](https://docs.scipy.org/doc/scipy/reference/interpolate.html).

### Contributors
Nicolas C. Jourdain (IGE, CNRS-UGA, Grenoble, France)

### Install
If you don't have a python environment yet, you can install anaconda following [this page](https://docs.anaconda.com/anaconda/install/).

This package may become an anaconda package at some stage, but it is not the case yet, and for now, here is the way to proceed:

If these modules are not installed (check with ```bash conda list```), install them:
```bash
conda install numpy xarray scipy
conda install dask pyproj gsw
```

Then, to enable the import of misomip2 functions from anywhere, do:

```bash
export MYPACK=/User/wmunk/MY_PACKAGES # to be adapted
cd $MYPACK

git clone https://github.com/nicojourdain/misomip2.git
# or git clone git@github.com:nicojourdain/misomip2.git

cat << EOF >> ~/.bashrc # or .bash_profile or .profile or equivalent
export PYTHONPATH="${MYPACK}:\$PYTHONPATH"
EOF
```

<br/><br/>
----------
----------

# Preprocessing
Contains scripts that facilitate interpolation and formatting to the MISOMIP2 standards.

To use the preprocessing tools, start by specifying:
```bash
import misomip2.preproc as mp
```
<br/><br/>

To generate the standard MISOMIP2 [lon,lat,depth] grids, use one of these fucntions:

### misomip2.preproc.grid\_bounds\_oce(region='Amundsen'):
> Gives minimum and maximum longitude and latitude for the common MISOMIP2 grid
>
>    region: 'Amundsen' (default), 'Weddell'
>
_Exemple_: 
```bash
[lonmin,lonmax,latmin,latmax] = grid_bounds_oce(region='Weddell')
```

### misomip2.preproc.generate\_3d\_grid\_oce(region='Amundsen'):
> Generates (longitude, latitude, depth) of the common MISOMIP2 3d grid
> 
>    region: 'Amundsen' (default), 'Weddell'
>
_Exemple_: 
```bash 
[lon,lat,depth]=mp.generate_3d_grid_oce(region='Weddell')
```
<br/>

### misomip2.preproc.generate\_section\_grid\_oce(region='Amundsen'):
> Generates (longitude, latitude, depth) of the common MISOMIP2 section
> 
>    region: 'Amundsen' (default), 'Weddell'
> 
_Exemple_: 
```bash
[lon,lat,depth]=mp.generate_section_grid_oce(region='Weddell')
```
<br/>

### misomip2.preproc.generate\_mooring\_grid\_oce(region='Amundsen'):
> Generates (longitude, latitude, depth) of the common MISOMIP2 mooring
>
>    region: 'Amundsen' (default), 'Weddell'
> 
_Exemple_:
```bash
[lon,lat,depth]=mp.generate_mooring_grid_oce(region='Weddell')
```
<br/><br/>

To put the MISOMIP2 standard attributes to the xarray dataset that will be saved as netcdf:

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
<br/><br/>

To load model outputs, either use an existing function or create a similar one if your model is not covered yet:

### misomip2.preproc.load\_oce\_mod\_mitgcm(files\_in='MITgcm\_output.nc', rho0=1026.0, teos10=False, region='Amundsen' ):
> Read MITgcm outputs and define an xarray dataset containing 
> all variables required in MISOMIP2. It automatically detects
> whether coordinates are stereographic or lon-lat.
>
>    files\_in: file or list of files containing all the variables.
>
>    rho0: reference volumic mass of seawater used in ocean model [kg m-3].
>
>    teos10=False -> assumes the mitgcm outputs are in potential temperature & practical salinity (EOS80).<br/>
>          =True  -> assumes the mitgcm outputs are in CT and AS and convert them to PT and PS.
>
_Example_:
```bash
dir= 'datadir/model/'
ff = [ dir+'MITgcm_y2009.nc', dir+'MITgcm_y2010.nc', dir+'MITgcm_y2011.nc' ]
ds = load_oce_mod_mitgcm(files_in=ff, rho0=1028.0, region='Weddell')
```
<br/>

### misomip2.preproc.load\_oce\_mod\_nemo(file\_mesh\_mask='mesh\_mask.nc', file\_bathy='bathy\_meter.nc', files\_gridT='nemo\_grid\_T.nc', files\_gridU='nemo\_grid\_U.nc', files\_gridV='nemo\_grid\_V.nc', files\_SBC='nemo\_flxT.nc', files\_ice='nemo\_icemod.nc', files\_BSF='nemo\_psi.nc', rho0=1026.0, teos10=False, region='Amundsen' ):
> Read NEMO outputs and define an xarray dataset containing 
> all variables required in MISOMIP2.
>
> Adapted to NEMO's C-grid and standard mesh/mask variables.
>
>    rho0 corresponds to rau0 value in NEMO's eosbn2.F90 i.e. reference volumic mass [kg m-3]
>
>    teos10=False -> assumes the nemo outputs are in potential temperature & practical salinity (EOS80).<br/>
>          =True  -> assumes the nemo outputs are in CT and AS and convert them to PT and PS.
>
>    region = 'Amundsen' (default) or 'Weddell'.
>
>    NB: files\_BSF contains the barotropic streamfunction calculated at U-points, e.g. using
>        the cdfpsi function which is part of the [CDFTOOLS](https://github.com/meom-group/CDFTOOLS).
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
<br/><br/>

The following functions are used for vertical and horizontal interpolation:

### misomip2.preproc.vertical\_interp(original\_depth,interpolated\_depth):
> Find upper and lower bound indices for simple vertical interpolation
>
>    original\_depth: 1d numpy array
>
>    interpolated\_depth: 1d numpy array
>
_Example_:
```bash
[kinf,ksup] = mp.vertical_interp(model_depth,dep_misomip)
```
<br/>

### misomip2.preproc.horizontal\_interp( lon\_in\_1d, lat\_in\_1d, mlat\_misomip, mlon\_misomip, lon\_out\_1d, lat\_out\_1d, var\_in\_1d ):
> Interpolates one-dimension data horizontally to a 2d numpy array reshaped to the misomip standard (lon,lat) format.
>
>    Method: triangular linear barycentryc interpolation, using nans (i.e. gives nan if any nan in the triangle)
>
>    lon\_in\_1d, lat\_in\_1d: 1d longitude and latitude of data to interpolate
>
>    var\_in\_1d: 1d input data (same dimension as lon\_in\_1d and lat\_in\_1d)
>
>    mlat\_misomip, mlon\_misomip: misomip grid size (nb points) alond latitude and longitude dimensions
>
>    lon\_out\_1d, lat\_out\_1d: 1d longitude and latitude of the target misomip grid
>
_Example_:
```bash
VAR_miso = mp.horizontal_interp( ds.lonT, ds.latT, mlat, mlon, lon_miso1d, lat_miso1d, ds.VAR )
```
<br/>

### misomip2.preproc.horizontal\_interp\_nonan( lon\_in\_1d, lat\_in\_1d, mlat\_misomip, mlon\_misomip, lon\_out\_1d, lat\_out\_1d, var\_in\_1d ):
> Interpolates one-dimension data horizontally to a 2d numpy array reshaped to the misomip standard (lon,lat) format.
>
>    Method: triangular linear barycentryc interpolation, NOT USING NANs (i.e. find triangle with non-nan values)
>            and nearest-neighbor interpolations for points not surrounded by 3 data points.
>
>    lon\_in\_1d, lat\_in\_1d: 1d longitude and latitude of data to interpolate
>
>    var\_in\_1d: 1d input data (same dimension as lon\_in\_1d and lat\_in\_1d)
>
>    mlat\_misomip, mlon\_misomip: misomip grid size (nb points) alond latitude and longitude dimensions
>
>    lon\_out\_1d, lat\_out\_1d: 1d longitude and latitude of the target misomip grid
>
_Example_:
```bash
VAR_miso = mp.horizontal_interp_nonan( ds.lonT, ds.latT, mlat, mlon, lon_miso1d, lat_miso1d, ds.VAR )
```
<br/><br/>

----------
----------

# Multi-model Analysis 
Contains scripts to quickly plot multi-model diagnostics. **TO BE COMPLETED**.

To use the analysis tools, start by specifying:
```bash
import misomip2.analysis as ma
```
<br/><br/>

----------
----------
# Examples

People who use the misomip2 package to interpolate their model results to the MISOMIP2 grids are invited to provide their scripts.

### interpolate\_to\_common\_grid.py
> This script has been used to interpolate NEMO and MITgcm outputs to the 3 misomip2 grids (3d grid, section, mooring). To use it, you need to adapt at least section 0 (General information), section 1 (Files and variables) and section 2 (Global attributes of output netcdf). For section1, ou may need to create or modify a load\_oce\_mod\_xxxx function similar to the one existing for NEMO and MITgcm if your model is not covered yet.
