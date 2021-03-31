# misomip2
A python package to postprocess model outputs to standard MISOMIP2 format and to analyse MISOMIP2 multi-model outputs.

This package is written for python3 and largely based on [xarray](http://xarray.pydata.org) and [scipy.interpolate](https://docs.scipy.org/doc/scipy/reference/interpolate.html).

### Contributors
* Nicolas C. Jourdain (IGE, CNRS-UGA, Grenoble, France)
* Jan De Rydt (U. Northumbria Newcastle, UK)

All MISOMIP participants who want to contribute to this package are invited to fork this github repo and make a pull request as [documented here](https://opensource.com/article/19/7/create-pull-request-github).

### Install
If you don't have a python environment yet, you can install anaconda following [this page](https://docs.anaconda.com/anaconda/install/).

This package may become an anaconda package at some stage, but it is not the case yet, and for now, here is the way to proceed:

If these modules are not installed (check with ```conda list```), install them:
```bash
conda install numpy xarray scipy
conda install dask pyproj
conda install -c conda-forge gsw 
# or conda install -c conda-forge/label/gcc7 gsw 
# or conda install -c conda-forge/label/cf201901 gsw 
# or conda install -c conda-forge/label/cf202003 gsw
```

Then, to clone the misomip2 package and enable the import of misomip2 functions from anywhere, do:

```bash
export MYPACK=/User/wmunk/MY_PACKAGES # to be adapted
cd $MYPACK

git clone https://github.com/nicojourdain/misomip2.git
# or git clone git@github.com:nicojourdain/misomip2.git

cat << EOF >> ~/.bashrc # or .bash_profile or .profile or equivalent
export PYTHONPATH="${MYPACK}:\$PYTHONPATH"
EOF
```

You can update the cloned directory anytime with ```git pull``` executed in that directory. 

<br/><br/>
----------
----------

# Examples

We provide a test case so that users can check that this package works well in their python environment before adapting it to their model outputs. People who use the misomip2 package to interpolate their model results to the MISOMIP2 grids are invited to provide their scripts, in particular if their model is not covered yet.

### Ocean test cases 

We provide 2 months of raw outputs from NEMO and MITGCM (Amundsen Sea configurations) in ```misomip2/examples/models/oce/```. To interpolate these model outputs to the standard MISOMIP2 grids, edit **interpolate_to_common_grid_oce.py** and select ```model='MITGCM_test'``` (you can try ```model='NEMO_test'``` just after). Execute it as:
```bash
python interpolate_to_common_grid_oce.py
```
This should create the following files:
* Oce3d\_MITGCM_test\_A1.nc
* OceMoor\_MITGCM\_test_A1.nc
* OceSec\_MITGCM\_test\_A1.nc

### Adapt to your own ocean configuration

You can copy interpolate_to_common_grid_oce.py and adapt it to your model. You need to adapt at least section 0 (General information), section 1 (Files and variables) and section 2 (Global attributes of output netcdf). For section1, ou may need to create or modify a load\_oce\_mod\_xxxx.py function similar to the one existing for NEMO and MITgcm if your model is not covered yet. If your model is already covered but variable names differ from what is understood by load\_oce\_mod\_xxxx.py, just add options for these variables in load\_oce\_mod\_xxxx.py and [make a pull request](https://opensource.com/article/19/7/create-pull-request-github) to upload it onto the official misomip2 repository (so that you keep this in case of updates).
 
### Ice test cases

**To be completed**

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
> Gives minimum and maximum longitude and latitude for the common MISOMIP2 ocean grid
>
>    region: 'Amundsen' (default), 'Weddell'
>
_Exemple_: 
```bash
[lonmin,lonmax,latmin,latmax] = grid_bounds_oce(region='Weddell')
```

### misomip2.preproc.generate\_3d\_grid\_oce(region='Amundsen'):
> Generates (longitude, latitude, depth) of the common MISOMIP2 3d ocean grid
> 
>    region: 'Amundsen' (default), 'Weddell'
>
_Exemple_: 
```bash 
[lon,lat,depth]=mp.generate_3d_grid_oce(region='Weddell')
```
<br/>

### misomip2.preproc.generate\_section\_grid\_oce(region='Amundsen'):
> Generates (longitude, latitude, depth) of the common MISOMIP2 ocean section
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

### misomip2.preproc.add\_standard\_attributes\_oce(ds,miss=9.969209968386869e36):
> Define standard netcdf attributes for ocean variables that are already present in the ds xarray dataset.
> (these variables must have the MISOMIP2 standard variable names)
> 
>    ds: xarray ocean dataset
>
>    miss: missing value (default=9.969209968386869e36)
>
_Example_:
```bash
mp.add_standard_attributes_oce(dsmiso,miss=1.e20)
```
<br/><br/>

To load model outputs, either use an existing function or create a similar one if your model is not covered yet:

### misomip2.preproc.load\_oce\_mod\_mitgcm(files\_T='MITgcm\_output.nc', files\_S, files\_U, files\_V, files\_I, files\_SRF, files\_M, rho0=1026.0, teos10=False, region='Amundsen' ):
> Read MITgcm outputs and define an xarray dataset containing 
> all variables required in MISOMIP2. It automatically detects
> whether coordinates are stereographic or lon-lat.
>
>    files\_T: file or list of files containing the temperature and related variables [default='MITgcm\_all.nc']
>
>    files\_S: file or list of files containing the salinity variable [optional, default=files\_T]
>
>    files\_U: file or list of files containing the x-velocity and related variables [optional, default=files\_T]
>
>    files\_V: file or list of files containing the y-velocity and related variables [optional, default=files\_T]
>
>    files\_I: file or list of files containing the sea-ice variables [optional, default=files\_T]
>
>    files\_SRF: file or list of files containing the surface fluxes variables [optional, default=files\_T]
>
>    files\_M: file or list of files containing grid/mesh variables [optional, default=files\_T]
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
ds = load_oce_mod_mitgcm(files_T=ff, rho0=1028.0, region='Weddell')
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

