# misomip2
A python package to postprocess model outputs to standard MISOMIP2 format and to analyse MISOMIP2 multi-model outputs.

This package is written for python3 and largely based on [xarray](http://xarray.pydata.org) and [scipy.interpolate](https://docs.scipy.org/doc/scipy/reference/interpolate.html).

This package contains:
* **preproc** : contains functions used to pre-process model outputs, i.e. to interpolate them to the MISOMIP2 grids and write files with standard attributes.
* **examples** : contains scripts that can be used on provided test cases (e.g. to check that your python environment works) and that can be adapted to your specific output files.
* **analysis** : (TO BE COMPLETED) contains functions used to analyse multiple model outputs that were previously interpolated to the common grids.

### Contributors
* Nicolas C. Jourdain (IGE, CNRS-UGA, Grenoble, France)
* Jan De Rydt (University Northumbria Newcastle, UK)
* Yoshihiro Nakayama (Hokkaido University, Japan)
* Ole Richter (Alfred Wegener Institute, Germany)

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

Then, clone the misomip2 package :

```bash
git clone https://github.com/nicojourdain/misomip2.git
# or git clone git@github.com:nicojourdain/misomip2.git
```

You can update the cloned directory anytime with ```git pull``` executed in that directory. 

<br/><br/>
----------
----------

# Examples of ocean interpolation to the MISOMIP2 grids

We provide a few test cases so that users can check that this package works well in their python environment before adapting it to their model outputs. These test cases are also used by the developers of this package to check that it is working with various types of model. 

For simplicity, we describe the example directly used in the misomip2 repositories, but the interpolation script can be run from anywhere.

```bash
cd misomip2/examples
vi interpolate_to_common_grid_oce.py # or using any other text editor than vi
```
The first thing to change in this file is the directory in which the misomip2 package has been cloned, i.e. change this path:
```python
sys.path.append("/Users/jourdain/MY_PACKAGES")
```
The second thing is to choose the test case you want to run by uncommenting one of these lines:
```python
test_case='NEMO_test'
#test_case='MITGCM_test'
#test_case='ROMS_test'
#test_case='eORCA025_test'
```
We recommend starting with 'NEMO_test' which is the smallest in size and fastest to run.

To make it work, you need to download the input files for the tes cases, which are all provided on [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4709851.svg)](https://doi.org/10.5281/zenodo.4709851). You can get the files directly using wget as:
```bash
cd test_cases/oce
wget https://zenodo.org/record/4709851/files/NEMO_test.zip
unzip NEMO_test.zip && rm NEMO_test.zip
```

Other test cases can be downloaded as follows (note the file sizes):
```bash
# MITGCM_test (242M):
wget https://zenodo.org/record/4709851/files/MITGCM_test.zip
unzip MITGCM_test.zip && rm MITGCM_test.zip
# eORCA025_test (7.8G):
wget https://zenodo.org/record/4709851/files/eORCA025_test.zip
unzip eORCA025_test.zip && rm eORCA025_test.zip
# ROMS_test (257M):
wget https://zenodo.org/record/4709851/files/ROMS_test.zip
unzip ROMS_test.zip && rm ROMS_test.zip
```

Then, execute the script as follows:
```bash
python interpolate_to_common_grid_oce.py
```
This should create the following files:
* Oce3d\_NEMO3.6-IGE-CNRS-UGA\_a\_Ocean-A1\_201001-201002.nc
* OceSec\_NEMO3.6-IGE-CNRS-UGA\_a\_Ocean-A1\_201001-201002.nc
* OceMoor\_NEMO3.6-IGE-CNRS-UGA\_a\_Ocean-A1\_201001-201002.nc

On a laptop (16Gb, 2GHz), the test cases took the following durations:
* NEMO\_test : 48s (smaller than MISOMIP2 domain; all variables calculated). 
* MITGCM_test : 3min 30s (smaller than MISOMIP2 domain; no sea-ice or surface fluxes).
* eORCA025_test : 2min 30s (global simulation; most variables calculated).
* ROMS_test : 6min 20s (circum-Antarctic simulation (sigma coordinates); no sea-ice or surface fluxes).

### Adapt to your own ocean configuration

You can copy interpolate_to_common_grid_oce.py and adapt it to your model. You need to adapt at least section 0 (General information), section 1 (Files and variables) and section 2 (Global attributes of output netcdf). For section 1, you may need to create or modify a load\_oce\_mod\_xxxx.py function similar to the one existing for NEMO, MITgcm and ROMS if your model is not covered yet. If your model is already covered but variable names differ from what is understood by load\_oce\_mod\_xxxx.py, just add options for these variables in load\_oce\_mod\_xxxx.py and [make a pull request](https://opensource.com/article/19/7/create-pull-request-github) to upload it onto the official misomip2 repository (so that you keep this in case of updates).

### Performance of ocean interpolation

The scripts were written in best effort by non-experts in a way that should work for various grids, and suggestions that may boost perfromances are welcome. We nonetheless found that the computing times were acceptable for several example of simulations, the idea being that the higher the original model resolutoin, the shorter the period to be computed in a single call of interpolate\_to\_common\_grid\_oce.py. Here are some example of computing time:

* Processing the "NEMO-AMUXL025" data (231x190x75 points, reduced to 203x107x75 when "loading nemo") on a 2.6 GHz computer required ~10 minutes and ~2.5Gb of memory per year of simulation (tests up to 10 years in a single call of interpolate\_to\_common\_grid\_oce.py show that both walltime and memory increase linearly). 

* Processing the "NEMO-AMUXL12" data (687x567x75 points, reduced to 623x286x75 when "loading nemo") on a 2.6 GHz computer required 1h22 minutes and a bit less than 16Gb of memory per year of simulation.
 
### Ice test cases

**To be completed**

<br/><br/>
----------
----------
**Below are described the misomip2 functions used either for pre-processing model outputs or for analysing multpile models on common grids**

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

### misomip2.preproc.add\_standard\_attributes\_oce(ds,miss=9.969209968386869e36,verbose=False):
> Define standard netcdf attributes for ocean variables that are already present in the ds xarray dataset.
> (these variables must have the MISOMIP2 standard variable names)
> 
>    ds: xarray ocean dataset
>
>    miss: missing value (default=9.969209968386869e36)
>
>    verbose: True or False
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
>  Input:  
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
>    parallel: If True, the open and preprocess steps of this function will be performed in parallel [optional, default=False]
>
>  Output: 
>
>    xarray dataset of coordinates ("time", "z", "sxy") (sxy is the one-dimensionalized horizontal space)
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
> all variables required in MISOMIP2. It is adapted to NEMO's 
> C-grid (including stretched grids like eORCA) and standard 
> mesh/mask variables.
>
>  Input:
>
>    file\_mesh\_mask: mesh mask file [default='mesh\_mask.nc']
>
>    file\_bathy: bathymetry file [optional, look for data in file\_mesh\_mask if not provided]
>
>    files\_gridT: file or list of files containing the temperature and salinity variables
>
>    files\_gridU: file or list of files containing the x-velocity and related variables
>
>    files\_gridV: file or list of files containing the y-velocity and related variables
>
>    files\_ice: file or list of files containing the sea-ice variables [look for data in files\_gridT if not provided]
>
>    files\_SBC: file or list of files containing the surface fluxes variables [look for data in files\_gridT if not provided]
>
>    files\_BSF: file or list of files containing the barotropic streamfunction calculated at U-points, e.g. 
>        using the cdfpsi function which is part of the [CDFTOOLS](https://github.com/meom-group/CDFTOOLS).
>        [look for data in files_gridU if not provided]
>
>    rho0 corresponds to rau0 value in NEMO's eosbn2.F90 i.e. reference volumic mass [kg m-3]
>
>    teos10=False -> assumes the nemo outputs are in potential temperature & practical salinity (EOS80).<br/>
>          =True  -> assumes the nemo outputs are in CT and AS and convert them to PT and PS.
>
>    region = 'Amundsen' (default) or 'Weddell'.
>
>    parallel: If True, the open and preprocess steps of this function will be performed in parallel [optional, default=False]
>
>  Output:
>
>    xarray dataset of coordinates ("time", "z", "sxy") (sxy is the one-dimensionalized horizontal space)
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

### misomip2.preproc.horizontal\_interp( lat\_in\_1d, lon\_in\_1d, mlat\_misomip, mlon\_misomip, lat\_out\_1d, lon\_out\_1d, var\_in\_1d, weight=[], threshold=1.e20, skipna=False, filnocvx=False ):
> Interpolates one-dimension data horizontally to a 2d numpy array reshaped to the misomip standard (lon,lat) format.
>
>    Method: triangular linear barycentryc interpolation.
>
>    Input:
>
>       * lon\_in\_1d, lat\_in\_1d: 1d longitude and latitude of data to interpolate [xarray 1d data array]
> 
>       * mlat\_misomip, mlon\_misomip: misomip grid size (nb points) alond latitude and longitude dimensions
> 
>       * lon\_out\_1d, lat\_out\_1d: 1d longitude and latitude of the target misomip grid [numpy 1d data array]
>
>       * var\_in\_1d: 1d input data (same dimension as lon\_in\_1d and lat\_in\_1d) [xarray 1d data array]
>
>       * skipna = False to keep nans in interpolation, i.e. gives nan if any triangle node is nan [default]
>
>                = True to find interpolation triangle nodes with non-nan values
>
>       * filnocvx = True to use nearest-neighbor to fill non-convex areas, i.e. for which no triangulation is possible [default]
>
>                  = False to fill non-convex areas with nans 
>
>       * weight = weights used for interpolation [optional, xarray data array]
>
>       * threshold = threshold below which weight value indicates a masked point [default=1.e20]
>
>    Output:
>
>       * numpy data array of dimension (mlat_misomip, mlon_misomip)
>
_Examples_:
```bash
VAR_miso = mp.horizontal_interp( ds.latT, ds.lonT, mlat, mlon, lat_miso1d, lon_miso1d, ds.VAR )
VAR_miso = mp.horizontal_interp( ds.latT, ds.lonT, mlat, mlon, lat_miso1d, lon_miso1d, ds.VAR, skipna=True, weight=ds.siconc, threshold=0.1 )
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

