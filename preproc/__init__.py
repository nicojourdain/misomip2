
# to enable direct call to functions:
from .def_grids import grid_bounds
from .def_grids import generate_3d_grid
from .def_grids import generate_section_grid
from .def_grids import generate_mooring_grid

from .def_attrs import add_standard_attributes

from .load_ocean_data_nemo import load_oce_mod_nemo
from .load_ocean_data_mitgcm_stereo import load_oce_mod_mitgcm_stereo

__version__ = '0.1'
