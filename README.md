# misomip2
A python package to postprocess model outputs to standard MISOMIP2 format and to analyse MISOMIP2 multi-model outputs.

## Contributors
Nicolas C. Jourdain (IGE, CNRS-UGA, Grenoble, France)


## Preprocessing
Contains scripts that facilitate interpolation and formatting to the MISOMIP2 standards.

To use the preprocessing tools, start by specifying:
```bash
import misomip2.preproc as misopre
```

**generate\_3d\_grid**(region='Amundsen'):
> Generates (longitude, latitude, depth) of the common MISOMIP2 3d grid
> 
>    region: 'Amundsen' (default), 'Weddell'
>
>    exemple: [lon,lat,depth]=generate_3d_grid(region='Amundsen')

**generate\_section\_grid**(region='Amundsen'):
> Generates (longitude, latitude, depth) of the common MISOMIP2 section
> 
>     region: 'Amundsen' (default), 'Weddell'
> 
>     exemple: [lon,lat,depth]=generate_section_grid(region='Amundsen')

**generate_mooring_grid**(region='Amundsen'):
> Generates (longitude, latitude, depth) of the common MISOMIP2 mooring
>
>    region: 'Amundsen' (default), 'Weddell'
> 
>    exemple: [lon,lat,depth]=generate_mooring_grid(region='Amundsen')


## Multi-model Analysis 
Contains scripts to quickly plot multi-model diagnostics.

To use the analysis tools, start by specifying:
```bash
import misomip2.analysis as misoan
```

## Examples
