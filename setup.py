#!/usr/bin/env python

import os
import re
from setuptools import setup, find_packages


here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'misomip2', '__init__.py')) as f:
    init_file = f.read()
version = re.search(r'{}\s*=\s*[(]([^)]*)[)]'.format('__version_info__'),
                    init_file).group(1).replace(', ', '.')

setup(name='misomip2',
      version=version,
      description='A python package to postprocess model outputs '
                  'to standard MISOMIP2 format and to analyse '
                  'MISOMIP2 multi-model outputs.',
      url='https://github.com/nicojourdain/misomip2',
      author='Nicolas C. Jourdain',
      author_email='nicolas.jourdain@univ-grenoble-alpes.fr',
      license='MIT',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: MIT License',
          'Operating System :: OS Independent',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering',
      ],
      packages=find_packages(),
      package_data={'misomip2': ['config.default']},
      install_requires=['numpy', 'scipy', 'matplotlib', 'netCDF4', 'xarray',
                        'dask', 'pandas', 'cmocean'],
      entry_points={'console_scripts':
                    ['misomip2 = misomip2.__main__:main']})
