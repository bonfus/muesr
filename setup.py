#!/usr/bin/env python

from setuptools import setup, Extension

#from distutils.command.install import INSTALL_SCHEMES
#for scheme in INSTALL_SCHEMES.values():
#    scheme['data'] = scheme['purelib']

## In case of missing numpy array
try:
    from numpy import get_include as numpy_get_include
    numpy_include_dir = [numpy_get_include()]
except:
    numpy_include_dir = []
    
setup(name='muesr',
      version='0.1',
      description='Magnetic structure and mUon Embedding Site Refinement',
      author='Pietro Bonfa',
      author_email='pietro.bonfa@fis.unipr.it',
      url='https://www.no.url',
      packages=['muesr',
                'muesr.tests',
                'muesr.tests.core',
                'muesr.tests.utilities',
                'muesr.core',
                'muesr.io', 
                'muesr.io.xsf',
                'muesr.io.cif',
                'muesr.engines',
                'muesr.utilities',
                ],
      include_package_data=True,
      ext_modules=[Extension('lfcext', sources = ['muesr/engines/LFCExt/ass.c', \
                                      'muesr/engines/LFCExt/rass.c', \
                                      'muesr/engines/LFCExt/incass.c', \
                                      'muesr/engines/LFCExt/vec3.c', \
                                      'muesr/engines/LFCExt/mat3.c', \
                                      'muesr/engines/LFCExt/pile.c', \
                                      'muesr/engines/LFCExt/dt.c',\
                                      'muesr/engines/LFCExt/LFCExt.c'],
                                      libraries=['m'],
                                      include_dirs=numpy_include_dir,
                                      extra_compile_args=['-std=c99',],
                                      define_macros=[('_EXTENSION',None),])],
     package_dir={'muesr': 'muesr' },
     install_requires=[
          'numpy',
     ],
     test_suite="muesr.tests",
     )
