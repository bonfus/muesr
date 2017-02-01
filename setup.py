#!/usr/bin/env python
import sys
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

openmp_compile_args = []
openmp_link_args = []


## In case of missing numpy array headers
try:
    from numpy import get_include as numpy_get_include
    numpy_include_dir = [numpy_get_include()]
except:
    numpy_include_dir = []

# Ugly hack to set compiler flags 
COMPILE_ARGS = {'msvc':[],'gcc':[],'unix':[]}
LINK_ARGS = {'msvc':[],'gcc':[],'unix':[]}

for compiler, args in [
        ('msvc', ['/EHsc', '/DHUNSPELL_STATIC']),
        ('gcc', ['-O3', '-g0', '-std=c99']),
        ('unix', ['-O3', '-g0', '-std=c99'])]:
    COMPILE_ARGS[compiler] += args
    

# Ugly hack to have openMP as option
if "--with-openmp" in sys.argv:
    for compiler, args in [
            ('msvc', ['/openmp']),
            ('unix', ['-fopenmp']),
            ('gcc', ['-fopenmp'])]:
        COMPILE_ARGS[compiler] += args    
    for compiler, args in [
            ('msvc', []),
            ('unix', ['-lgomp']),
            ('gcc', ['-lgomp'])]:
        LINK_ARGS[compiler] += args    

    sys.argv.remove("--with-openmp")

class build_ext_compiler_check(build_ext):
    def build_extensions(self):
        
        compiler = self.compiler.compiler_type
        cargs = COMPILE_ARGS[compiler]
        for ext in self.extensions:
            ext.extra_compile_args = cargs
            
        largs = LINK_ARGS[compiler]
        for ext in self.extensions:
            ext.extra_link_args = largs
        
        build_ext.build_extensions(self)


    
setup(name='muesr',
      version='0.1',
      description='Magnetic structure and mUon Embedding Site Refinement',
      author='Pietro Bonfa',
      author_email='pietro.bonfa@fis.unipr.it',
      url='https://github.com/bonfus/muesr',
      packages=['muesr',
                'muesr.tests',
                'muesr.tests.core',
                'muesr.tests.utilities',
                'muesr.core',
                'muesr.i_o', 
                'muesr.i_o.xsf',
                'muesr.i_o.cif',
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
                                      define_macros=[('_EXTENSION',None),])],
     package_dir={'muesr': 'muesr' },
     install_requires=[
          'numpy',
     ],
     test_suite="muesr.tests",
     cmdclass={ 'build_ext': build_ext_compiler_check }
     )
