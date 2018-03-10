#!/usr/bin/env python
from setuptools import setup
long_description ="""
Muesr is a python package that performs
dipole-field and dipole-tensor sums using real space algorithms for
Muon Spin Rotation and Relaxation Spectroscopy.
"""
setup(name='muesr',
      version='0.1.2',
      description='Magnetic structure and mUon Embedding Site Refinement',
      long_description=long_description,
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
      package_dir={'muesr': 'muesr' },
      install_requires=[
            'numpy >= 1.6',
      ],
      test_suite="muesr.tests",
     )
