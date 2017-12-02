Magnetic structure and mUon Embedding Site Refinement
=====================================================

Muesr is a tool for quickly evaluating local fields at the muon sites in Muon Spin Rotation and Relaxation experiments (muSR).

Requirements
------------

Muesr works with the following minimal requirements:

| package | version    | url        |
|---------|------------|------------|
| python  | 2.7+ or 3  | python.org |
| numpy   | 1.6        | numpy.org  |
| lfclib  | 0.1        | [muLFC](http://www.github.com/bonfus/muLFC) |

Optional dependencies are:

| package  | version    | url        | provides |
|----------|------------|------------|----------|
| spglib   | 1.8        | [Spglib](http://atztogo.github.io/spglib) |  library for finding lattice structure symmetries |
| XCrysden | any        | [XCrysden](http://www.xcrysden.org) | tool for showing lattice and magnetic structures |
| sympy    | 1.0        | [Sympy](http://sympy.org) | for symbolic Fourier components definition |


Install and Usage
-----------------

For a detailed description see the [documentation](http://muesr.readthedocs.io/en/latest/Install.html).

The easiest way to install `lfclib` is using pre-baked wheels.
They are hosted in a test pip repository

    pip install --index-url https://testpypi.python.org/pypi mulfc
    
or at [packagecloud.io/muLFC/wheels](https://packagecloud.io/muLFC/wheels).

In the latter case choose the right package for your python installation. 
For example, `LFC-0.1-cp36-cp36m-manylinux1_i686.whl` is compatible with
Python 3.6 (cp36) and a 32 bit system, while `LFC-0.1-cp34-cp34m-manylinux1_x86_64.whl`
is for Python 3.4 installed on a 64 bit system.
See [documentation](http://muesr.readthedocs.io/en/latest/Install.html#system-wide-installation-with-wheels) for further details.

Known problems
--------------

- Non-standard spacegroup settings can cause some tedious problems when 
  using the spglib functions. Pay attention!
- OpenMP implementation is EXPERIMENTAL (passes tests but needs more work)


Please note that the code is still under heavy development. 
You'll probably find bugs so please report them!

Authors & Contributors
----------------------

Pietro Bonfa', Ifeanyi John Onuorah, Anthony Lim and Roberto De Renzi.


Notes
-----

Part of the code in this repository is from the ASE
(https://wiki.fysik.dtu.dk/ase/index.html) project. 
