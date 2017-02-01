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


Optional dependencies are:

| package  | version    | url        | provides |
|----------|------------|------------|----------|
| spglib   | 1.8        | [Spglib](http://atztogo.github.io/spglib) |  library for finding lattice structure symmetries |
| XCrysden | any        | [XCrysden](http://www.xcrysden.org) | tool for showing lattice and magnetic structures |
| sympy    | 1.0        | [Sympy](http://sympy.org) | for symbolic Fourier components definition |


Install and Usage
-----------------

See the documentation at http://muesr.readthedocs.io

Known problems
--------------

- Non-standard spacegroup settings can cause some tedious problems when 
  using the spglib functions. Pay attention!
- OpenMP implementation is EXPERIMENTAL (passes tests but needs more work)


Please not that the code is still under heavy development. 
You'll probably find bugs so please report them!

Authors
-------

Pietro Bonfa', Ifeanyi John Onuorah and Roberto De Renzi.

Notes
-----

Part of the code in this repository is from the ASE
(https://wiki.fysik.dtu.dk/ase/index.html) project. 
