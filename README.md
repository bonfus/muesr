Magnetic structure and mUon Embedding Site Refinement
=====================================================

Muesr is a tool for quickly evaluating local fields at the muon sites in Muon Spin Rotation and Relaxation experiments (muSR).

Requirements
------------

Muesr works with the following minimal requirements:

| package | version    | url        |
|---------|------------|------------|
| python  | 2.7+ or 3  | python.org |
| numpy   | 1.6+       | numpy.org  |
| lfclib  | 0.1        | [muLFC](http://www.github.com/bonfus/muLFC) |

Optional dependencies are:

| package  | version    | url        | provides |
|----------|------------|------------|----------|
| spglib   | 1.8+       | [Spglib](http://atztogo.github.io/spglib) |  library for finding lattice structure symmetries |
| XCrysden | any        | [XCrysden](http://www.xcrysden.org) | tool for showing lattice and magnetic structures |
| sympy    | 1.0+       | [Sympy](http://sympy.org) | for symbolic Fourier components definition |

Install and Usage
-----------------

The easiest way to install `muesr` is using pip:

    pip install -r requirements.txt muesr

For a detailed description see the [documentation](http://muesr.readthedocs.io/en/latest/Install.html).

Known problems
--------------

- Non-standard spacegroup settings can cause some tedious problems when 
  using the spglib functions. Pay attention!

Please note that the code is still under development. 
You'll probably find bugs, please report them.

Authors & Contributors
----------------------

Pietro Bonfa', Ifeanyi John Onuorah, Anthony Lim and Roberto De Renzi.

Notes and License
-----------------

Part of the code in this repository comes from the [ASE](https://wiki.fysik.dtu.dk/ase/index.html) project.
Where not differently specifies, the source code is provided under the GPLv3 license.
