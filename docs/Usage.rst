Usage
=====

Defining a sample
-----------------

This is very easy! You just do: ::


    >>> from muesr.core.sample import Sample
    >>> smp = Sample()
    >>> smp.name = "My very interesting experiment on ..."

the :py:attr:`~muesr.core.sample.Sample.name` property is optional.


Defining a lattice structure
----------------------------

This can be done at code level by using the
:py:attr:`~muesr.core.sample.Sample.cell` property.
For example ::

    >>> from muesr.core.sample import Sample
    >>> from muesr.core.atoms import Atoms
    >>> smp = Sample()
    >>> smp.cell = Atoms(...)

where the dots replace the :py:class:`~muesr.core.atoms.Atoms` 
class initialization arguments.

However this can be tedious and error prone. Therefore you are strongly
suggested to load the lattice information from a file using use one of
these functions:

 - :py:func:`~muesr.i_o.cif.cif.load_cif` for Crystallographic Information Files
 - :py:func:`~muesr.i_o.cif.cif.load_mcif` for Magnetic Crystallographic Information Files
 - :py:func:`~muesr.i_o.xsf.xsf.load_xsf` for XCrysDen files.

.. note ::
   For CIF files the symmetry is automatically parsed and set from the file.
   For MCIF and XSF files it is not set and must be defined by hand or 
   with the :py:func:`~muesr.utilities.symsearch.symsearch` function 
   (only available if `spglib` is installed).


Defining a magnetic structure
-----------------------------

In :mod:`muesr` the magnetic structure is defined by the propagation 
vector `k`, the Fourier components and the phases 
(see :ref:`intro_description_of_magnetic_structures`)

Propagation vector
++++++++++++++++++

In :mod:`muesr`, the propagation vector is always specified 
in **reciprocal lattice units** with the 
:py:attr:`~muesr.core.magmodel.MM.k` property.

Fourier components
++++++++++++++++++

There are `many <http://magcryst.org/resources/magnetic-coordinates/>`_
conventions to specify the Fourier components of a magnetic structures.
At the current stage :mod:`muesr` supports the following coordinate 
systems:

0. Cartesian system
    the same used to define the lattice parameters which is implicitly 
    defined by the lattice vectors.

1. Bohr Magneton/Angstrom units, with x||a, y||b and z||c
    This is the reduced lattice coordinate system, where the magnetic 
    metric tensor (M) is the same metric used for inter-atomic distances 
    (G).

2. Bohr Magneton units, with x||a, y||b and z||c
    This is the crystal-axis coordinate system, where components of the
    moment are defined by their projections along the lattice basis
    vectors.
    If we define L = {{a,0,0},{0,b,0},{0,0,c}}, then the magnetic metric
    tensor is M = L.G.L^(-1), which is unit-less.

Here's a table connecting the three possible input and related functions

.. table ::

   ======================== ============================================ ===========================================================
   coordinate system number :py:mod:`~muesr.core.magmodel.MM` property   String version (used in YAML files and in helper functions)
   ======================== ============================================ ===========================================================
   0                        :py:attr:`~muesr.core.magmodel.MM.fc`        'bohr-cartesian' or 'b-c' (case insensitive)
                            :py:attr:`~muesr.core.magmodel.MM.fcCart`    
   1                        :py:attr:`~muesr.core.magmodel.MM.fcLattBMA` 'bohr/angstrom-lattice' or 'b/a-l' (case insensitive)
   2                        :py:attr:`~muesr.core.magmodel.MM.fcLattBM`  'bohr-lattice' or 'b-l' (case insensitive)
   ======================== ============================================ ===========================================================

Quick overview
++++++++++++++

To define a new magnetic structure just do ::

    >>> from muesr.core.sample import Sample
    >>> from muesr.core.magmodel import MM
    >>> smp = Sample()
    >>> 
    >>> # load a lattice structure!
    >>>
    >>> smp.new_mm()

The newly created magnetic structure is automatically selected as the 
current magnetic model and can be obtained with the 
:py:attr:`~muesr.core.sample.Sample.mm` property.
From that you can access and define all the properties of the magnetic definition ::

    >>> smp.mm.k
    ... array([0, 0, 0])
    

The three fundamental properties of a magnetic model are:

  - :py:attr:`~muesr.core.magmodel.MM.fc`
  - :py:attr:`~muesr.core.magmodel.MM.k`
  - :py:attr:`~muesr.core.magmodel.MM.phi`

Please see the :py:mod:`muesr.core.magmodel.MM` documentation for the 
details.

To simplify the definition of the magnetic structure, the 
:py:func:`~muesr.utilities.ms.mago_add` helper function is available in
the :py:mod:`muesr.utilities.ms` module.

It prompts an interactive interface like the one shown below
(for a Ti2O3 structure): ::

    >>> from muesr.utilities.ms import mago_add
    >>> mago_add(smp,coordinates='bohr-lattice')
    ...      Propagation vector (w.r.t. conv. rec. cell): 0 0 0
    ... 	 Magnetic moments in bohr magnetons and lattice coordinates.
    ... 	 Which atom? (enter for all)Ti
    ... 	 Lattice vectors:
    ... 	   a    5.149000000000000    0.000000000000000    0.000000000000000
    ... 	   b   -2.574499999999999    4.459164804086075    0.000000000000000
    ... 	   c    0.000000000000001    0.000000000000001   13.641999999999999
    ... 	 Atomic positions (fractional):
    ... 	     1 Ti  0.00000000000000  0.00000000000000  0.34500000000000  47.867
    ... 	     2 Ti  0.66666666666667  0.33333333333333  0.67833333333333  47.867
    ... 	     3 Ti  0.33333333333333  0.66666666666667  0.01166666666667  47.867
    ... 	     4 Ti  0.00000000000000  0.00000000000000  0.84500000000000  47.867
    ... 	     5 Ti  0.66666666666667  0.33333333333333  0.17833333333333  47.867
    ... 	     6 Ti  0.33333333333333  0.66666666666667  0.51166666666667  47.867
    ... 	     7 Ti  0.00000000000000  0.00000000000000  0.15500000000000  47.867
    ... 	     8 Ti  0.66666666666667  0.33333333333333  0.48833333333333  47.867
    ... 	     9 Ti  0.33333333333333  0.66666666666667  0.82166666666667  47.867
    ... 	    10 Ti  0.00000000000000  0.00000000000000  0.65500000000000  47.867
    ... 	    11 Ti  0.66666666666667  0.33333333333333  0.98833333333333  47.867
    ... 	    12 Ti  0.33333333333333  0.66666666666667  0.32166666666667  47.867
    ... 	 FC for atom 1 Ti (3 real, [3 imag]): 0 1 0
    ... 	 FC for atom 2 Ti (3 real, [3 imag]): 0 1 0
    ... 	 FC for atom 3 Ti (3 real, [3 imag]): 0 1 0
    ... 	 FC for atom 4 Ti (3 real, [3 imag]): 0 1 0
    ... 	 FC for atom 5 Ti (3 real, [3 imag]): 0 1 0
    ... 	 FC for atom 6 Ti (3 real, [3 imag]): 0 1 0
    ... 	 FC for atom 7 Ti (3 real, [3 imag]): 0 1 0
    ... 	 FC for atom 8 Ti (3 real, [3 imag]): 0 1 0
    ... 	 FC for atom 9 Ti (3 real, [3 imag]): 0 1 0
    ... 	 FC for atom 10 Ti (3 real, [3 imag]): 0 1 0
    ... 	 FC for atom 11 Ti (3 real, [3 imag]): 0 1 0
    ... 	 FC for atom 12 Ti (3 real, [3 imag]): 0 1 0
    ... 
    
This produces the following Fourier components in Cartesian coordinates ::

    >>> smp.mm.fc
    ... array([[-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [-0.5000000+0.j,  0.8660254+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j],
    ...        [ 0.0000000+0.j,  0.0000000+0.j,  0.0000000+0.j]])

Which are indeed: ::

    >>> smp.mm.fcLattBM
    ... array([[ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  1.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j],
    ...        [ 0.+0.j,  0.+0.j,  0.+0.j]])


The zeros in the Fourier components are from the atoms different from 
`Ti`.

  .. note::
     The phases can only be set with the :py:attr:`~muesr.core.magmodel.MM.phi`
     property.
     

Useful readings
+++++++++++++++
 - http://www.neutron-sciences.org/articles/sfn/pdf/2014/01/sfn201402001.pdf


Setting the muon position
-------------------------

The muon position can be easily set with the 
:py:attr:`~muesr.core.sample.Sample.add_muon` method.

If symmetry is defined, equivalent muon positions can be obtained with 
the function :py:func:`~muesr.utilities.muon.muon_find_equiv` in the 
:py:mod:`muesr.utilities.muon` module.

Calculate local fields 
------------------------

The function simulating local fields at the muon site is 
:py:func:`~muesr.engines.clfc.locfield`. 

There are three type of simulations which are targeted to different
types of problems:

- `sum`: a simple sum of all the magnetic moments in the Lorentz sphere.
- `rotate`: rotates the local moments around a given axis and perform the
  sum. This function offer great flexibility in the way local moments
  are rotated but is not computationally efficient. For incommensurate
  magnetic orders the following function is much more efficient.
- `incommensurate`: Fast version of 'rotate' which exploits the method 
  discussed in Phys. Rev. B **93**, 174405 (2016).



Calculate the dipolar tensor
----------------------------

The function providing the dipolar tensor at the muon site is 
:py:func:`~muesr.engines.clfc.dipten`.

To use it you have to specify a (arbitrary) value for the Fourier components
of the magnetic atoms that you want to include in the sum. The specified 
value has no meaning only the 0 vs different from zero has.

.. note::
   Results are provided in Angstrom^-3 !


Generate grid of interstitial points for DFT simulations
---------------------------------------------------------

A useful function to prepare the input for DFT simulations is 
:py:func:`~muesr.utilities.dft_grid.build`.

The function provides a set of symmetry inequivalent interstitial 
positions with the additional constraint of being sufficiently separated
from the atoms of the hosting system.


Understanding errors
--------------------

:mod:`muesr` raises the conventional python exceptions (mainly ValueError and
TypeError) or other 4 specific Exceptions:

 - :py:class:`~muesr.core.sampleErrors.CellError`
 - :py:class:`~muesr.core.sampleErrors.MuonError`
 - :py:class:`~muesr.core.sampleErrors.MagDefError`
 - :py:class:`~muesr.core.sampleErrors.SymmetryError`
 
To see their meaning follow the links.

N.B.: the utility functions are mainly intended for interactive usage 
and therefore report problems by printing error messages on the screen.
Exceptions are only raised in core components.


Saving and loading sample details to/from file
----------------------------------------------

To save a sample use :py:func:`~muesr.i_o.sampleIO.save_sample`. To load
a saved sample use :py:func:`~muesr.i_o.sampleIO.load_sample`.

Data is stored in an YAML file. It is possible (but error prone) to write
an input file by hand. When loaded, the file will undergo a minimal 
validation. Identifying the errors is not so easy so the best method to specify
the sample details is probably using the various functions discussed in 
this manual.
