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
suggested to load the lattice informations from a file using use one of
these functions:

 - :py:func:`~muesr.io.cif.cif.load_cif` for Crystallographic Information Files
 - :py:func:`~muesr.io.cif.cif.load_mcif` for Magnetic Crystallographic Information Files
 - :py:func:`~muesr.io.xsf.xsf.load_xsf` for XCrysDen files.

.. note ::
   For CIF files the symmetry is automatically parsed and set from the file.
   For MCIF and XSF files it is not set and must be defined by hand or 
   with the :py:func:`~muesr.utilities.symsearch.symsearch` function 
   (only available is `spglib` is installed).


Defining a magnetic structure
-----------------------------

In :mod:`muesr` the magnetic structure is defined by the propagation 
vector `k`, the Fouerier components and the phases 
(see :ref:`intro_description_of_magnetic_structures`)

The propagation vector 

There are `many <http://magcryst.org/resources/magnetic-coordinates/>`_
conventions to specify the Fourier components of a magnetic structures:



At the current stage :mod:`muesr` supports the following coordinate 
systems:

0. Cartesian system
    the same used to define the lattice parameters which is implicitly 
    defined by the lattice vectors.

1. Bohr Magneton/Angstrom units, with x||a, y||b and z||c
    This is the reduced lattice coordinate system, where the magnetic 
    metric tensor (M) is the same metric used for interatomic distances 
    (G).

2. Bohr Magneton units, with x||a, y||b and z||c
    This is the crystal-axis coordinate system, where components of the
    moment are defined by their projections along the lattice basis
    vectors.
    If we define L = {{a,0,0},{0,b,0},{0,0,c}}, then the magnetic metric
    tensor is M = L.G.L^(-1), which is unitless.
    



Useful readings
+++++++++++++++
 - http://www.neutron-sciences.org/articles/sfn/pdf/2014/01/sfn201402001.pdf


Calculate local fields 
------------------------

The function simulating local fields at the muon site is 
:py:func:`~muesr.engines.clfc.locfield`. 

There are three type of simulations which are targeted to different
types of problems:

- `sum`: a simple sum of all the dipols in the Lorentz sphere.
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

 - CellError
 
To see their meaning follow the links.

The utility functions are mainly intented for interactive usage and report
problems 
