Tutorial
========

With the help of a few examples, we will show how to quickly evaluate local fields in Muesr.

To use :py:mod:`~muesr` you must be familiar with python. An interactive shell like ipython or jupyter can 
help a lot but it is not needed.

First steps with muesr
---------------------------

Definig the sample
+++++++++++++++++++++++++++++++++

The fundamental component of muesr is the :py:class:`muesr.core.sample.Sample` object.
You can import and instantiate it like this:

.. code-block:: python
    
    >>> from muesr.core import Sample
    >>>
    >>> mysample = Sample()

Specifying the atomic structure
++++++++++++++++++++++++++++++++++++

The first thing that must be defined is the atomic structure. Muesr uses 
the ASE :class:`~Atoms` class. You can declare your own, fore example
like this (Copper in simple cubic lattice)

.. code-block:: python
    
    >>> import numpy
    >>> from muesr.core.atoms import Atoms
    >>> 
    >>> atms = Atoms(symbols=['Cu'], scaled_positions=[[0.,0.,0.]], cell=numpy.diag([3.,3.,3.]), pbc=True)
    >>> 
    >>> mysample.cell = atms
    
However this is quite tedious and error prone so it is much better to use some
builtin functions to parse crystallographic files.

At the moment musr can parse XCrysDen files (xsf), CIF (cif) and mCIF (mcif)
files. Here's a few examples:

.. code-block:: python
    
    >>> # load data from XCrysden .xsf file
    >>> from muesr.i_o.xsf.xsf import load_xsf
    >>> 
    >>> load_xsf(mysample, "/path/to/file.xsf")
    >>> 
    >>> 
    >>> # load data from .cif file
    >>> from muesr.i_o.cif.cif import load_cif
    >>> 
    >>> load_cif(mysample, "/path/to/file.cif")
    >>> 
    >>> 
    >>> # load data from .mcif file
    >>> from muesr.i_o.cif.cif import load_mcif
    >>> 
    >>> load_mcif(mysample, "/path/to/file.mcif")


The :py:func:`~muesr.i_o.cif.cif.load_cif` function will also load symmetry information. 
Please note that only a single lattice structure at a time can be
defined so each load function will remove the previous lattice structure
definition.

Setting muon positions
++++++++++++++++++++++

When the lattice structure is defined it is possible to specify the
muon position and the magnetic orders.

To specify the muon position, simply do:

.. code-block:: python
    
    >>> mysample.add_muon([0.1,0,0])
    
positions are assumed to be in fractional coordinates. Cartesian coordinates
can be specified as

.. code-block:: python
    
    >>> mysample.add_muon([0.1,0,0], cartesian=True)

If proper symmetry of the sample is present in the sample definition, it
is usually usefull to get symmetry equivalent sites.
This can be done with the utility function :py:func:`~muesr.utilities.muon.muon_find_equiv`.

.. code-block:: python
    
    >>> from muesr.utilities import muon_find_equiv
    >>> muon_find_equiv(mysample)


Defining a magnetic structure
++++++++++++++++++++++++++++++

The next step is the definition of a magnetic structure. To do so one 
must specify the propagation vector and the Fourier components and, 
optionally, the phases.
A quick way to do that is using the helper function :py:func:`~muesr.utilities.ms.mago_add` from
:py:mod:`~muesr.utilities.ms`. 

.. code-block:: python
    
    >>> from muesr.utilities.ms import mago_add
    >>> 
    >>> mago_add(mysample)
    
You will be asked the propagation vector and the Fourier coefficients
for the specified atomic symbol. By default the Fourier components are
specified in **Cartesian** coordinates. You can use the keyword argument
`inputConvention` to change this behaviour.
Here's an example::

     >>> mago_add(a)
        Propagation vector (w.r.t. conv. rec. cell): 0 0 0
        Magnetic moments in bohr magnetons and cartesian coordinates.
        Which atom? (enter for all)Cu
        Lattice vectors:
            a    5.000000000000000    0.000000000000000    0.000000000000000
            b    0.000000000000000    5.000000000000000    0.000000000000000
            c    0.000000000000000    0.000000000000000    5.000000000000000
        Atomic positions (fractional):
            1 Cu  0.00000000000000  0.00000000000000  0.00000000000000  63.546
        FC for atom 1 Cu (3 real, [3 imag]): 0 0 1
        
The same can be achieved in a more pythonic way like this:

.. code-block:: python
    
    >>> mysample.new_mm()
    >>> mysample.k = numpy.array([ 0.,  0.,  0.])
    >>> mysample.fc = numpy.array([[ 0.+0.j,  0.+0.j,  1.+0.j]])

.. note::
   In this method each atom must have a Fourier component! For a 8 atoms
   unit cell the numpy array specifying the value must be a 8 x 3 complex
   array!
   

It is possible to specify multiple magnetic structure for the same lattice
structure. Each time a new magnetic structure is added or set the 
previously specified magnetic orders are kept.


Checking the magnetic structure
+++++++++++++++++++++++++++++++

The Fourier components are complex vector and therefore not so easy to 
visualize. There are two ways to actually see the magnetic moment 
defined in the system. One is to generate a (possibly trivial) supercell 
and visualize it in XCrysDen. The other is to use FPStudio.


Evaluating the local field
++++++++++++++++++++++++++

Once you are done with the definition of the sample details it's time to
crunch some numbers!
To evaluate the local fields at the muon site :py:mod:`~muesr` uses a 
python extension written in C in order to get decent performances.
You can load a simple wrapper to the extension as providing local fields
with the following command ::

    >>> from muesr.engines.clfc import locfield

A detailed description of the possible calculatros is given in the 
:py:func:`~muesr.engines.clfc.locfield` documentation.

Let's go strainght to the local field evaluation which is obtained by 
running the command: ::

    >>> results = locfield(mysample, 'sum', [30, 30, 30] , 100)

Let's review this command in details. The first argument is just the 
sample object that we considered till now.
The second argument tells the code to simpli sum all magnetic moments
in a supercell generated by the expansion of the unitc cell 30x30x30 
times along the lattice vectors (third argument of the function).
The fourth argument is the radius of the Lorentz sphere considered.
All magnetic moments outside the Lorentx sphere are ignored.
The muon is automatically placed in the center of the supercell. 

.. note::
   To get an estimate of the largest radius that you can use to avoid 
   sampling outside the supercell size you can use :py:func:`~muesr.engines.clfc.find_largest_sphere`

The results variable now contains a list of 
:py:class:`~muesr.core.magmodel.LocalField` objects.
However if you print the results you'll see something which looks like
a numpy array: ::

    >>> print(results)
    [array([  3.00964434e-06,  -8.91586975e-20,  -7.39755731e+00])]
    
the numbers shown here are the total field for the magnetic structure 
discussed above. To access the various componenets you do: ::

    >>> results[0].Lorentz
    array([ 0.        ,  0.        ,  0.02303293])
    
    >>> results[0].Dipolar
    array([  3.00964434e-06,  -8.91586975e-20,  -7.42059024e+00])
    
    >>> results[0].Contact
    array([ 0.,  0.,  0.])


And you are done! Remember that all results are in Tesla units.

In the next tutorial we will discuss the Hyperfine Contact Field.

The Contact field contribution
------------------------------

TODO






