LiFePO4
=======

In this example we will go through the relevant lines of the script
`run_example.py` in the `LiFePO4` directory of the examples.
This example shows how to calculate the local fields at the muon sites in
LiFePO4 and loosely follows the description of Ref. [Sugiyama2011]_.

Interactive magnetic order definition
-------------------------------------

In this section we will prepare a script that asks the user to specify
the magnetic order at runtime and prints the output on screen.

Let's first import the necessary python objects and functions of :py:mod:`muesr`.

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 14-17
   :lineno-start: 14
   :language: python


To create a new sample definition and initialize it, just do:

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 28
   :lineno-start: 28
   :language: python

The sample object holds the following information:
  * Lattice structure
  * Magnetic Orders
  * Muon positions
  * Symmetry (optional)

The lattice structure is orthorhombic, with *Pnma* symmetry and
lattice parameters :math:`a=10.3244(2), b=6.0064(3), c=4.6901(5)`.
We can import these data from a CIF file using the function
:py:func:`~muesr.i_o.cif.cif.load_cif`:

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 31
   :lineno-start: 31
   :language: python
   
   
The next step is setting the the definition of the magnetic order. This 
can be done interactively during the script execution or programmatically.
Let's discuss the former case first. The function :py:func:`~muesr.utilities.ms.mago_add`
will prompt the input and let you specify the Fourier components and the
propagation vector.
Here we want to describe an anti-ferromagnetic order in which the
iron moments lie along the :math:`y` axis and are
4.19 :math:`\mu_{\mathrm{B}}` in magnitude.

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 43-61
   :lineno-start: 43
   :emphasize-lines: 1
   :language: python

In order to complete the setup we specify the muon positions. 
The method :py:attr:`~muesr.core.sample.Sample.add_muon` inputs lattice
(fractional) coordinates by default. Four muon sites are discussed by 
Jun Sugiyama `et al.` [Sugiyama2011]_.

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 67-70
   :lineno-start: 67
   :language: python
   

.. note:: The muon sites will be automatically placed in the central 
          unit cell of the supercell that will be used for the subsequent
          calculations.

The local fields are finally calculated with the :py:func:`~muesr.engines.clfc.locfield`
command.

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 80-81
   :emphasize-lines: 2
   :lineno-start: 80
   :language: python

.. warning:: * Always check convergence against the supercell size.
             * Always use a Lorentz sphere that can be inscribed in the 
               selected supercell.
   

Please see the documentation of the function :py:func:`~muesr.engines.clfc.locfield`
to see a description of the input parameters. As it is immediately evident,
the supercell size and the Lorentz radius are not the ideal choices.
Improving this parameters is left as an exercise to the reader.

The results are stored in a list of 
:py:class:`~muesr.core.magmodel.LocalField` objects which, for each muon
site, contain the Total, Dipolar, Lorentz and Contact contributions in
**Tesla**.

The next four lines of code print the results of the simulation in T/:math:`\mu_{\mathrm{B}}` 

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 84-87
   :lineno-start: 84
   :language: python
   
   
.. note::   **The results are always** reported in the **Cartesian** 
            coordinate system defined by the lattice vectors of the
            crystal.


Programmatic magnetic order definition
---------------------------------------

As already mentioned above, it's also possible to scepcify a magnetic 
order programmatically. This can be done with the help of the methods
:py:attr:`~muesr.core.sample.Sample.new_mm`, 
:py:attr:`~muesr.core.magmodel.MM.k` and :py:attr:`~muesr.core.magmodel.MM.fc`.

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 95-141
   :emphasize-lines: 2,9,13,47
   :lineno-start: 95
   :language: python

Please remember to specify the FC in a 2D array with 3 columns and :math:`N_{\mathrm{Atoms}}`
rows of complex values.

.. note::  The **propagation vector** is **always** specified in **recirocal lattice units**. 
           On the other hand, the **Fourier components** can be specified with **three
           different coordinate system and units**: 
           
              1. Bohr Magnetons in Cartesian coordinates (Cartesian vector notation)
              2. Bohr Magnetons/Angstrom along the lattice vectors (Lattice vector notation)
              3. Modulus (in Bohr Magnetons) along the lattice vectors.

The results can be retrived as described above.

.. [Sugiyama2011] Jun Sugiyama `et al.`, Phys. Rev. B 84, 054430 (2011)
