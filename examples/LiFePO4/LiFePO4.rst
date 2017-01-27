LiFePO4
-------

In this example we will go through the relevant lines of the script
`run_example.py` in the `LiFePO4` directory of the examples.
This example shows how to calculate the local fields at the muon site in
ferromagnetic LiFePO4.

To prepare the environment, first load the necessary python objects and
functions of `mod`:muesr.

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 14-17
   :lineno-start: 14
   :language: python


The sample object holds the following information:
  * Lattice structure
  * Magnetic Orders
  * Muon positions
  * Symmetry (optional)
  
To create a new sample definition and initialize it, just do:

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 28
   :lineno-start: 28
   :language: python

the lattice structure can be imported from a CIF file using the function
:py:func:`~muesr.io.cif.cif.load_cif`:

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 31
   :lineno-start: 31
   :language: python
   
   
The next step is setting the the definition of the magnetic order. This 
can be done interactively during the script execution or programmatically.
Let's discuss the former case first. The function :py:func:`~muesr.utilities.ms.mago_add`
will prompt the input and let you specify the Fourier components and the
propagation vector.
In this case the iron moments lie along the y axis and are 4.19 Bohr
magnetons in magnitude.

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
   

The local fields are finally calculated with the :py:func:`~muesr.engines.clfc.locfield`
command.

.. literalinclude:: ../examples/LiFePO4/run_example.py
   :lines: 80-81
   :emphasize-lines: 2
   :lineno-start: 80
   :language: python
   

The results are stored in a list of 
:py:class:`~muesr.core.magmodel.LocalField` objects which, for each muon
site, contain the Total, Dipolar, Lorentz and Contact contributions in
*Tesla*.



.. [Sugiyama2011] Jun Sugiyama `et al.`, Phys. Rev. B 84, 054430 (2011)
