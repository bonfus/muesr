MnSi
====

This example will guide you through the experimental investigation conducted by Amato et al. and Dalmas de Réotier et. al reported in the following journal articles [Amato2014]_, [DalmasdeReotier2016]_.

.. note::  This is an advanced example. It is assumed that you already
           familiarized with Muesr by following one of the other examples
           or the tutorial.


Scientific background
----------------------

MnSi has a cubic lattice structure with lattice constant 4.558 Å. 
It has P213 (No. 198) space group symmetry and the Mn-ion occupy the 
position (0.138,0.138,0.138), while the Si-ion the position (0.845,0.845,0.845).


At T ~ 0 K, the magnetic structure of MnSi is characterized by spins forming a 
left-handed incommensurate helix with a propagation vector k≃0.036 Å^−1
in the [111] direction [5–7]. The static Mn moments ( ∼0.4μB for T→0 K)
point in a plane perpendicular to the propagation vector. 

Dipolar Tensor
--------------

Given the number of oscillations that are found in the experiment, only a Wyckoff site of type 4a is possible.
We therefore proceed by inspecting the dipolar tensor for all points the points of type 4a which happen
to be in the 111 (and equivalent) direction of the cubic cell.
As a first step, we will identify the number of parameters that characterize the dipolar tensor numerically by 
visually checking the dipolar tensor elements for all the equivalent sites. This can be done much more accurately
analytically but a numerical check can come in handy.

.. literalinclude:: ../examples/MnSi/run_example.py
   :lines: 1-78
   :emphasize-lines: 39,42,46-53,56-59,62
   :lineno-start: 1
   :language: python
   
   
In order to compare with the experiment, we will evaluate the dipolar 
tensor for 100 values along the 111 direction. 

.. literalinclude:: ../examples/MnSi/run_example.py
   :lines: 83-97
   :emphasize-lines: 3,6
   :lineno-start: 83
   :language: python


Lines 99-110 produce the following figure: 
|dipten|

.. |dipten| image:: ../examples/MnSi/reference/DipolarTensor.png


Local Fields
-------------
In order to calculate the local fields at the muon site in the helical state, we will first define the magnetic order and then evaluate the local fields with a optimized algorithm for incommensurate magnetic structures.

Let us first create a few useful variables for the definition of the Fourier components (FC).

.. literalinclude:: ../examples/MnSi/run_example.py
   :lines: 115-146
   :lineno-start: 115
   :language: python

We can now define both a left handed and a right handed helix.

.. literalinclude:: ../examples/MnSi/run_example.py
   :lines: 150-169
   :lineno-start: 150
   :language: python

We finally add the muon positions and the two magnetic orders with the 
commands

.. literalinclude:: ../examples/MnSi/run_example.py
   :lines: 175-188
   :lineno-start: 175
   :language: python

.. note:: When a new magnetic order is added, the index is automatically
          incremented and the new entry is immediately selected



In order to get the local contributions to the magnetic field at the 
muon site we use the function locfield function and specify the 'i' sum
type.

.. literalinclude:: ../examples/MnSi/run_example.py
   :lines: 190-195
   :lineno-start: 190
   :language: python

By default, the contact coupling is 0 for all sites. In order to have a 
contact term different from 0 we have to set the parameter ACont
for all the muon sites that we defined.

.. literalinclude:: ../examples/MnSi/run_example.py
   :lines: 198-200
   :lineno-start: 190
   :language: python

The lines 203-262 produce the following pictures:

|locf|
|hist|

.. |locf| image:: ../examples/MnSi/reference/TotalFields.png
.. |hist| image:: ../examples/MnSi/reference/Histogram.png



Interestingly, left-handed and right-handed orders produce different local field. This is due to the lack of inversion symmetry.

The phase
---------

TODO

Bibliography
------------

.. [Amato2014] Phys. Rev. B 89, 184425
.. [DalmasdeReotier2016] Phys. Rev. B 93, 144419
