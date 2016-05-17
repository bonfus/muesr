Usage
=====

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


