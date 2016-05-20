Contact Hyperfine Fields - Some notes
=====================================

Muesr assumes that hyperfine fields are isotropic. The contact term is obtained by specifying two terms:
1) the number of nearest neighbors magnetic atoms to the muon
2) a rcont parameters which governs the maximum radius of the interaction.

Finally, the LocalField object has a ACont property which is an effective hyperfine contact coupling term.
The total field is therefore obtained as

.. math::

   \mathbf{B_T} = \mathbf{B_D} + \mathbf{B_L} + ACont \cdot \frac{2 \mu _0}{3} \sum _i ^N \frac{r _i ^{-3} }{\sum _i ^N r_i ^{-3}}  \mathbf{m}_i

The value of N can strongly impact on the results. The best approximation strongly depend on the simulated system and must be considered case by case.


From CGS to SI
--------------

Hyperfine couplings are often reported in mol/emu while Muesr uses :math:`{\buildrel _{\circ} \over {\mathrm{A}}}^{-3}`.
Here's how to convert the former into the latter.

.. math::

   B_c = -\frac{8}{3} \pi  |\Psi (0)|^2  \boldsymbol{\mu}_e \qquad \mbox{(c.g.i)} \\
   B_c = -\frac{2}{3} \mu_0  |\Psi (0)|^2 \boldsymbol{\mu}_e \qquad \mbox{(S.I.)}
   

Assuming the following formula for the hyperfine contact field produced by the nearest neighboring magnetic atom

.. math::

   B_c = N_A A_{cont} \boldsymbol{\mu}_e
   
where :math:`A_{cont}` is expressed in mol/emu and  :math:`\boldsymbol{\mu}_e` is in emu.
As a consequence:
    
.. math::

   A_{cont} = \frac{8 \pi |\Psi (0)|^2}{3 N_A}

Assuming 1 mol/emu, the value for A in :math:`{\buildrel _{\circ} \over {\mathrm{A}}}^{-3}` is

.. math::

   A_{cont} = \frac{1 \mathrm{mol/emu} * 3 N_A}{8 \pi} = 7.188E22 cm^{-3} = 0.071884019 {\buildrel _{\circ} \over {\mathrm{A}}}^{-3}
