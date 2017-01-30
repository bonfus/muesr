Basic theory
====================

Local fields in MuSR
---------------------

Muon spin rotation and relaxation spectroscopy is mainly used to probe 
magnetic materials.
We briefly describe here the interactions involved between the muon and 
the electrons of the hosting system which produce a local magnetic field
at the muon site in magnetically ordered samples.

Dipolar Field
+++++++++++++

The dipolar field is produced by the magnetic dipolar interaction between
the polarized electronic orbitals and the muon spin.
Even though the interaction is best described with quantum mechanics, 
for the sake of simplicity, here we approximate the polarized electronic
orbitals with classical dipoles centered at the nuclei. This 
approximation is also implicit in the code and works rather well in many
cases.

The dipolar field is given by

.. math::

   \mathbf{B_{\mathrm{dip}}^\prime} = \frac{\mu_0}{4 \pi} \sum _{i=1} ^N \left( -\frac{\mathbf{m}}{r^3 _i} + \frac{3 (\mathbf{m}_i \cdot \mathbf{r}_i)\mathbf{r}_i }{r^5 _i} \right)

where, from a quantum perspective, :math:`\mathbf{m}_i = -g_i \mu_\mathrm{B} \mathbf{J}_i`
and :math:`\mathbf{J}_i` is the total angular momentum of the i-th atom.
Finally, the radius :math:`\mathbf{r}_i` is the distance between the muon
and the i-th atom of N magnetic ions in the sample.

When the above sum is performed in real space, it is customary to 
select a spherical portion of the sample (smaller than a magnetic domain)
centered at the muon site and subdivide :math:`\mathbf{B_{\mathrm{dip}}}` in
three terms:

.. math::

   \mathbf{B_{\mathrm{dip}}^\prime} = \mathbf{B_{\mathrm{dip}}} + \mathbf{B_{\mathrm{Lor}}} + \mathbf{B_{\mathrm{dem}}}

The first term originates from the magnetic moments inside the sphere of 
radius :math:`R_\mathrm{sphere}`, i.e.:

.. math::

   \mathbf{B_{\mathrm{dip}}} = \frac{\mu_0}{4 \pi} \sum _{r_i<R_\mathrm{sphere}} \left( -\frac{\mathbf{m}}{r^3 _i} + \frac{3 (\mathbf{m}_i \cdot \mathbf{r}_i)\mathbf{r}_i }{r^5 _i} \right)


The second and the term originate from magnetic moments outside the
sphere and are evaluated in the continuum approximation.
They are

.. math::

   \mathbf{B_{\mathrm{Lor}}} = \frac{\mu_0}{3} \mathbf{M}_{\mathrm{Lor}} = \frac{\mu_0}{3 V_\mathrm{sphere}} \sum _{r_i < R_\mathrm{sphere}} \mathbf{m}_i
   

.. math::

    \mathbf{B_{\mathrm{dem}}} = - \mu_0 \mathbf{N} \mathbf{M}_\mathrm{meas}
    
where :math:`\mathbf{N}` is the demagnetization tensor and :math:`\mathbf{M}_\mathrm{meas}`
is the **bulk** magnetization of the sample. 


.. note::
  :py:mod:`muesr` only estimates :math:`\mathbf{B}_\mathrm{dip}` and 
  :math:`\mathbf{B}_\mathrm{Lor}`.
  The demagnetisation field depends on both the sample details and the 
  experiment details and must be evaluated case by case.



Contact Hyperfine field
+++++++++++++++++++++++


There is another source of local magnetic field at the muon site
which is referred to as Fermi contact hyperfine field.
It originates from the direct interaction between the muon and polarized
electrons at the muon site.
For a polarized spherical electronic cloud surrounding the muon one has

.. math::

   \mathbf{B_{\mathrm{cont}}} = \frac{2 \mu_0}{3} \vert \psi_s (\mathbf{r}_\mu) \vert ^2 \mathbf{m}_e ^s
   
In :py:mod:`muesr`, only a scalar relation between :math:`\mathbf{B_{\mathrm{cont}}}` and 
:math:`\mathbf{m}_e` is allowed and is expressed as :math:`\vert \psi_s (\mathbf{r}_\mu) \vert ^2`.

There is another important point which strongly affects the hyperfine 
field results: the number of nearest neighbours considered in the above sum.
The importance of this term is a direct consequence of the strong 
approximations that we are introducing the the current version of :py:mod:`muesr`.
The contact hyperfine interaction is a purely quantistic phenomenon and
an accurate description would require the knowledge of the electronic
distribution at the muon site.
This is **very badly** approximated by considerig that each magnetic 
atom while contribute to the total hyperfine field by an amount which is
inversely propostional to the cube of its distance from the muon. The 
total is then scaled by the facotr ACont.

[TODO]

Improve discussion about effective nature of the contact term used in muesr!!!!


.. _intro_description_of_magnetic_structures:

Description of Magnetic Structures
-----------------------------------

There are two possibilities to describe a magnetic structure: using the
colored group theory or with the propagation vector and Fourier 
coefficients formalism. :mod:`muesr` opts for the latter.
A magnetic structure is defined as

.. math::

   \mathbf{\mu_{n \nu}} = \sum _{\mathbf{k}} \mathbf{m}_{\nu \mathbf{k}} e ^{- 2 \pi i \mathbf{k} \cdot \mathbf{R}_n}
   
where :math:`\nu` runs over the atoms of the unit cell and :math:`n` 
identifies the n-th cell where atomic positions :math:`\mathbf{R}_{n\nu}` 
are obtained according to

.. math::

   \mathbf{R}_{n\nu} = \mathbf{R}_{n} + \mathbf{r_\nu}
   
with :math:`\mathbf{R}_{n} = n_a \mathbf{a} + n_b \mathbf{b} + n_c \mathbf{c}` 
and :math:`\mathbf{r}_\nu = x_\nu \mathbf{a} + y_\nu \mathbf{b} + z_\nu \mathbf{c}`.

The fourier coefficients :math:`m_{\nu \mathbf{k}}` are three dimensional
complex vectors. They are related to the  irreducible representations 
of the so called "little groups" i.e. the subgroup of the crystallographic space 
group formed by the operators leaving invariant the propagation vector.



[TODO] Discuss the phase!


:mod:`muesr` can only handle 1-k magnetic structures.
However, since local field are linear in the magnetic moment, the
results for multiple-k magnetic orders can be obtained by performing 
multiple simulations for each of the k vectors and Fourier components
which describe the system and summing the results.

Implementation details
----------------------------

:mod:`muesr` is a tool to analyze muon sites and local field contributions
generated by a known magnetic structure. It is intended to be used in an 
interactive python environment such as `IPython <http://ipython.org>`_ or `Jupyter <http://jupyter.org>`_ notebooks.

Internally, muesr uses Tesla units and Angstrom for lengths if not 
specified. Magnetic moments are specified in units of Bohr magnetons.


