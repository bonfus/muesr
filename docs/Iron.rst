BCC Iron
---------

This example shows how to calculate the local field at the tetrahedral 
site(s) in bcc-Fe as described in Ref. [Schmolz1986]_.

Let's first load some useful tools.

.. literalinclude:: ../examples/Fe_bcc/run_example.py
   :lines: 6-14
   :lineno-start: 6
   :language: python

The lattice structure and one of the twelve tetrahedral muon sites are
added to our sample object:

.. literalinclude:: ../examples/Fe_bcc/run_example.py
   :lines: 20-24
   :lineno-start: 20
   :language: python

Since the lattice structure was parsed from a CIF file, symmetry information
are already present in the Sample object. This allows to obtain all the
twelve tetrahedral sites with the command

.. literalinclude:: ../examples/Fe_bcc/run_example.py
   :lines: 29
   :lineno-start: 29
   :language: python

Finally, we define the ferromagnetic structure with local moments parallel
to z in Cartesian coordinates.

.. literalinclude:: ../examples/Fe_bcc/run_example.py
   :lines: 32-37
   :lineno-start: 32
   :language: python

The local fields at the muon sites are obtained with the following commands

.. literalinclude:: ../examples/Fe_bcc/run_example.py
   :lines: 77-78
   :lineno-start: 77
   :language: python

By default, for each muon site the contact coupling is set to 0. To obtain
the total contribution at the muon sites the parameter ACont must be set 
for all the twelse muon sites (can be different in general).

.. literalinclude:: ../examples/Fe_bcc/run_example.py
   :lines: 92-102
   :lineno-start: 92
   :emphasize-lines: 9
   :language: python

In the above lines we collected the results in numpy arrays to simplify the
output statements.


As discussed in [Schmolz1986]_, "The field contribution at each equivalent site
is either parallel or antiparallel to the magnetization of the domains" such that
B_dip(parallel)=-2B_dip(antiparallel) "the average of the dipolar field at these three sites vanishes "

The final print statements shows that this result is actually verified 
by our calculations.

.. literalinclude:: ../examples/Fe_bcc/run_example.py
   :lines: 115-123
   :lineno-start: 115
   :language: python



.. [Schmolz1986] M. Schmolz et.al, Hyperfine Interactions 31 199 (1986)
