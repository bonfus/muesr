Symmetry
========

MuESR tries to grab symmetry information from either the CIF files or 
using the `spglib` routines.
If a magnetic cif (*.mcif) is loaded, the symmetry is disregarded. 
The symmetry equivalent sites depend indeed only on the symmetry of the 
parent cell. If one consider the symmetry of the magnetic structure and 
searches for the equivalent muon sites, many of them will be missing.
The symmetry of every cell can always be identified with the help of 
the `sym_search` command.
Note however that, when dealing with supercells, some of the symmetry 
equivalent muon positions must be specified by hand.
