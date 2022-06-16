"""
This file was adapted from ASE.
Module to read and write atoms in cif file format.

See http://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for a
description of the file format.  STAR extensions as save frames,
global blocks, nested loops and multi-data values are not supported.
"""

import re, os
import shlex
import warnings

import numpy as np

from muesr.core.magmodel import MM
#from muesr.core.spg import Spacegroup
from muesr.core.nprint import nprint, nprintmsg
from muesr.core.isstr import isstr

from ase.spacegroup import get_spacegroup

#from muesr.i_o.cif.crystal import crystal
from ase.spacegroup import crystal
from ase.spacegroup.spacegroup import spacegroup_from_data
#from muesr.i_o.cif.cell import cellpar_to_cell
from ase.io import read

# Old conventions:
old_spacegroup_names = {'Abm2': 'Aem2',
                        'Aba2': 'Aea2',
                        'Cmca': 'Cmce',
                        'Cmma': 'Cmme',
                        'Ccca': 'Ccc1'}


def load_cif(sample, filename, reset_muon=True, reset_sym=True):
    """

    Loads the structural information from a Cif file.


    :param muesr.core.sample.Sample sample: the sample object.
    :param str filename:   the cif file path (path + filename).
    :param bool reset_muon: if true the muon positions is reinitialized.
    :param bool reset_sym:  if true the symmetry is reinitialized. ALWAYS TRUE
    :returns: True if succesfull, false otherwise
    :rtype: bool
    """

    filename = str(filename)

    filename = os.path.expanduser(filename)

        #nprint ("Problems with file name...file not found!", 'warn')
        #return False

    #print(read_cif(f,0))
    atoms = read(filename) # selectd index 0
    #f.close()
    if atoms:
        sample._reset(cell=True,sym=True,magdefs=True,muon=reset_muon)
        sample.cell = atoms
        #try:
        sample.sym = get_spacegroup(atoms)
        #except:
        #    nprint ("Symmetry not loaded!", 'warn')
        return True
    else:
        nprint ("Atoms not loaded!", 'warn')
        return False

def load_mcif(sample, filename, reset_muon=True, reset_sym=True):
    """
    Loads both the crystalline structure and the magnetic order from a
    mcif file.
    N.B.: This function is EXPERIMENTAL.

    .. note::
       Only the first lattice and magnetic structure is loaded.
       mCif files hosting multiple structures are not supported.

    :param muesr.core.sample.Sample sample: the sample object.
    :param str filename:   the mcif file path (path + filename).
    :param bool reset_muon: if true the muon positions is reinitialized.
    :param bool reset_sym:  if true the symmetry is reinitialized.
    :returns: True if succesfull, false otherwise
    :rtype: bool
    """
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor

    s = Structure.from_file(filename)

    sample._reset(muon=reset_muon,sym=reset_sym)
    # symmetry is already introduced when parsing magnetism
    sample.cell=AseAtomsAdaptor.get_atoms(s)

    # initialization needs the number of atoms in the unit cell
    nmm=MM(len(s),sample._cell.cell.array)
    # magnetic moments are specified in cartesian coordinates.
    # position in crystal coordinates...not nice but simple!
    nmm.fc_set(np.array([x.moment for x in s.site_properties['magmom']], dtype=np.complex))
    # propagation vector is 0 since the cell "contains" the magnetic
    #   structure, no incommensurate structures, sorry!
    nmm.k=np.array([0.,0.,0.])
    sample.mm=nmm

    return True
