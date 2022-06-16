import os
import numpy as np
from warnings import warn
from copy import deepcopy

from muesr.settings import config
from ase.io import read, write
from ase.build import make_supercell

from muesr.core.sampleErrors import MuonError
from muesr.core.parsers import *
from muesr.core.ninput import ninput
from muesr.core.nprint import nprintmsg


def load_xsf(sample, filename):
    """
    Loads structural data from Xcrysden Files in sample object.

    WARNING: this will reset all current muon positions, magnetic
    definitions, and the symmetry definition.

    :param sample: a sample object.
    :param str filename: the filename
    :return: True if succeful, False otherwise
    :rtype: bool
    :raises: TypeError

    """
    try:
        fname = str(filename)
    except:
        raise TypeError

    atoms = read(fname)
    if atoms:
        sample._reset(cell=True,sym=True,magdefs=True,muon=True)
        sample.cell = atoms
        return True
    else:
        nprint ("Atoms not loaded!", 'warn')
        return False


def save_xsf(sample, filename, supercell=[1,1,1], addMuon=True):
    """
    Export structure to XCrysDen.

    :param sample: a sample object.
    :param str filename: path of the destination file.
    :param list supercell: a list containig the number of replica along the three lattice parameters
    :param bool addMuon: if true, adds the muon positions (if any) in the central unit cell
    :return: True if succeful, False otherwise
    :rtype: bool
    :raises: CellError, TypeError
    """


    if len(filename) == 0:
        raise ValueError("Invalid filename.")

    sample._check_lattice()

    sc = make_supercell(sample._cell, np.diag(supercell))

    if (not sc is None) and addMuon:
        try:
            for m in sample.muons:
                spos_m = [ (m[0] + int(supercell[0]/2)) / supercell[0],
                                                (m[1] + int(supercell[1]/2)) / supercell[1],
                                                (m[2] + int(supercell[2]/2)) / supercell[2] ]
                sc.append('H')
                sc.positions[-1] = sc.cell.cartesian_positions(spos_m)
        except MuonError:
            pass

    if not sc is None:
        write(filename,sc)
        return True
    else:
        raise RuntimeError("Cannot build (super)cell for display.")
        return False
