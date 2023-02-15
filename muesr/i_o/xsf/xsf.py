import os
import numpy as np
from warnings import warn
from copy import deepcopy

from muesr.settings import config
from ase.io import read, write
from ase.atoms import Atoms

from muesr.core.sampleErrors import MuonError, MagDefError
from muesr.core.parsers import *
from muesr.core.ninput import ninput
from muesr.core.nprint import nprintmsg


from ase.calculators.calculator import Calculator, all_changes


class Fake(Calculator):

    implemented_properties = ['energy', 'forces']
    nolabel = True

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)


    def calculate(
        self,
        atoms=None,
        properties=None,
        system_changes=all_changes,
    ):
        if properties is None:
            properties = self.implemented_properties

        Calculator.calculate(self, atoms, properties, system_changes)

        self.results['forces'] = atoms.get_initial_magnetic_moments()
        self.results['energy'] = 0


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


def save_xsf(sample, filename, supercell=[1,1,1], boundary=[0,0,0], fields = None):
    """
    Export structure to XCrysDen.

    :param sample: a sample object.
    :param str filename: path of the destination file.
    :param list supercell: a list containig the number of replica along the three lattice parameters
    :return: True if succeful, False otherwise
    :rtype: bool
    :raises: CellError, TypeError
    """
    from muesr.engines.clfc import moments

    if len(filename) == 0:
        raise ValueError("Invalid filename.")

    sample._check_lattice()

    sc_array, p, sc_symbols, m = moments(sample, supercell, rsc=boundary, muon_local_field=fields)
    sc = Atoms(cell=sc_array,pbc=True, magmoms=m, positions=p,symbols=sc_symbols)

    # Convert moments to forces
    C = Fake()
    C.calculate(atoms=sc)
    sc.set_calculator(C)

    if not sc is None:
        write(filename,sc)
        return True
    else:
        raise RuntimeError("Cannot build (super)cell for display.")
        return False
