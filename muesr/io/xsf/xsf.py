import os
import numpy as np
from copy import deepcopy

from muesr.settings import config

from muesr.io.xsf.xsfio import *
from muesr.io.xsf.xsfrun import *

from muesr.core.parsers import *
from muesr.core.ninput import ninput
from muesr.core.nprint import nprintmsg

from muesr.core.cells import get_simple_supercell


    
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
        
    atoms = read_xsf_file(fname)
    if atoms:
        sample._reset(cell=True,sym=True,magdefs=True,muon=True)
        sample.cell = atoms
        return True
    else:
        nprint ("Atoms not loaded!", 'warn')
        return False
            

    
def show_supercell(sample, supercell=[1,1,1]):
    """
    Creates a supercell and shows it in XCrysDen.
    This function only works interactively.
    WARNING: this can potentially harm your system!
        
    :param sample: a sample object.
    :param list supercell: a list containig the number of replica along the three lattice parameters
    :return: True if succeful, False otherwise 
    :rtype: bool
    :raises: TypeError        
    """
    
    
    
    
    try:
        ans = ninput('Are you sure?! [y/N]', parse_bool)
    except EOFError:
        return False
        
        
    if ans:
        sc = get_simple_supercell(sample, supercell)

        if not sc is None:
            write_xsf(os.path.join(config.XCrysTmp,'preview.xsf'),sc)
            run_xcrysden(os.path.join(config.XCrysTmp,'preview.xsf'))
            return True
        else:
            nprintmsg('esupcell')
            return False
            
        
def show_cell(sample):
    """
    Open XCrysDen with specified structure.
        
    :param sample: a sample object.
    :return: True if succeful, False otherwise 
    :rtype: bool
    """
    
    
    cell = sample.cell
    
    if not cell is None:
        for m in sample.muons:
            cell.extend(symbol="mu",scaled_position=m)
            
        write_xsf(os.path.join(config.XCrysTmp,'preview.xsf'),cell)
        run_xcrysden(os.path.join(config.XCrysTmp,'preview.xsf'))
        return True
        
    else:
        nprintmsg('ecrystal')
        return False
        
