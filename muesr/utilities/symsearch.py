import warnings

from muesr.core.spg import spacegroup_from_data
from muesr.core.parsers import *
#from muesr.core.symmetry import Symmetry
from muesr.core.nprint  import nprint, nprintmsg

have_spg = True

try:
    import spglib as spg
except ImportError:
    try:
        from pyspglib import spglib as spg
    except ImportError:
        nprint("Spg Library not loaded", "warn")
        have_spg = False



        
def symsearch(sample, precision=1e-4):
    """
    Identifies symmetry operations of the unit cell using spglib and 
    update the sample definition.
    
    :param sample: A sample object.
    :param float precision: atoms are assumed equivalent if distances are 
                            smaller than precision. In Angstrom.
    :returns: True if succesful, False otherwise.
    :rtype: bool
    """
    
    if sample._check_lattice() and have_spg:
        symbol,number = spg.get_spacegroup(sample.cell, symprec=precision).split()
        operations    = spg.get_symmetry(sample.cell, symprec=precision)
        number = int(''.join([c for c in number if c.isdigit()]))
        sample.sym = spacegroup_from_data(number,symbol,rotations=operations['rotations'],translations=operations['translations']) ## 
        
        warnings.warn("Warning, only works for standard settings!")
        
        return True
    elif not have_spg:
        warnings.warn("Warning, spglib not found. This function will always return False.")
        return False
    else:
        return False
     
