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
    :returns: True if successful, False otherwise.
    :rtype: bool
    """
    
    if sample._check_lattice() and have_spg:
        # This may not identify rotation and translations correctly
        #sample.sym = spg.get_spacegroup(sample.cell, symprec=precision)
        
        dataset = spg.get_symmetry_dataset(sample.cell, symprec=precision)
        # spglib 1.9.9 returns a long int in python2, this causes problems
        # for spg.py, simple workaround is to convert it even if not necessary
        sample.sym = spacegroup_from_data(no=int(dataset['number']),
                                          rotations=dataset['rotations'],
                                          translations=dataset['translations'],
                                          subtrans=dataset['origin_shift'])
        
        warnings.warn("Warning, information regarding spacegroup setting might be wrong!",RuntimeWarning, stacklevel=0)
        if sample.sym:
            return True
        else:
            return False
    else:
        return False
     
