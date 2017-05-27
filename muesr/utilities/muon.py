from muesr.core.parsers import *
from muesr.core.nprint  import nprint, nprintmsg
from muesr.core.ninput  import ninput

import numpy as np


        
def muon_set_frac(sample, arg = None):
    """
    Set muon position in lattice coordinates.
    
    
    returns: None
    """
    if sample._check_lattice():
        
        try:
            if arg is None:
                pos=ninput('Muon position (frac coord): ', parse_vector)
            else:
                if isinstance(arg, np.ndarray) and ( arg.shape == (3,)):
                    pos = arg
                else:
                    pos = parse_vector(arg,3)
                    
        except EOFError:
            nprint("Ok.")
            return False
        except TypeError:
            nprint("Cannot parse position.",'warn')
            return False
            
        if (not isinstance(arg, np.ndarray)):
            pos = np.array(pos)
        
        sample.add_muon(pos)
        return True
    else:
        #nprintmsg('NS')
        return False
        
        
    
    
def muon_find_equiv(sample, eps=1.e-3):
    """
    Given the unit cell symmetry, finds the equivalent muon sites.
    Magnetic order is NOT considered
    :params: eps, number of decimal figures for comparison
    """
    sample._check_sym()
    sample._check_muon()
    
    eqpoints = sample._sym.equivalent_sites(np.array(sample.muons), \
                                                symprec=eps, \
                                                onduplicates='warn')[0]
    sample._reset(muon=True)
    for p in eqpoints:
        sample.add_muon(np.array(p))
    return True
   
def muon_reset(sample):
    """
    Removes all previously set muon positions.
    returns: True
    """
    sample._reset(muon=True)
    return True
