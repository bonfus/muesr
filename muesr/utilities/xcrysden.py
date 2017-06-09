from muesr.settings import config
from muesr.i_o.xsf import save_xsf
import os, subprocess, warnings


def show_structure(sample, supercell=[1,1,1], askConfirm=True, block=True):
    """
    Creates a supercell and shows it in XCrysDen.
    This function only works interactively.
    WARNING: this can potentially harm your system!
        
    :param sample: a sample object.
    :param list supercell: a list containig the number of replica along the three lattice parameters
    :param bool askConfirm: for CLI, ask for confirmation before running XCrysDen
    :param bool block: XCrysDen is forked if block=False, only works on unix
    :return: True if succeful, False otherwise 
    :rtype: bool
    :raises: TypeError        
    """
    
    if type(supercell) != list:
        raise ValueError("Supercell must be list")
    if type(supercell[0]) != int or \
        type(supercell[1]) != int or \
        type(supercell[2]) != int:
        warnings.warn("supercell parameters must be integers. Conversion forced.")
        supercell[0] = int(supercell[0])
        supercell[1] = int(supercell[1])
        supercell[2] = int(supercell[2])
        
    if supercell[0] <= 0 or supercell[1] <= 0 or supercell[2] <= 0:
        raise ValueError("supercell dimension must be a positive integer.")
    
    ans = True
    if askConfirm and supercell[0]>5 or \
        supercell[1]>5 or supercell[2]>5:
        try:
            ans = ninput('Are you sure?! [y/N]', parse_bool)
        except EOFError:
            return False
    
    if ans:
        save_xsf(sample, os.path.join(config.XCrysTmp,'preview.xsf'),supercell=supercell)
        return run_xcrysden(os.path.join(config.XCrysTmp,'preview.xsf'),block)
    else:
        return ans

def run_xcrysden(fname, block=True):
    if config.XCrysExec == None:
        warnings.warn("XCrysDen executable not found. Check configs.")
        return False
    
    
    spargs = dict(
        args   = [config.XCrysExec, "--xsf", fname], 
        stdout = subprocess.PIPE, 
        stderr = subprocess.PIPE
    )
    
    if not block:
        if os.name == 'posix':
            spargs['preexec_fn'] = os.setpgrp
        elif os.name == 'nt':
            spargs['creationflags'] = subprocess.CREATE_NEW_PROCESS_GROUP
        
        
    p = subprocess.Popen(**spargs)
    
    if block:
        out, err = p.communicate()
    return True
