from muesr.settings import config
from muesr.i_o.xsf import save_xsf
import os, subprocess, warnings


def show_structure(sample, supercell=[1,1,1], boundary=[0,0,0] , fields=None, askConfirm=True, block=True, visualizationTool=''):
    """
    Creates a supercell and shows it in XCrysDen or VESTA.
    This function only works interactively.
    WARNING: this can potentially harm your system!

    :param sample: a sample object.
    :param list supercell: a list containig the number of replica along the three lattice parameters
    :param bool askConfirm: for CLI, ask for confirmation before running XCrysDen
    :param bool block: the visualization program is forked if block=False, only works on unix
    :param str  visualizationTool: the preferred application for visualization. Valid values are 'X', 'XCrySDen', 'V' or 'VESTA'.
                      The selection is case insensitive. If missing, any available application will be used.
    :return: True if succeful, False otherwise
    :rtype: bool
    :raises: TypeError
    """

    if supercell[0] <= 0 or supercell[1] <= 0 or supercell[2] <= 0:
        raise ValueError("supercell dimension must be a positive integer.")

    if visualizationTool=='':
        visualizationTool = config.DefaultVisualizationApp

    save_xsf(sample, os.path.join(config.XCrysTmp,'preview.xsf'),supercell=supercell, boundary=boundary, fields=fields)

    if visualizationTool.lower() in ['x', 'xcrysden']:
        return True if run_xcrysden(os.path.join(config.XCrysTmp,'preview.xsf'),block) else \
                        run_vesta(os.path.join(config.XCrysTmp,'preview.xsf'),block)
    elif visualizationTool.lower() in ['v', 'vesta']:
        return True if run_vesta(os.path.join(config.XCrysTmp,'preview.xsf'),block) else \
                        run_xcrysden(os.path.join(config.XCrysTmp,'preview.xsf'),block)
    else:
        warnings.warn("Visualization tool not found.")
        return False

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

def run_vesta(fname, block=True):
    if config.VESTAExec == None:
        warnings.warn("VESTA executable not found. Check configs.")
        return False


    spargs = dict(
        args   = [config.VESTAExec, fname],
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
