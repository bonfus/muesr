import os
import numpy as np
from muesr.core.parsers  import parse_bool
from muesr.core.nprint   import nprint
from muesr.core.cells    import get_cell_parameters,get_angles
from muesr.core.atoms    import Atoms
from muesr.core.spg      import spacegroup_from_data
from muesr.core.magmodel import MM



def export_fpstudio(sample, filename):
    """
    Exports magnetic structure in FullProfStudio format.
    
    :param sample: the sample object
    :param str filename: the filename for the FP Studio file
    :return: None.
    :rtype: None

    
    """
    
    
    outbuffer = ""
    
    if sample._check_lattice():
        lattice = sample._cell.get_cell()
        
        outbuffer += "CELL "
        outbuffer += " ".join([str(x) for x in get_cell_parameters( lattice )])
        outbuffer += " "
        outbuffer += " ".join([str(x) for x in get_angles( lattice )])
        outbuffer += "\n"
        
        for i, atom in enumerate(sample._cell):
            outbuffer += "ATOM " #example Ce1 Ce 0 0 0
            outbuffer += atom[0]+str(i) + " "
            outbuffer += atom[0] + " "
            outbuffer += " ".join([str(x) for x in atom[2]])
            outbuffer += "\n"
        

        outbuffer += "SPACEG P 1\n"

    if sample._check_magdefs():
        outbuffer += "{\n"
        outbuffer += "K " + " ".join([str(x) for x in sample.mm.k]) + "\n"
        outbuffer += "LATTICE P\n"
        outbuffer += "SYMM x,y,z\n"
        outbuffer += "MSYM u,v,w,0.0\n"
        for i, atom in enumerate(sample.cell):
            outbuffer += "MATOM "
            outbuffer += atom[0]+str(i) + " "
            outbuffer += atom[0] + " "  
            outbuffer += " ".join([str(x) for x in atom[2]])                              
            outbuffer += " SCALE 0.8 ENVELOP\n"
            
            outbuffer += "skp 1 1 "
            outbuffer += " ".join([str(x) for x in sample.mm.fcLattBM[i].real]) + " "
            outbuffer += " ".join([str(x) for x in sample.mm.fcLattBM[i].imag]) + " "
            outbuffer += str(sample.mm.phi[i]) + "\n"
        
        

        outbuffer += "}\n"
        
    # TODO: check overwrite
    with open(filename,'w') as f:
        f.write(outbuffer)
            
        
        
        
        
