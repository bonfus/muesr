# http://magcryst.org/resources/magnetic-coordinates/

from muesr.core.magmodel import MM, have_sympy
if have_sympy:
    from muesr.core.magmodel import SMM

from muesr.core.parsers import *
from muesr.core.nprint  import nprint, nprintmsg, print_cell
from muesr.core.ninput  import *
from muesr.core.cells  import get_cell_parameters


import numpy as np

   

def mago_set_k(sample, kvalue=None, mm=None):
    
    """
    Set propagation with respect to the conventional reciprocal cell
    kavalue : propagation vector (either string or tuple or list), if None, input is prompetd
    mm      : Magnetic Model to be used. If None the current magnetic model is used
    returns: None
    """
    
    
    
    if mm is None:
        smm = sample.mm
    else:
        smm = mm
        
        
    
    if sample._check_lattice():
        if not kvalue is None:
            if isinstance(kvalue, np.ndarray) and ( kvalue.shape == (3,)):
                smm.k = kvalue
                return True
            else:
                raise TypeError("Invalid type for kvalue. Must be 3D vector defined as numpy array.")
        else:
            try:
                if kvalue is None:
                    kval=ninput('Propagation vector (w.r.t. conv. rec. cell): ', parse_vector)
                else:
                    kval = parse_vector(kvalue,3)

            except EOFError:
                nprint("Ok.")
                return False
            except TypeError:
                nprint("Cannot parse position.",'warn')
                return
            smm.k=np.array(kval,dtype=np.float)
            return True


def mago_add(sample, coordinates='b-c', fcs=None, kvalue=None):
    """
    Adds a magnetic model (fourier components and K vector).
    The order is automatically selected if succesfully added.
    
    :param sample: A sample object.
    :param string coordinates: coordinates system and units of the Fourier components.
                               Options are: 'b-c', 'b/a-l', 'b-l'.
                               
                               b-c  : Fourier components in Bohr magnetons and Cartesian coordinates.
                               
                               b/a-l: Fourier components in Bohr magnetons/Angstrom and lattice coordinates.
                               
                               b-l  : Fourier components in Bohr magnetons and in lattice coordinates.
                               
    :param np.complex fcs: Fourier components in coordinate system 
                          (default: Bohr magnetoc/ Cartesian coordinates)
    :param np.ndarray kvalue: Propagation vector in r.l.u.
    :returns: True if successful, False otherwise.
    :rtype: bool
    
    """        
    
    nmm = MM(sample.cell.get_number_of_atoms(),sample.cell.get_cell())
    
    ret = mago_set_k(sample, mm=nmm, kvalue=kvalue)
    if not ret:
        return False
        
    ret = mago_set_FC(sample, mm=nmm, inputConvention=coordinates, fcs=fcs)
    if not ret:
        return False
        
    sample.mm = nmm
    return True

        

def mago_set_FC(sample, fcs = None, atoms_types = None, mm=None, inputConvention='b-c'):
    """
    Defines fourier components for the unit cell.
    
    """


    sample._check_lattice()
    
    
    if mm is None:
        smm = sample.mm
    else:
        smm = mm
    

    
    inputConvEnum = -1
        
    if inputConvention.lower() in ['bohr-cartesian', 'b-c']:
        nprint('Fourier components in Bohr magnetons and Cartesian coordinates.')
        inputConvEnum = 0
        
    elif inputConvention.lower() in ['bohr/angstrom-lattice', 'b/a-l']:
        nprint('Fourier components in Bohr magnetons/Angstrom and lattice coordinates.')
        inputConvEnum = 1
        
    elif inputConvention.lower() in ['bohr-lattice','b-l']:
        nprint('Fourier components in Bohr magnetons and lattice coordinates.')
        inputConvEnum = 2
        
    else:
        nprint('Invalid Fourier Description method: ' + inputConvention.lower() ,'error.')
        return False
    


    
    if not fcs is None:
        smm.fc_set(fcs, inputConvEnum)
        return True
            


    unit_cell = sample._cell    
    lattice = unit_cell.get_cell()
    
    
    if atoms_types is None:
        atoms_types=ninput('Which atom? (enter for all): ').split()
    
    if atoms_types == []:
        set_this_fc = lambda x: True
    else:
        atoms_types = [x.lower() for x in atoms_types] #lowercase
        set_this_fc = lambda x: x.lower() in atoms_types
        # Print where they are
        print_cell(unit_cell,set_this_fc)
    
    

    
    fcs = np.zeros([unit_cell.get_number_of_atoms(),3],dtype=np.complex)
    
    gotEOS = False
    for i, atom in enumerate(unit_cell):                
        try:
            if set_this_fc(atom[0]): # check if we should set this fc
                FCin = ninput_mt("FC for atom %i %s (3 real, [3 imag]): " % (i+1,atom[0]), parse_complex_vector)
                fcs[i]=[FCin[j]+1.j*FCin[j+3] for j in range(3)]
                    
            else:
                # 0 is always 0, no matter the input format! ;)
                fcs[i]=([0.+0.j,0.+0.j,0.+0.j])
                
        except EOFError:
            gotEOS = True
            break

    if gotEOS:
        nprint('Nothing set')
        return False
    else:
        smm.fc_set(fcs, inputConvEnum)
        return True



