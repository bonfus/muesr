
import numpy as np
import warnings

from copy import deepcopy

from muesr.core.sampleErrors import *
from muesr.core.nprint import cstring
from muesr.core.spg  import Spacegroup
from muesr.core.atoms  import Atoms
from muesr.core.isstr  import isstr
from muesr.core.magmodel  import MM, have_sympy



if have_sympy:
    from muesr.core.magmodel  import SMM




class Sample(object):
    """
    This object contains all the information about the sample under
    investigation.
    It includes a definition of the lattice parameters and the atomic
    positions, the symmetry, the positions of the muons and the magnetic
    structure.
    
    In each sample there can only be a single lattice structure.
    Therefore, also a single symmetry is defined.
    
    Multiple magnetic structures can be defined instead.
    Every time a magnetic structure is defined, it automatically becomes
    the current magnetic structure.
    To select a different magnetic structure the user can use the 
    property :py:attr:`~current_mm_idx`.
    
    >>> yoursample = Sample()  #initialize the sample
    >>> yoursample.name = "Test"
    
    """
    

    __isfrozen = False       # used to prevent adding stuff

    def __init__(self):
        self._name    = "No name"     # description of the sample
        self._muon    = []            # list containing the muon pos (frac. coord)
        self._magdefs = []            # list containing MM objects
        self._sym     = None          # contains symmetry object
        self._cell    = None          # contains an Atoms object
        self._selected_mm = -1        # the selected magnetic structure.        
        self._freeze()
    
    def __setattr__(self, key, value):
        #check if attribute is present
        _hasattr = False
        # this is a work around till I find a better way to check 
        # properties that raises exceptions both in python2 and python3
        if key in ['name', 'mm', 'current_mm_idx','sym','cell']:
            _hasattr = True
        else:
            try:
                _hasattr = hasattr(self, key)
            except SampleException:
                _hasattr = True
        
        if self.__isfrozen and not _hasattr:
            raise TypeError( "Cannot set attributes in this class" )
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True    
    
    @property
    def name(self):
        """
        The sample name.
                
        :getter: Returns the sample name
        :setter: Sets the sample name
        :type: str
        """
        return self._name
        
    @name.setter
    def name(self, value):
        if isstr(value):
            self._name = value
        else:
            raise TypeError('Invalid type for \'name\' property. Must be string.')
        
    @property
    def muons(self):
        """
        Muon positions in fractional coordinates.
        
        :getter: Returns a list of numpy array of shape (3,).
        """
        self._check_muon()
        return deepcopy(self._muon)
        
    def add_muon(self, position, cartesian=False):
        """
        Adds a muon position.
        
        :param list position: list of three float or numpy array of 
                              shape (3,) describing the position of the
                              muon.
        :param bool cartesian: if True, the position os assumed to be in 
                               Cartesian coordinates. If False (default)
                               the position is assumed to be in fractional
                               coordinates.
        :returns: None
        :rtype: None
        :raises: MuonError                       
        """
        # validation of input
        self._check_lattice()
        
        if (type(position) is list):
            position = np.array(position)
        
        if (type(position) is np.ndarray):
            if not position.shape == (3,):
                raise ValueError('Invalid shape for muon position. Must be 3D vector.')
        else:
            raise TypeError('Invalid input for muon position. Must be list of numpy array.')
        
        
        if cartesian:
            #go to reduced lattice coordinates...check this
            self._muon.append(np.dot(position,
                                        np.linalg.inv(self._cell.cell)))
        else:
            self._muon.append(1.0*position) # make floats from ints
    
    
    @property
    def mm_count(self):
        """
        Count the loaded Magnetic Model(s).
        
        :getter: Returns the current Magnetic Model (a MM object)
        :setter: Read-only
        :type: :py:class:`~MM` object
        :raises: None
        """
        
        return len(self._magdefs)
    
    @property
    def mm(self):
        """
        Current Magnetic Model.
        
        :getter: Returns the current Magnetic Model (a MM object)
        :setter: Adds and select a new Magnetic Model (MM object)
        :type: :py:class:`~MM` object
        :raises: TypeError, MagDefError
        """
        self._check_magdefs()
        return self._magdefs[self._selected_mm]
            
    @mm.setter
    def mm(self, value):
        
        self._check_lattice()
        
        if isinstance(value, MM):
            if self._cell.get_number_of_atoms() == len(value.fc):
                self._magdefs.append(value)
                self._selected_mm = len(self._magdefs) - 1
            else:
                raise MagDefError('Number of fourier components does not match number of atoms')
        else:
            raise TypeError("MM or SMM instance expected")
            

    def new_mm(self):
        """
        Creates a new empty magnetic model (everything is set to zero)
        
        :returns: None
        :rtype: None        
        """        
        self._check_lattice()
        self._magdefs.append(MM(self._cell.get_number_of_atoms(), 
                                self._cell.get_cell()))

                                
        self._selected_mm = len(self._magdefs) - 1        
    
    def new_smm(self, symbolic_parameters):
        """
        Creates a new empty Symbolic Magnetic Model (everything is set to zero)
        
        :param str symbolic_parameters: String containing a list of symbolic parameters. Example "x,y,z"
        :returns: None
        :rtype: None   
        :raises: ImportError if sympy is not installed. CellError if lattice is not defined.
        """         
        if not have_sympy:
            raise ImportError('Missing sympy dependency! Install it to use this fuction')
        
        self._check_lattice()       
        self._magdefs.append(SMM(self._cell.get_number_of_atoms(), 
                                number_of_symbolic_parameters, 
                                self._cell.get_cell()))
        
        # select last model
        self._selected_mm = len(self._magdefs) - 1 
            
    @property
    def current_mm_idx(self):
        """
        Index of the currently selected magnetic model (starting from 0)
        
        :getter: returns the index of the currently selected Magnetic Model
        :setter: Sets the current Magnetic Model from the index.
        :type: int
        :raises: MagDefError        
        """
        self._check_magdefs()
        return self._selected_mm
    

    @current_mm_idx.setter
    def current_mm_idx(self, i):
        self._check_magdefs()
        
        if i < len(self._magdefs) and i >= 0:
            self._selected_mm = i
        else:
            raise IndexError('Only %d magnetic structures defined' % len(self._magdefs))
        
    @property
    def sym(self): 
        """
        Symmetry of the system, i.e. a Spacegroup Object.
        
        :getter: returns a Spacegroup Object.
        :setter: Sets the symmetry of the system from a Spacegroup Object.
        :type: int
        :raises: TypeError, SymError
        """
        self._check_sym()
        return self._sym
        
    @sym.setter
    def sym(self,value):
        if not (isinstance(value,Spacegroup)):
            raise TypeError('Symmetry is invalid.')	
        
        self._sym = value
        
    @property
    def cell(self): 
        """
        Returns the atomic structure definition, i.e. an Atoms object.

        :getter: returns a Atoms Object.
        :setter: Sets the atomic structure definition from a Atoms object.
        :type: int
        :raises: TypeError, CellError        
        """
        self._check_lattice()
        return deepcopy(self._cell)
        
    @cell.setter
    def cell(self, value): 
        if not (isinstance(value,Atoms)):
            raise TypeError('Cell is invalid.')
                
        self._cell = value
    
        
    def __repr__(self):
        """
        Print current status of sample definition on print() command.
        """
        
        has_mag= False
        
        status = "Sample status: \n\n"
        
        # Check crystal structure
        status += " Crystal structure:".ljust(30)
        if self._cell is None:
            status += cstring('No','warn')
        else:
            status += cstring('Yes','ok')
        status += '\n'

        # Check magnetic structure
        status += " Magnetic structure:".ljust(30)
        if not self._magdefs:
            status += cstring('No','warn')
        else:
            status += cstring('Yes','ok')
        status += '\n'

        # Check muon position
        status += " Muon position(s):".ljust(30)
        if not self._muon:
            status += cstring('No','warn')
        else:
            status += cstring(str(len(self._muon))+' site(s)','ok')
        status += '\n'

        # Check symmetry
        status += " Symmetry data:".ljust(30)
        if not self._sym:
            status += cstring('No','warn')
        else:
            status += cstring('Yes','ok')
        status += '\n'

        if self._magdefs:
            status += "\nMagnetic orders available ('*' means selected)\n\n"
            status += " Idx | Sel | Desc. \n"
            for i, m in enumerate(self._magdefs):
                fmt_args = (i, \
                            '*' if (i == self._selected_mm) else ' ', \
                            m.desc)
                status += " {:2d}  |  {}  | {}\n".format(*fmt_args)
            
            

        return status
    
    def check_status(self, cell=False, magdefs=False, muon=False, sym=False):
        ok = True
        if cell:
            try:
                self._check_lattice()
            except CellError:
                ok = False
                
        if magdefs:
            try:
                self._check_magdefs()
            except MagDefError:
                ok = False
                
        if muon:
            try:
                self._check_muon()
            except MuonError:
                ok = False

        if sym:
            try:
                self._check_sym()
            except SymmetryError:
                ok = False
                
        return ok
                
                

    def _reset(self,cell=False,magdefs=False,muon=False,sym=False):
        """
        Reset some fields of the sample definitions to initial value.
        """
        if cell:
            self._cell = None
            if not muon:
                warnings.warn("Resetting cell but preserving muon " +
                              "position! Are you sure?", RuntimeWarning)
            if not magdefs:
                warnings.warn("Resetting cell but preserving magnetic" +
                              "structures! Are you sure?", RuntimeWarning)
            if not sym:
                warnings.warn("Resetting cell but preserving symmetry" +
                              "details! Are you sure?", RuntimeWarning)
            
        if magdefs:
            self._magdefs = []
            self._selected_mm = -1
        if muon:
            self._muon = []
        if sym:
            self._sym = None

    def _check_sym(self):
        if self._sym is None:
            raise SymmetryError('Symmetry is not defined.')
        if not (isinstance(self._sym,Spacegroup)):
            raise SymmetryError('Symmetry is invalid.')	
        return True

    def _check_lattice(self):
        if self._cell is None:
            raise CellError('Crystal structure not defined')
        else:
            if (isinstance(self._cell,Atoms)):
                return True
            else:
                raise CellError('Crystal structure type is wrong!')
        
    def _check_magdefs(self):
        if type(self._magdefs) is list :
            if len(self._magdefs) > 0:
                return True
            else:
                raise MagDefError('Magnetic structure not defined')
        else:
            raise MagDefError('Magnetic structure type is wrong!')
        
    def _check_muon(self):
        if type(self._muon) is list:
            if len(self._muon) > 0:
                for e in self._muon:
                    if not(type(e) is np.ndarray):
                        raise MuonError('Muon position not defined or wrong type')
                return True
            else:
                raise MuonError('Muon position not defined')


        
