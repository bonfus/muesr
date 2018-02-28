# http://magcryst.org/resources/magnetic-coordinates/

import numpy as np
from muesr.core.cells import get_cell_parameters
from muesr.core.isstr import isstr

have_sympy = True
try:
    import sympy as sy
except:
    have_sympy = False


class MM(object):
    """
    Magnetic model class.
    
    The magnetic structures are defined in terms of a single propagation
    vector, complex Fourier components and phases.
    
    
    
    :param int cell_size: number of Fourier components=number of atoms
    :param latt_vects: lattice vectors as numpy 3x3 ndarray. 
        If None, the Fourier componentscan only be specified in cartesian coordinates.        
    :raises: TypeError
    """
    
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "Cannot set attributes in this class" )
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True
            
    def __init__(self, cell_size, latt_vects=None):
        
        self._size = int(cell_size)
        
        try:
            assert(self._size>0)
        except:
            ValueError("Cannot parse size of mag model." + 
                      "Must be (strictly) positive int.")

        self._description = "No title"
                
        self._validFormats = [0]
        if isinstance(latt_vects, np.ndarray):
            if latt_vects.shape == (3,3):
                self._latt = latt_vects
                self._rlatt = np.dot(np.diag(np.divide([1.,1.,1.],
                            get_cell_parameters(self._latt))),self._latt)
                self._validFormats += [1,2]
            else:
                raise TypeError("Cannot parse lattice vectors.")
        elif latt_vects is None:
            self._latt = None
            self._rlatt = None
        else:
            raise TypeError("Lattice vectors must be numpy array.")
        
        
        self._fc = np.zeros([cell_size,3],dtype=np.complex)
        self._k = np.array([0,0,0])
        self._phi = np.zeros(cell_size,dtype=np.float)
        
        self._freeze() # no new attributes after this point. 

    @property
    def size(self):
        """
        Number of Fourier components.
        """
        return self._size
        
    @property
    def lattice_params(self):
        """
        Lattice parameters used to convert the Fourier components to the 
        various coordinates systems.
        """        
        return self._latt.copy()

    @property
    def k(self):
        """
        The propagation vector in reciprocal lattice units.
                
        :getter: Returns a numpy array of shape (3,) 
                 with the propagation vector.
        :setter: Sets the propagation vector.
        :type: numpy ndarray of shape (3,)         
        """
        return self._k.copy()

    @k.setter
    def k(self, value):
        if isinstance(value, list):
            try:
                value = np.asarray(value,np.float)
                if value.shape == (3,):
                    self._k=value
                else:
                    raise ValueError('Array must be a single 3D vector.')
            except:
                raise ValueError('Array must be a single 3D vector.')
                
        elif isinstance(value, np.ndarray):
            if value.shape == (3,):
                self._k=np.asarray(value,np.float)
            else:
                raise ValueError('Array must be a single 3D vector.')
        else:
            raise TypeError('k value must be a 3D numpy array or list.')
            
    @property
    def fc(self):
        """
        Fourier components  in Cartesian coordinates. 
        Same as :py:attr:`~fcCart`
        """
        return self.fc_get()
        
    @fc.setter
    def fc(self, value):
        return self.fc_set(value)


    @property
    def fcCart(self):
        """
        Get Fourier components in Cartesian coordinates. 
                
        :getter: Returns a numpy array of size (:py:attr:`~size`,3) 
                 with the fourier components.
        :setter: Sets the fourier compinents
        :type: numpy ndarray of size (:py:attr:`~size`,3)         
        """
        return self.fc_get()

    @fcCart.setter
    def fcCart(self,value):
        return self.fc_set(value)

    @property
    def fcLattBMA(self):
        """
        Get Fourier components in Bohr Magneton/Angstrom units, with x||a, y||b and z||c
        
        :getter: Returns a numpy array of size (:py:attr:`~size`,3) 
                 with the fourier components.
        :setter: Sets the fourier compinents
        :type: numpy ndarray of size (:py:attr:`~size`,3)         
        """
        return self.fc_get(1)

    @fcLattBMA.setter
    def fcLattBMA(self,value):
        return self.fc_set(value, 1)

    @property
    def fcLattBM(self):
        """
        Fourier components in Bohr Magneton units, with x||a, y||b and z||c
        
        :getter: Returns a numpy array of size (:py:attr:`~size`,3) 
                 with the fourier components.
        :setter: Sets the fourier compinents
        :type: numpy ndarray of size (:py:attr:`~size`,3) 
        """
        return self.fc_get(2)

    @fcLattBM.setter
    def fcLattBM(self, value):
        return self.fc_set(value,2)

    def fc_get(self, coord_system=0):
        """
        Retrives Fourier components.
        
        :params int coord_system: requested coordinate system.
                
            * 0 means Cartesian coordinates
            * 1 means Lattice coordinates (values in units of Bohr Magneton/Angstrom)
            * 2 means Bohr magnetons in each lattice direction (values in units Bohr Magneton)
            
            See http://magcryst.org/resources/magnetic-coordinates/    
                
        :raises: ValueError
        """
        
        # check if the lattice is defined, otherwise no conversion :/
        coord_system = int(coord_system)
        if not (coord_system in self._validFormats):
            raise ValueError("Invalid/unsupported input type. Have you provided lattice cell at instantiation?")
        
        if coord_system == 0:
            return np.copy(self._fc)
        elif coord_system == 1:
            return np.copy(np.dot(self._fc,np.linalg.inv(self._latt)))
        elif coord_system == 2:
            return np.copy(np.dot(self._fc, np.linalg.inv(self._rlatt)))
        

        
    
    def fc_set(self, value, coord_system=0):
        """
        Sets Fourier components.
        
        :param value: numpy array containing the fourier components for all the atoms.
        :param int coord_system: requested coordinate system.
                
            * 0 means Cartesian coordinates
            * 1 means Lattice coordinates (values in units of Bohr Magneton/Angstrom)
            * 2 means Bohr magnetons in each lattice direction (values in units Bohr Magneton)
            
            See http://magcryst.org/resources/magnetic-coordinates/    
                
        :raises: ValueError, TypeError 
        """       
        
        validFCS = False
        validType = False
        
        # validate FCs
        if isinstance(value, np.ndarray):
           
            if value.dtype != np.complex:
                raise ValueError("Fourier components must be a complex array!")
            
            #check that number of FCs is the same as atoms
            if not (value.shape == self._fc.shape):
                raise ValueError("Invalid shape: {}".format(value.shape) + " instead of {}".format(self._fc.shape))
        else:
            raise TypeError("Value must be numpy array.")
        
        #validate coord_system
        if type(coord_system) != int:
            raise TypeError("coord_system must be int.")
        
        # check if the lattice is defined, otherwise no conversion :/
        if not (coord_system in self._validFormats):
            raise ValueError("Invalid/unsupported input type. Have you provided lattice cell at instantiation?")
        
        if coord_system==0:
            # no conversion needed
            self._fc = value
            
        elif coord_system==1:
            self._fc = np.dot(value, self._latt)
            
        elif coord_system==2:                           
            # we get cartesian coordinates as for atoms.
            self._fc = np.dot(value, self._rlatt)


    @property
    def phi(self):
        """
        The phase for each fourier component in units of 2 PI.

        :getter: Returns a numpy array of size :py:attr:`~size` with the 
                 phases.
        :setter: Sets the phases. Type can be list of numpy array.
        :type: numpy ndarray        
        """
        return self._phi.copy()

    @phi.setter
    def phi(self, value):
        if isinstance(value, list):
            value = np.array(value,dtype=np.float)
        
        if isinstance(value, np.ndarray):
           
            if value.shape == self._phi.shape:
                self._phi=np.asarray(value,np.float)
            else:
                raise ValueError("Incorrect size of array/list. Must " +
                                  "be a 1D list of " + str(self._size) + 
                                  "phases.")
        else:
            raise TypeError("Must be numpy array or list.")

    @property
    def desc(self):
        """
        Description of the magnetic structure
        
        :getter: Returns the description
        :setter: Sets the description
        :type: str

        """           
        return self._description
    
    @desc.setter
    def desc(self, value):
        if isstr(value):
            try:
                value = str(value)
                self._description = value
            except:
                raise TypeError("Description type must be type string, got {} instead.".format(type(value)))
        else:
            raise TypeError("Description type must be type string, got {} instead.".format(type(value)))
            

    @property
    def isSymbolic(self):
        """
        True if Fourier components are defined as sympy symbols, 
        False otherwise.
        
        ALWAYS False for MM objects.
        """        
        return False


#def set_k(mm, kval):
#    '''
#    Sets propagation vector of a MagneticModel object.
#    '''
#    if (type(kval) is np.ndarray):
#        mm.k = kval
#        return True
#        
#    elif (type(kval) is list):
#        mm.k = np.array(kval)
#        return True
#        
#    else:
#        return False
#    



if have_sympy:
    class SMM(MM):
        """
        This class defines a symbolic magnetic order in which Fuorier 
        components are symbols and not numbers.
        
        :param int cell_size: number of Fourier components=number of atoms
        :param str sparams: comma separated symbols used in sympy.symbols. Example: "x,y,z"
        :param latt_vects: lattice vectors as numpy 3x3 ndarray. 
            If None, the Fourier componentscan only be specified in cartesian coordinates.        
        
        """
        
        
        def __init__(self, cell_size, sparams, latt_vects=None):
            """
            Initialize symbolic magnetic structure.
            

            
            """
            
            self._symbols = sy.symbols(sparams)
            self._symFCexpr = None
            self._symFClambda = None
            self._inputType = -1

            MM.__init__(self, cell_size, latt_vects)
        
        
        def set_symFC(self, value, coord_system=0):
            """
            This function is used to set the symbolic value of the 
            Fourier components for each atom in the system (non magnetic
            atoms should have 0 value)
            
            :param str value: a string passed to sympy.sympify. Must be a 2D list. Example "[[0,0,0],[x,y,0]]" or "[[0,0,z],[I,I,0]]". See sympy documentation for more details
            :param int coord_system: coordinate definition for Fourier components. 
                
                * 0 means Cartesian coordinates
                * 1 means Lattice coordinates (values in units of Bohr Magneton/Angstrom)
                * 2 means Bohr magnetons in each lattice direction (values in units Bohr Magneton)
            
                See http://magcryst.org/resources/magnetic-coordinates/
            
            :return: None
            :rtype: None
            :raises: TypeError, ValueError
            """
            #TODO: check that only declared symbols are present
            
            if type(coord_system) != int:
                raise TypeError("Cannot parse inputType for fourier componets definition.")
            else:
                self._inputType = coord_system

            
            if not isstr(value):
                raise TypeError("value type must be str")
           
            
            
            self._symFCexpr = sy.sympify(value) # "[[0,0,0],[x,y,0]]"
            if self._symFCexpr is None:
                raise ValueError("Cannot sympyfy expression! See sympy manual.")
                
            if len(self._symFCexpr) != self.size:
                size_given = len(self._symFCexpr)
                self._symFCexpr = None
                self._symFClambda = None
                raise ValueError("Invalid shape size: " + str(size_given)
                                    + " instead of " + str(self.size))

                
            self._symFClambda = sy.lambdify(self._symbols,
                                               self._symFCexpr,
                                               modules='numpy')
        
        def set_params(self, values):
            """
            Coverts the symbolic values into numerical values.
            The order of the parameters must be the same as the one given
            in the class instantiaion.
            
            :param list values: a list containing the values for the sympy symbols.            
            :return: None
            :rtype: None
            :raises: RuntimeError           
            """


            
            if self._symFClambda != None:
                self.fc_set( np.array(self._symFClambda(*values), \
                                      dtype=np.complex), \
                             self._inputType)
            else:
                raise RuntimeError("Symbolic FC not defined!")

        @property
        def isSymbolic(self):
            """
            True if Fourier components are defined as sympy symbols, 
            False otherwise.
            
            ALWAYS True for SMM objects.
            """
            return True
    

