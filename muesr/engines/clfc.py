import os
import numpy as np
from copy import deepcopy

from muesr.core.sample import Sample
from muesr.core.occupations import Occupations
from muesr.core.isstr import isstr

import lfclib as lfcext
try:
    from LFN import Simple
except:
    print("Numba version not available.")
    pass

class LocalFields(object):
    """
    Object containing local field components.

    The three local field components are Lorentz field, Dipolar field and Contact hyperfine field.
    These are stored in numpy ndarrays in Tesla units.
    The arrays can have shape (3,) for simulations of a single angle
    (i.e. without rotation of the local moments) or (n,3) for simulations where local moments are rotated n times.

    The total field is the sum of

    .. math::

        \\mathbf{B}_L + \\mathbf{B}_D + \\mathrm{ACont} \\cdot \\frac{2 \\mu_0}{3} \\sum_{r_i < r_{cont}} W(r_i) m

    where  :math:`m_i` are the local magnetic moments and :math:`r_{cont}` is defined during the simulation and cannot be changed.
    :math:`W(r)` is a weighting factor defined as

    .. math::

        W(r) = \\frac{1/r^3}{\\sum_{r < r_{cont}} 1/r_i^3}

    where  :math:`r_i` are the distances of the local moments from the muon.
    By default :math:`\\mathrm{ACont} = 0` so the contact field is zero.
    To introduce the contact term the user must modify :math:`\\mathrm{ACont} = 0`
    using the property :py:attr:`~ACont`.

    The object is initialized as LocalFields(BCont, BDip, BLor, ACont=0.).

    """

    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "Cannot set attributes in this class" )
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True

    def __init__(self, BCont, BDip, BLor, ACont=0., nMuons=1, reshape=True):

        try:
            assert(type(BLor) is np.ndarray)
            assert(type(BDip) is np.ndarray)
            assert(type(BCont) is np.ndarray)
        except AssertionError:
            raise TypeError("Must be numpy arrays!")

        try:
            assert (BDip.shape == BCont.shape == BLor.shape)
        except AssertionError:
            raise ValueError("Must have the same shape!")

        nAngs = np.size(BDip) / nMuons / 3
        if nAngs.is_integer():
            nAngs = int(nAngs)
        else:
            raise ValueError("Something went horribly wrong!")

        self._BLor = np.asarray(BLor,np.float).reshape(nMuons, nAngs, 3)
        self._BDip = np.asarray(BDip,np.float).reshape(nMuons, nAngs, 3)
        self._BCont = np.asarray(BCont,np.float).reshape(nMuons, nAngs, 3)

        try:
            self._ACont = np.asarray(ACont,np.float)
        except:
            raise TypeError( "Cannot set value for ACont. Must be float." )

        if (nMuons <= 0):
            raise ValueError("At least one muon site!")

        if (reshape):
            self._BLor  = np.squeeze(self._BLor)
            self._BDip  = np.squeeze(self._BDip)
            self._BCont = np.squeeze(self._BCont)

        self._nMuons = nMuons
        self._freeze() # no new attributes after this point.

    def __repr__(self):
        return (self._BLor + self._BDip + self._scale(self._ACont,self._BCont)).__repr__()

    def _scale(self, alpha, V):
        return (V.T * alpha).T

    @property
    def L(self):
        """
        Same as :py:attr:`~Lorentz`
        """
        return self._BLor

    @property
    def Lorentz(self):
        """
        Lorentz field in Tesla units.

        :getter: Returns a numpy ndarray containing the Lorentz field.
        """
        return self._BLor

    @property
    def D(self):
        """
        Same as :py:attr:`~Dipolar`
        """
        return self._BDip

    @property
    def Dipolar(self):
        """
        Dipolar field in Tesla units.

        :getter: Returns a numpy ndarray containing the dipolar field.
        """
        return self._BDip

    @property
    def C(self):
        """
        Same as :py:attr:`~Contact`
        """
        return self._scale(self._ACont,self._BCont)

    @property
    def Contact(self):
        """
        Contact hyperfine field in Tesla units.

        :getter: Returns a numpy ndarray containing the contact hyperfine field i.e. :math:`\\mathbf{B}_{Cont} = \\mathrm{ACont} \\cdot \\frac{2 \\mu_0}{3} \\sum_{r_{cont}} m`
        """
        return self._scale(self._ACont,self._BCont)

    @property
    def T(self):
        """
        Same as :py:attr:`~Total`
        """
        return self._BLor + self._BDip + self._scale(self._ACont,self._BCont)

    @property
    def Total(self):
        """
        Total field in Tesla units.

        :getter: Returns a numpy ndarray containing the total field i.e. :math:`\\mathbf{B}_L + \\mathbf{B}_D + \\mathbf{B}_{Cont}`
        """
        return self._BLor + self._BDip + self._scale(self._ACont,self._BCont)

    @property
    def ACont(self):
        """
        Contact hyperfine coupling. Units are :math:`Angstrom^{-3}`

        :getter: Returns the coupling constant
        :setter: Sets the coupling constant
        :type: float
        """
        return self._ACont

    @ACont.setter
    def ACont(self,value):

        try:
            v = np.asarray(value)
        except:
            raise TypeError( "Cannot set value for ACont" )
        if (self._nMuons == 1):
            if not (v.shape == () or v.shape == (1,)):
                raise ValueError("ACont must be set for "+str(self._nMuons)+ " muons")
        else:
            if not (v.shape == (self._nMuons,)):
                raise ValueError("ACont must be set for "+str(self._nMuons)+ " muons")

        self._ACont = v


def find_largest_sphere(sample, supercell):
    """
    Simple function to evaluate the shortest distance between the
    muon position in the supercell and one of the lattice planes
    defined by the lattice coordinates.

    :param sample: the sample object
    :param list supercell: the size of the supercell along the lattice coordinates.
    :return: the radius of the biggest sphere that can be inscribed or None.
    :rtype: float
    :raises: ValueError, TypeError
    """

    # http://mathworld.wolfram.com/Point-PlaneDistance.html
    # https://en.wikipedia.org/wiki/Plane_%28geometry%29

    # check current status is ok
    sample._check_lattice()
    sample._check_muon()

    if (type(supercell) is list):
        if len(supercell) != 3:
            raise ValueError("Wrong supercell definition")

    elif (type(supercell) is np.ndarray):
        if supercell.spahe != (3,):
            raise ValueError("Wrong supercell definition")
    else:
        raise TypeError("Argument supercell must be list of numpy array")

    if np.min(supercell) < 1:
        raise ValueError("Unit cell repetitions must be strictly positive!")

    # build supercell, only diagonal
    cell = sample._cell.cell.array
    scell = np.dot(cell,np.diag(supercell))


    for mu in sample.muons:

        distances = []

        # muon position in Cartesian coordnates
        mup = np.dot((mu + np.floor(np.array(supercell,dtype=np.float64)/2.))/np.array(supercell,dtype=np.float64), scell)
        # define planes using point and normal vector

        for i,j,k in [[1,1,0],[1,0,1],[0,1,1]]:

            p1 = float(i)*scell[0,:]
            p2 = float(j)*scell[1,:]
            p3 = float(k)*scell[2,:]
            # parallel plane shifted by one lattice vector
            sp = float(1-i)*scell[0,:] + float(1-j)*scell[1,:] + float(1-k)*scell[2,:]


            n = np.cross((p2-p1),(p3-p1))
            n /= np.linalg.norm(n)
            # r = p1 , is included in the formula below.

            D = np.dot(n,mup-p1)
            distances.append(np.abs(D))

            # plane shifted by one lattice vector
            D = np.dot(n,mup-p1-sp)
            distances.append(np.abs(D))

    #nprint("WARNING: this is and experimental function!",'warn')
    return np.min(distances)

def locfield(sample, ctype, supercellsize, radius, nnn = 2, rcont = 10.0, nangles = 1, axis = np.zeros(3), constraints=None, backend='available'):
    """
    Evaluates local fields at the muon site.

    This function gives access to three types of calculations specfified
    with the option `ctype`:

        * 'sum': magnetic moments are just summed up to the Lorentz sphere
        * 'rotate': magnetic moments are rotate `nangles` times around the `axis` axis.
        * 'incommensurate': this function is particulary useful for
          incommensurate structures. It performs the sum with the method discussed in PRB XX XXXXXX


    :param sample: the sample object
    :param str ctype: calculation type. Can be 'sum', 'rotate' or 'incommensurate' (or abbreviations 's', 'r', 'i').
    :param list supercellsize: the size of the supercell along the lattice coordinates.
    :param float radius: the radius of the sphere used to evaluate the dipolar tensor.
    :param int nnn: number of local moments nearest neighbours of the muon considered for the contact hyperfine field estimation. Default 2.
    :param float rcont: maximum radius used to search for local moments close to the muon in the contact hyperfine field estimation in Angstrom. Default 10 Angstrom.
    :param int nangles: for 'rotate' and 'incommensurate' simulations, a nangles number of  estimation will perfomed on local moments incrementally rotated by 360/nangles.
    :param list axis: for 'rotate' simulations, axis used to perform the rotation. In 'incommensurate' simulations the axis is defined as the perpendicular vector to the real and the imaginary parts of the fourier componts (warnings will be printed if this vector is not well defined).
    :return: a list of :py:class:`~LocalFields` containing the local field components for each muon site defined in the sample.
    :rtype: list
    :raises: TypeError, ValueError

    """

    # check sample is a Sample object
    if not isinstance(sample, Sample):
        raise TypeError("sample must be a Sample instance.")


    if not isstr(ctype):
        raise TypeError("ctype must be a of type str")

    # validate input
    if ctype != 's' and ctype != 'sum' and \
        ctype != 'r' and ctype != 'rotate' and  \
        ctype != 'i' and ctype != 'incommmensurate' and\
        ctype != 'rnd' and ctype != 'random' :
        raise ValueError("Invalid calculation type.")

    # if 'i', nangles must be defined
    if ctype == 'i' or ctype == 'incommmensurate' or \
        ctype == 'r' or ctype == 'rotate':
        if nangles < 1:
            raise ValueError("Number of angles must be >=1.")
        try:
            nangles = int(nangles)
        except:
            raise ValueError("Cannot convert number of angles to int.")
    if ctype == 'r' or ctype == 'rotate':
        if axis is np.zeros(3):
            raise ValueError("Axis for rotation must be specified.")
        try:
            axis = np.array(axis)
            axis = axis/np.linalg.norm(axis)
        except:
            raise ValueError("Cannot convert axis for rotation to np.ndarray.")

    if ctype == 'rnd' or ctype == 'random':
        if constraints is None:
            raise ValueError("constraints must be specified.")
        if True:
            constr = np.zeros([len(constraints),2], dtype=np.float)
            constr_grp = np.zeros([len(constraints),sample.cell.get_number_of_atoms()], dtype=np.int32)
            for i, c in enumerate(constraints):
                constr_grp[i, c[0]] = 1
                constr[i,0] = c[1][0]
                constr[i,1] = c[1][1]
        #except:
        #    raise ValueError("Cannot convert constraints to np.ndarray.")
    else:
        if constraints is None:
            constraints = 0.
        else:
            constraints = float(constraints)

    try:
        sc = np.array(supercellsize, dtype=np.int32)
    except:
        raise TypeError("Cannot convert supercellsize to NumPy array.")

    if (np.min(sc) <= 0):
        raise ValueError("Supercellsize must be strictly positive.")


    if sc.shape != (3,):
        raise ValueError("Propagation vector has the wrong shape.")

    try:
        r= float(radius) # Lorentz radius (in A)
    except:
        raise TypeError("Cannot convert radius to float.")

    try:
        nnn = int(nnn)
    except:
        raise TypeError("Cannot convert nnn to int.")

    if nnn<0:
        raise ValueError("nnn must be positive.")

    rc=0
    try:
        rc = float(rcont)
    except:
        raise TypeError("Cannot convert rcont to float.")

    if rc<0:
        raise ValueError("rcont must be positive.")

    # check current status is ok
    sample._check_lattice()
    sample._check_magdefs()

    # Prepare input for C extension
    unitcell = sample._cell
    p = unitcell.get_scaled_positions()
    if isinstance(unitcell, Occupations):
        occ  = unitcell.get_occupations()
        occg = unitcell.get_occupations_groups()
        corr = unitcell.get_sites_correlation()
    else:
        occ  = np.ones(len(p))
        occg = np.arange(len(p), dtype=np.int32)
        corr = np.zeros([len(p), len(p)])

    latpar = unitcell.cell.array
    fc = sample.mm.fc   # Fourier components in Bohr Cartesian
    phi = sample.mm.phi # phase in magnetic order definition
    k = sample.mm.k

    if backend == 'available':
        backend = 'clfc2' if lfcext.__version__ == '0.0.3' else 'clfc'

    # Use Numba backend
    if backend == 'nlfc':
        res = []

        if ctype == 's' or ctype == 'sum':
            for mu in sample.muons:
                res.append(LocalFields(*Simple(p, fc,k,phi,mu,sc,latpar,r,nnn,rc, occ, rmin=constraints)))
        else:
            raise NotImplemented("Calculation "+ctype+" not implemented by numba backend")
        return res

    if backend == 'clfc2':
        res = LocalFields(*lfcext.Fields(ctype, p, fc, k,phi,sample.muons,sc,latpar,r,nnn,rc,nangles,axis),
                            nMuons=len(sample.muons))

    else:
        res = []
        # Use CLFC backend
        for mu in sample.muons:
            if ctype == 's' or ctype == 'sum':
                #res.append(LocalFields(*lfcext.Fields(ctype, p,fc,k,phi,mu,sc,latpar,r,nnn,rc)))
                res.append(LocalFields(*lfcext.Simple(p, fc,k,phi,mu,sc,latpar,r,nnn,rc, occ, occg, corr)))
            elif ctype == 'i' or ctype == 'incommensurate':
                res.append(LocalFields(*lfcext.Fields(ctype, p,fc,k,phi,mu,sc,latpar,r,nnn,rc,nangles)))
            elif ctype == 'r' or ctype == 'rotate':
                res.append(LocalFields(*lfcext.Fields(ctype, p,fc,k,phi,mu,sc,latpar,r,nnn,rc,nangles,axis)))

    return res


def dipten(sample, supercellsize, radius):
    """
    Calculates dipolar tensor for given muon sites.

    A fake magnetic order must be defined. This is only used to
    distinguish the ''magnetically active'' atoms for which the
    calculation is perfomred.

    The results are provided in 1/Angstrom^3.
    The conversion factor to mol/emu is 1.6605389 which is given by

    .. math::

        (1 angstrom^{-3} )/(1cm^{-3}) = 1e24

        1/N_A = 1.6605E-24 mole

    :param sample: the sample object
    :param list supercell: the size of the supercell along the lattice coordinates.
    :param float radius: the radius of the sphere used to evaluate the dipolar tensor.
    :return: a list of numpy ndarray containing the dipolar tensor for each muon site defined in the sample.
    :rtype: list
    :raises: TypeError, ValueError: when radius cannot be converted to float or when radius is negative.
    """

    # check current status is ok
    sample._check_lattice()
    sample._check_magdefs()

    try:
        r= float(radius) # Lorentz radius (in A)
        if r<0:
            raise ValueError("Lorentz radius must be greater or equal to 0.")
    except:
        raise TypeError("Cannot convert radius to float.")


    try:
        sc = np.array(supercellsize, dtype=np.int32)
        if sc.shape != (3,):
            raise ValueError("Supercellsize has wrong shape.")
    except:
        raise TypeError("Cannot convert supercellsize to NumPy array.")

    # Remove non magnetic atoms from list

    unitcell = sample.cell

    positions = unitcell.get_scaled_positions()
    latpar = unitcell.cell.array

    ufc = sample.mm.fc

    magnetic_atoms=[]
    for i, e in enumerate(ufc):
        if not np.allclose(e,np.zeros(3,dtype=np.complex)):
            magnetic_atoms.append(i)

    p = positions[magnetic_atoms,:]

    res = []
    for mu in sample.muons:
        res.append(lfcext.DipolarTensor(p,mu,sc,latpar,r))

    return res

