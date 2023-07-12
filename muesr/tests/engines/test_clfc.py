#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy as np

from muesr.core.atoms import Atoms
from muesr.core.sample import Sample
from muesr.core.sampleErrors import *
from muesr.core.magmodel import MM, have_sympy
if have_sympy:
    from muesr.core.magmodel import SMM

have_lfclib = False
try:
    import lfclib as lfcext
    have_lfclib = True
except Exception as e:
    print(e)


if have_lfclib:
    from muesr.engines.clfc import LocalFields, find_largest_sphere, locfield

class TestLocalFields(unittest.TestCase):
    def setUp(self):
        self.assertTrue(have_lfclib,"lfclib not found! Sorry, it's a mandatory")
        
    def test_init(self):

        with self.assertRaises(TypeError):
            LocalFields('a',2,None)
            
        with self.assertRaises(ValueError):
            LocalFields(np.zeros(2),np.zeros(3),np.zeros(3))


        with self.assertRaises(ValueError):
            LocalFields(np.array(['a','a']),np.array(['a','a']),np.array(['a','a']))
            
            
        lf = LocalFields(np.zeros(3),np.ones(3),2.*np.ones(3))
        np.testing.assert_array_equal(lf._BLor,2.*np.ones(3))
        np.testing.assert_array_equal(lf._BCont,0.*np.ones(3))
        np.testing.assert_array_equal(lf._BDip,np.ones(3))
        
        
        with self.assertRaises(TypeError):
            lf.bubu = 1
            
    def test_D(self):
        lf = LocalFields(np.ones(3),2.*np.ones(3),3.*np.ones(3))
        np.testing.assert_array_equal(lf.D,2.*np.ones(3))
        
    def test_Dipolar(self):
        lf = LocalFields(np.ones(3),2.*np.ones(3),3.*np.ones(3))
        np.testing.assert_array_equal(lf.Dipolar,2.*np.ones(3))
        
    def test_C(self):
        lf = LocalFields(np.ones(3),2.*np.ones(3),3.*np.ones(3))
        np.testing.assert_array_equal(lf.C,0.*np.ones(3))
        lf.ACont=1.
        np.testing.assert_array_equal(lf.C,np.ones(3))
        
    def test_Contact(self):
        lf = LocalFields(np.ones(3),2.*np.ones(3),3.*np.ones(3))
        np.testing.assert_array_equal(lf.Contact,0.*np.ones(3))
        lf.ACont=1.
        np.testing.assert_array_equal(lf.Contact,1.*np.ones(3))

    def test_L(self):
        lf = LocalFields(np.ones(3),2.*np.ones(3),3.*np.ones(3))
        np.testing.assert_array_equal(lf.L,3.*np.ones(3))
        
    def test_Lorentz(self):
        lf = LocalFields(np.ones(3),2.*np.ones(3),3.*np.ones(3))
        np.testing.assert_array_equal(lf.Lorentz,3.*np.ones(3))
        
    def test_T(self):
        lf = LocalFields(np.ones(3),2.*np.ones(3),3.*np.ones(3))
        np.testing.assert_array_equal(lf.T,5.*np.ones(3))
        lf.ACont = 2
        np.testing.assert_array_equal(lf.T,7.*np.ones(3))
        
    def test_Total(self):
        lf = LocalFields(np.ones(3),2.*np.ones(3),3.*np.ones(3))
        lf.ACont = 3
        np.testing.assert_array_equal(lf.Total,8.*np.ones(3))
        
        
    def test_ACont(self):
        lf = LocalFields(np.ones(3),2.*np.ones(3),3.*np.ones(3))
        lf.ACont = -5.
        self.assertEqual(lf.ACont,-5.)
        np.testing.assert_array_equal(lf.Total,0.*np.ones(3))
        np.testing.assert_array_equal(lf.C,-5.*np.ones(3))

class TestCLFC(unittest.TestCase):
 
    def setUp(self):
        self.assertTrue(have_lfclib)
        self.sample = Sample()

    def _set_a_cell(self):
        self.sample.cell = Atoms(symbols=['Co'],
                                  scaled_positions=[[0,0,0]],
                                  cell=[[3.,0,0],
                                        [0,3.,0],
                                        [0,0,3.]])
        
    def test_find_largest_sphere(self):
        
        self.sample._reset(cell=True,muon=True,sym=True,magdefs=True)
        
        with self.assertRaises(CellError):
            find_largest_sphere(self.sample,[3,3,3])
        
        
        
        self._set_a_cell() 

        with self.assertRaises(MuonError):
            find_largest_sphere(self.sample,[3,3,3])
        
        self.sample.add_muon([0.5,0.5,0.5])
        
        
        with self.assertRaises(TypeError):
            find_largest_sphere(self.sample,1)

        with self.assertRaises(ValueError):
            find_largest_sphere(self.sample,[0,1,2])
            
        self.assertEqual(find_largest_sphere(self.sample,[1,1,1]),1.5)
        self.sample.add_muon([1.,1.,1.],cartesian=True)
        self.assertEqual(find_largest_sphere(self.sample,[1,1,1]),1.)
        
    def test_locfield(self):
        
        self.sample._reset(cell=True,muon=True,sym=True,magdefs=True)
      
        # sample, ctype, supercellsize, radius, nnn = 2, rcont = 10.0, nangles = None, axis = None):
        with self.assertRaises(TypeError):
            locfield(None, None, None, None)
        
        with self.assertRaises(TypeError):
            locfield(self.sample, None, None, None)

        with self.assertRaises(TypeError):
            locfield(self.sample, 's', None, None)
            
        with self.assertRaises(TypeError):
            locfield(self.sample, 's', [2,2,2], None)

        #with self.assertRaises(TypeError):
        #    locfield(self.sample, 's', [2,2,2], 3.)

        with self.assertRaises(TypeError):
            locfield(self.sample, 's', [2,2,2], 3., '2a')

        with self.assertRaises(TypeError):
            locfield(self.sample, 's', [2,2,2], 3., '2', '10b')

        with self.assertRaises(ValueError):
            locfield(self.sample, 'a', [2,2,2], 3.)
            
        with self.assertRaises(ValueError):
            locfield(self.sample, 'i', [2,2,2], 3.)

        with self.assertRaises(ValueError):
            locfield(self.sample, 's', [0,2,2], 3.)
            
        
            
        
        

# http://stackoverflow.com/a/6802723
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    naxis = axis/np.linalg.norm(axis)
    a = np.cos(theta/2.0)
    b, c, d = -naxis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

        
class TestLFCExtension(unittest.TestCase):
    
    def setUp(self):
        self.assertTrue(have_lfclib)
        
    
    def test_one_over_r_cube(self):
        
        p  = np.array([[0.,0.,0.]])
        fc = np.array([[0.,0.,1.]],dtype=np.complex_)
        k  = np.array([0.,0.,0.])
        
        phi= np.array([0.,])
        
        mu = np.array([0.5,0.5,0.5])
        
        sc = np.array([1,1,1],dtype=np.int32)
        latpar = np.diag([2.,2.,2.])
        
        r = 10.
        nnn = 0
        rc=1.
        
        c,d,l = lfcext.Fields('s', p,fc,k,phi,mu,sc,latpar,r,nnn,rc)
        
        # zero hyperfine field since zero nnn
        np.testing.assert_array_equal(c,np.zeros(3))
        
        # this is 0.3333333333⋅magnetic_constant⋅(1 bohr_magneton/(4/3⋅pi⋅(10 angstrom)^3))=9.2740095E-4 tesla
        np.testing.assert_array_almost_equal(l,np.array([0,0,9.2740095E-4]))
        
        mu = np.array([0.5,0.,0.])
        
        c,d,l = lfcext.Fields('s', p,fc,k,phi,mu,sc,latpar,r,nnn,rc)
        
        # (1/(4pi))magnetic_constant⋅(1 bohr_magneton/(1 angstrom^3)) = 0.92740095 tesla
        np.testing.assert_array_almost_equal(d, np.array([0,0,-0.92740095]) )
        
        mu2 = np.array([0.654,0.,0.])

        c,d2,l = lfcext.Fields('s', p,fc,k,phi,mu2,sc,latpar,r,nnn,rc)
        
        # ratios must be like 1/r^3
        np.testing.assert_array_almost_equal(d2, np.array([0,0,-0.92740095])*(1./(np.linalg.norm(mu2*2.))**3) )
        
    
    def test_rotation_of_cart_coord(self):
        
        p  = np.array([[0.1,0.2,0.3]])
        fc = np.array([[0.2,0.4,1.]],dtype=np.complex_)
        k  = np.array([0.2,0.3,0.4])
        
        phi= np.array([0.,])
        
        mu = np.array([0.5,0.5,0.5])
        
        sc = np.array([6,6,6],dtype=np.int32)
        latpar = np.array([[5.7369999886, 0.0, 0.0],
                           [2.2372645948,8.3929280278, 0.0],
                           [1.9062265066,0.8564261924,10.7293797745]])

        
        r = 100000.  #sum everything
        nnn = 4
        rc=6.
        
        
        c,d,l = lfcext.Fields('s', p,fc,k,phi,mu,sc,latpar,r,nnn,rc)
        
        #now rotate lattice system and recalculate
        rmat = rotation_matrix([3.,4.,5.],1.86)
        mrmat = rotation_matrix([3.,4.,5.],-1.86)
        
        rlatpar = np.zeros_like(latpar)
        rlatpar[0,:] = np.dot(rmat,latpar[0,:])
        rlatpar[1,:] = np.dot(rmat,latpar[1,:])
        rlatpar[2,:] = np.dot(rmat,latpar[2,:])
        
        rfc = np.zeros_like(fc)
        rfc[0] = np.dot(rmat, fc[0])
        
        cr,dr,lr = lfcext.Fields('s', p,rfc,k,phi,mu,sc,rlatpar,r,nnn,rc)

        np.testing.assert_array_almost_equal(c,np.dot(mrmat,cr))
        np.testing.assert_array_almost_equal(d,np.dot(mrmat,dr))
        np.testing.assert_array_almost_equal(l,np.dot(mrmat,lr))
        
    
    def test_rotate1(self):
        p  = np.array([[0.,0.,0.]])
        fc = np.array([[0.,0.,1.]],dtype=np.complex_)
        k  = np.array([0.,0.,0.])
        
        phi= np.array([0.,])
        
        mu = np.array([0.5,0.5,0.5])
        
        sc = np.array([1,1,1],dtype=np.int32)
        latpar = np.diag([2.,2.,2.])
        
        r = 10.
        nnn = 0
        rc=1.        
        
        nangles=10
        axis=np.array([0,0,1.])
        
        
        # rotation with axis parallel to local moment
        c,d,l = lfcext.Fields('r', p,fc,k,phi,mu,sc,latpar,r,nnn,rc,nangles,axis)
        
        # this tests that nothing changes since we are rotating with an
        # axis parallel to the moment direction
        np.testing.assert_array_almost_equal(np.diff(d,axis=0), np.zeros([9,3]))
        
        
        mu = np.array([0.5,0.,0.])
        axis=np.array([1.,0.,0.])
        nangles=4

        
        c,d,l = lfcext.Fields('r', p,fc,k,phi,mu,sc,latpar,r,nnn,rc,nangles,axis)
        
        
        np.testing.assert_array_almost_equal(d, np.array([[0,0,-0.92740095],
                                                          [0,0.92740095,0],
                                                          [0,0,0.92740095],
                                                          [0,-0.92740095,0]]))
    
    
    def test_icommensurate(self):
        p  = np.array([[0.,0.,0.]])
        fc = np.array([[0.,1.j,1.]],dtype=np.complex_)
        k  = np.array([0.,0.,0.0])
        
        phi= np.array([0.,])
        
        mu = np.array([0.5,0.,0.])
        
        sc = np.array([1,1,1],dtype=np.int32)
        latpar = np.diag([2.,2.,2.])
        
        r = 10. #all atoms, i.e. 1 atom
        nnn = 1
        rc=10.        
        
        nangles=4
        
        c,d,l = lfcext.Fields('i', p,fc,k,phi,mu,sc,latpar,r,nnn,rc,nangles)
        
        # N.B.: this rotater with the opposite angle with respect to the
        #       rotate function!
        np.testing.assert_array_almost_equal(d, np.array([[0,0,-0.92740095],
                                                          [0,-0.92740095,0],
                                                          [0,0,0.92740095],
                                                          [0,0.92740095,0]]))
        # we have one dipole so Lorentz opposite to dipolar
        np.testing.assert_array_almost_equal(l, np.array([[0,0,0.92740095E-3],
                                                          [0,0.92740095E-3,0],
                                                          [0,0,-0.92740095E-3],
                                                          [0,-0.92740095E-3,0]]))
    
    
    def test_phase(self):
        p  = np.array([[0.,0.,0.]])
        fc = np.array([[0.,0.,1.]],dtype=np.complex_)
        k  = np.array([0.,0.,0.0])
        
        phi= np.array([0.,])
        
        mu = np.array([0.5,0.4,0.3])
        
        sc = np.array([10,10,10],dtype=np.int32)
        latpar = np.diag([2.,2.,2.])
        
        r = 10.
        nnn = 1
        rc=10.
        
        #### simple tests with phase
        refc,refd,refl = lfcext.Fields('s', p,fc,k,phi,mu,sc,latpar,r,nnn,rc)
        refcr,refdr,reflr = lfcext.Fields('r', p,fc,k,phi,mu,sc,latpar,r,nnn,rc,10,np.array([0,1.,0]))
        
        
        phi= np.array([0.25,])
        
        c,d,l = lfcext.Fields('s', p,fc,k,phi,mu,sc,latpar,r,nnn,rc)
        
        np.testing.assert_array_almost_equal(c, np.zeros(3))
        np.testing.assert_array_almost_equal(d, np.zeros(3))
        np.testing.assert_array_almost_equal(l, np.zeros(3))

        c,d,l = lfcext.Fields('r', p,fc,k,phi,mu,sc,latpar,r,nnn,rc,10,np.array([0,1.,0]))
        
        np.testing.assert_array_almost_equal(c, np.zeros([10,3]))
        np.testing.assert_array_almost_equal(d, np.zeros([10,3]))
        np.testing.assert_array_almost_equal(l, np.zeros([10,3]))
        

        phi= np.array([0.5,])
       
        c,d,l = lfcext.Fields('s', p,fc,k,phi,mu,sc,latpar,r,nnn,rc)
        
        np.testing.assert_array_almost_equal(c, -refc)
        np.testing.assert_array_almost_equal(d, -refd)
        np.testing.assert_array_almost_equal(l, -refl)
        
        
        c,d,l = lfcext.Fields('r', p,fc,k,phi,mu,sc,latpar,r,nnn,rc,10,np.array([0,1.,0]))
        np.testing.assert_array_almost_equal(c, -refcr)
        np.testing.assert_array_almost_equal(d, -refdr)
        np.testing.assert_array_almost_equal(l, -reflr)
        
        #### Check spin density dephased
        ##   
        ##   depending on the phase different magnetic structures are 
        ##   obtained. With k=0. only the size is modulated.
        ##   If k = 0.25 one can obtain orders like (for example)
        ##   + 0 - 0     or    + + - - 
        
        p  = np.array([[0.,0.,0.]])
        fc = np.array([[0.,0.,1.j]],dtype=np.complex_)
        k  = np.array([0.,0.,0.0])
        
        phi= np.array([0.125,])
        
        mu = np.array([0.5,0.4,0.])
        
        sc = np.array([10,10,10],dtype=np.int32)
        latpar = np.diag([2.,2.,2.])
        
        r = 10.
        nnn = 1
        rc=10.
        
        refc,refd,refl = lfcext.Fields('s', p,fc,k,phi,mu,sc,latpar,r,nnn,rc)
        
        # no use the equivalent order without phase
        phi= np.array([0.,])
        #pi/4 is sqrt(2)/2
        fc = np.array([[0.,0.,np.sqrt(2.)/2.]],dtype=np.complex_)
        
        c,d,l = lfcext.Fields('s', p,fc,k,phi,mu,sc,latpar,r,nnn,rc)
        
        np.testing.assert_array_almost_equal(c, refc)
        np.testing.assert_array_almost_equal(d, refd)
        np.testing.assert_array_almost_equal(l, refl)
        
        #### now test incommensurate function
        p  = np.array([[0.,0.,0.]])
        fc = np.array([[0.,1.j,1.]],dtype=np.complex_)
        k  = np.array([0.1,0.,0.0])
        
        phi= np.array([0.125,])
        
        mu = np.array([0.5,0.4,0.1])
        
        sc = np.array([20,20,20],dtype=np.int32)
        latpar = np.diag([2.,2.,2.])
        
        r = 20.
        nnn = 1
        rc=5.
        
        refc,refd,refl = lfcext.Fields('i', p,fc,k,phi,mu,sc,latpar,r,nnn,rc,8)
        
        # no use the equivalent order without phase
        phi= np.array([0.,])
        
        c,d,l = lfcext.Fields('i', p,fc,k,phi,mu,sc,latpar,r,nnn,rc,8)
        
        # shifted by one angle (that is the phase)
        np.testing.assert_array_almost_equal(np.take(c,range(1,9),mode='wrap',axis=0), refc)
        np.testing.assert_array_almost_equal(np.take(l,range(1,9),mode='wrap',axis=0), refl)
        np.testing.assert_array_almost_equal(np.take(d,range(1,9),mode='wrap',axis=0), refd)
        

    def test_null_by_symmetry(self):
        p  = np.array([[0.,0.,0.]])
        fc = np.array([[0.,0.,1.]],dtype=np.complex_)
        k  = np.array([0.,0.,0.0])
        
        phi= np.array([0.,])
        
        mu = np.array([0.5,0.5,0.5])
        
        sc = np.array([10,10,10],dtype=np.int32)
        latpar = np.diag([2.,2.,2.])
        
        r = 10.
        nnn = 2
        rc=10.
        
        #### simple tests with phase
        c,d,l = lfcext.Fields('s', p,fc,k,phi,mu,sc,latpar,r,nnn,rc)
        np.testing.assert_array_almost_equal(d, np.zeros(3))
        
        # this works with a large grid but then unit testing is slow!
        #  by the way, 1.9098593 is the ratio between cube and sphere
        #  volumes.
        #np.testing.assert_array_almost_equal(l, np.array([0.,0.,0.92740095/1.9098593]))
        #
        # the comparison done below is with the unconverged results that 
        # is obtained in the 10x10x10 supercell
        np.testing.assert_array_almost_equal(l, np.array([0.,0.,0.46741]),decimal=4)
        
        # (2 magnetic_constant/3)⋅1bohr_magneton   = ((2 ⋅ magnetic_constant) ∕ 3) ⋅ (1 ⋅ bohr_magneton)
        # ≈ 7.769376E-27((g⋅m^3) ∕ (A⋅s^2))
        # ≈ 7.769376 T⋅Å^3
        #
        np.testing.assert_array_almost_equal(c, np.array([0,0,7.769376]))
        
        
        
    def test_dipolar_tensor(self):
        # initial stupid test...
        ###### TODO : do a reasonable test!!!  ######
        p  = np.array([[0.,0.,0.]])
        fc = np.array([[0.,0.,1.]],dtype=np.complex_)
        k  = np.array([0.,0.,0.0])
        
        phi= np.array([0.,])
        
        mu = np.array([0.5,0.5,0.5])
        
        sc = np.array([10,10,10],dtype=np.int32)
        latpar = np.diag([2.,2.,2.])
        
        r = 10.
        res = lfcext.DipolarTensor(p,mu,sc,latpar,r)
        np.testing.assert_array_almost_equal(res, np.zeros([3,3]))
        
        mu = np.array([0.25,0.25,0.25])
        res = lfcext.DipolarTensor(p,mu,sc,latpar,r)
        np.testing.assert_array_almost_equal(np.trace(res), np.zeros([3]))
        np.testing.assert_array_almost_equal(res, res.copy().T)
        
    
if __name__ == '__main__':
    unittest.main()
