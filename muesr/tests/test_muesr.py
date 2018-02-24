#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Note: this is not really a unit test but rather a comprehensive
#       test for the main functionalities of the package.
#       It can also be useful to understand some features.

import unittest
import os
import numpy as np

#check extension has been installed
have_lfclib = False
try:
    import lfclib as lfcext
    have_lfclib = True
except:
    pass

from muesr.core.sample import Sample
from muesr.core.magmodel import MM, have_sympy
if have_sympy:
    from muesr.core.magmodel import SMM
from muesr.i_o.xsf.xsf import load_xsf
from muesr.i_o.cif.cif import load_mcif

if have_lfclib:
    from muesr.engines.clfc import locfield
from muesr.utilities.muon import muon_reset, muon_set_frac

#from muesr.core.magmodel import MM

class TestMuesr(unittest.TestCase):
 
    def setUp(self):
        self.assertTrue(have_lfclib,"lfclib not found! Sorry, it's a mandatory")
        cdir = os.path.dirname(__file__)
        self._stdir = os.path.join(cdir,'structures')
    
    def oom(self, v):
        # returns order of magnitude
        return int(np.max(np.log10(np.abs(v))))
 
    def test_one_dipole_outside(self):
        m = Sample()
        load_xsf(m, os.path.join(self._stdir,'crys.xsf'))
        mm = MM(1)
        mm.k=np.array([0.,0.,0.])
        mm.fc_set(np.array([[1.+0j,0.+1.j,0.+1.j]]))
        m.mm = mm
        m.add_muon([0.5,0.5,0.5])
        # the diagonal is 8.6603
        r = locfield(m, 's',[1,1,1],8.66/2.)[0]
        self.assertEqual( np.sum(np.abs(r.D)), 0.)

    def test_one_dipole_inside(self):
        m = Sample()
        load_xsf(m, os.path.join(self._stdir,'crys.xsf'))
        mm = MM(1)
        mm.k=np.array([0.,0.,0.])
        mm.fc_set(np.array([[1.+0j,0.+0.j,0.+0.j]]))
        m.mm = mm
        m.add_muon([0.5,0.5,0.5])
        r = locfield(m,'s',[1,1,1],8.67/2.)[0]
        
        # BDip should be: [0 tesla, 11.4226 milliteslas, 11.4226 milliteslas]
        #
        #    (magnetic_constant/4pi)⋅((1/(1angstrom⋅sqrt(2.5^2+2.5^2+2.5^2)))^3)⋅(3⋅(1 bohr_magneton)⋅cos(atan(sqrt(2)))⋅([1,1,1]/sqrt(3))−1bohr_magneton⋅[1,0,0])
        #
        # BLor should be: [ 0.01138414 T, 0 , 0 ]
        #
        #    0.3333333333⋅magnetic_constant (1 bohr_magneton)/((4/3)pi(angstrom 8.67/2.)^3)
        #
        # BCont should be (Assuming ACont = 1 AA^-3 )
        #
        #	 (2/3)⋅magnetic_constant (1 bohr_magneton)⋅(1angstrom^−3) = 7.769376 T
        #
        
        self.assertAlmostEqual( r.D[0], 0.,places=7, msg=None, delta=None)
        self.assertAlmostEqual( r.D[1], r.D[2], places=10, msg=None, delta=None)
        self.assertAlmostEqual( r.D[1], 11.4226e-3, places=7, msg=None, delta=None)
        
        self.assertAlmostEqual( r.L[0], 0.01138414  ,places=7, msg=None, delta=None)
        self.assertAlmostEqual( r.L[1], 0.0  ,places=7, msg=None, delta=None)
        self.assertAlmostEqual( r.L[2], 0.0  ,places=7, msg=None, delta=None)
        
        r.ACont = 1. # Ang^-3
        
        self.assertAlmostEqual( r.C[0], 7.769376  ,places=7, msg=None, delta=None)
        self.assertAlmostEqual( r.C[1], 0.0  ,places=7, msg=None, delta=None)
        self.assertAlmostEqual( r.C[2], 0.0  ,places=7, msg=None, delta=None)
        
    def test_cubic_sym1(self):
        m = Sample()
        load_xsf(m, os.path.join(self._stdir,'crys.xsf'))
        mm = MM(1)
        mm.k=np.array([0.,0.,0.])
        mm.fc_set(np.array([[1.+0j,0.+0.j,0.+0.j]]))
        m.mm = mm
        m.add_muon([0.5,0.5,0.5])
        r = locfield(m, 's',[100,100,100],240)[0]
        # Bdip should be: [0 tesla, 11.4226 milliteslas, 11.4226 milliteslas]
        #
        #    (magnetic_constant/4pi)⋅((1/(1angstrom⋅sqrt(2.5^2+2.5^2+2.5^2)))^3)⋅(3⋅(1 bohr_magneton)⋅cos(atan(sqrt(2)))⋅([1,1,1]/sqrt(3))−1bohr_magneton⋅[1,0,0])
        #
        #
        self.assertAlmostEqual( r.D[0], r.D[1],places=10, msg=None, delta=None)
        self.assertAlmostEqual( r.D[0], r.D[2], places=10, msg=None, delta=None)
        self.assertAlmostEqual( r.D[0], 0., places=7, msg=None, delta=None)    
        
    def test_cubic_sym2(self):
        m = Sample()
        load_xsf(m,os.path.join(self._stdir,'crys.xsf'))
        mm = MM(1)
        mm.k=np.array([0.,0.,0.])
        mm.fc_set(np.array([[1.+0j,0.+1.j,0.+1.j]]))
        m.mm=mm
        m.add_muon([0.5,0.5,0.5])
        r = locfield(m, 's',[100,100,100],240)[0]
        # r should be: [0 tesla, 11.4226 milliteslas, 11.4226 milliteslas]
        #
        #    (magnetic_constant/4pi)⋅((1/(1angstrom⋅sqrt(2.5^2+2.5^2+2.5^2)))^3)⋅(3⋅(1 bohr_magneton)⋅cos(atan(sqrt(2)))⋅([1,1,1]/sqrt(3))−1bohr_magneton⋅[1,0,0])
        #
        #
        self.assertAlmostEqual( r.D[0], r.D[1],places=10, msg=None, delta=None)
        self.assertAlmostEqual( r.D[0], r.D[2], places=10, msg=None, delta=None)
        self.assertAlmostEqual( r.D[0], 0., places=7, msg=None, delta=None)

    def test_rotation(self):
        m = Sample()
        load_xsf(m, os.path.join(self._stdir,'crys3.xsf'))
        mm = MM(1)
        mm.k=np.array([0.,0.,0.])
        mm.fc_set(np.array([[0.+0j,1.+0.j,0.+0.j]]))
        m.mm = mm
        m.add_muon([0.0,0.001,0.0])
        r = locfield(m, 'r',[1,1,1],1.2,nnn=0,nangles=300,axis=[1,0,0])[0]

        # r for angle = 0 should be: [0 T, 1.8548 T, 0 T]
        #
        #    (magnetic_constant/4pi)⋅((1/(1angstrom⋅sqrt(1^2)))^3)⋅(3⋅(1 bohr_magneton)⋅cos(0)⋅([0,1,0]/sqrt(1))−1bohr_magneton⋅[0,1,0])
        #
        #
        # r for angle = 90 should be: [0 T, 0 T, −0.927401 T]
        # 
        #    (magnetic_constant/4pi)⋅((1/(1angstrom⋅sqrt(1^2)))^3)⋅(3⋅(1 bohr_magneton)⋅cos(pi/2)⋅([0,1,0]/sqrt(1))−1bohr_magneton⋅[0,0,1])
        #
        # Opposite moments produce opposite fields. 150 is half of 300
        np.testing.assert_array_almost_equal(r.D[0],-r.D[150],decimal=7)
        np.testing.assert_array_almost_equal(r.D[0],np.array([0., 1.8548018,0.]),decimal=7)
        np.testing.assert_array_almost_equal(r.D[75],np.array([0., 0,-0.9274009]),decimal=6)
        
        rnorms = np.apply_along_axis(np.linalg.norm,1,r.T-r.L)
        self.assertAlmostEqual( np.min(rnorms), 0.9274009 ,places=6, msg=None, delta=None)
        self.assertAlmostEqual( np.max(rnorms), 1.8548018 ,places=6, msg=None, delta=None)

        r = locfield(m, 'r',[5,1,1],1.5,nnn=0,nangles=300,axis=[1,0,0])[0]
        np.testing.assert_array_almost_equal(r.D[0],-r.D[150],decimal=7)
        np.testing.assert_array_almost_equal(r.D[0],np.array([0., 2.18268753,0.]),decimal=7)
        np.testing.assert_array_almost_equal(r.D[75],np.array([0., 0,-1.58317237]),decimal=6)
        
        rnorms = np.apply_along_axis(np.linalg.norm,1,r.T-r.L)
        self.assertAlmostEqual( np.min(rnorms), 1.58317237 ,places=6, msg=None, delta=None)
        self.assertAlmostEqual( np.max(rnorms), 2.18268753 ,places=6, msg=None, delta=None)
        
        # Now test incomm
        mm = MM(1)
        mm.k=np.array([0.,0.,0.])
        mm.fc_set(np.array([[0.+0j,1.+0.j,0.+1.j]]))
        m.mm = mm
        i = locfield(m, 'i',[1,1,1],1.2,nnn=0,nangles=300)[0]
        
        inorms = np.apply_along_axis(np.linalg.norm,1,i.T-i.L)
        
        self.assertAlmostEqual( np.min(inorms), 0.9274009 ,places=6, msg=None, delta=None)
        self.assertAlmostEqual( np.max(inorms), 1.8548018 ,places=6, msg=None, delta=None)
        
        i = locfield(m, 'i',[5,1,1],1.5,nnn=0,nangles=300)[0]
        
        rnorms = np.apply_along_axis(np.linalg.norm,1,i.T-i.L)
        self.assertAlmostEqual( np.min(rnorms), 1.58317237 ,places=6, msg=None, delta=None)
        self.assertAlmostEqual( np.max(rnorms), 2.18268753 ,places=6, msg=None, delta=None)        
        
        
        
        
    def test_double_moment(self):
        m = Sample()
        load_xsf(m, os.path.join(self._stdir,'crys.xsf'))
        mm = MM(1)
        mm.k=np.array([0.,0.,0.])
        mm.fc_set(np.array([[0.+0j,0.,1.+0.j]]))
        m.mm = mm
        m.add_muon([0.1,0.1,0.1])
        r = locfield(m, 's',[1,1,1],4)[0]
        m.mm.fc_set(np.array([[0.+0j,0.+0.j,2.+0.j]]))
        r2 = locfield(m, 's',[1,1,1],4)[0]
        # r should be: [0 tesla, 11.4226 milliteslas, 11.4226 milliteslas]
        #
        #    (magnetic_constant/4pi)⋅((1/(1angstrom⋅sqrt(2.5^2+2.5^2+2.5^2)))^3)⋅(3⋅(1 bohr_magneton)⋅cos(atan(sqrt(2)))⋅([1,1,1]/sqrt(3))−1bohr_magneton⋅[1,0,0])
        #
        #
        self.assertAlmostEqual( r.D[0], 0.5*r2.D[0], places=10, msg=None, delta=None)
        self.assertAlmostEqual( r.D[1], 0.5*r2.D[1], places=10, msg=None, delta=None)
        self.assertAlmostEqual( r.D[2], 0.5*r2.D[2], places=7, msg=None, delta=None)              

    def test_compare_ass_and_rass(self):
        m = Sample()
        load_xsf(m, os.path.join(self._stdir,'crys2.xsf'))
        mm = MM(9)
        mm.k=np.array(np.random.rand(3))
        mm.fc_set((np.random.rand(27)+1.j*np.random.rand(27)).reshape(9,3))
        m.mm = mm
        m.add_muon(np.random.rand(3))
        m.add_muon(np.random.rand(3))
        r1 = locfield(m, 's',[10,10,10],20)
        r2 = locfield(m, 'r',[10,10,10],20,axis=[1,1,1],nangles=1)

        r1[0].ACont = 1.0
        r2[0].ACont = 1.0
        r1[1].ACont = 1.0
        r2[1].ACont = 1.0
        
        min_oom = min(self.oom(r1[0].T), self.oom(r2[0].T[0]))
        np.testing.assert_array_almost_equal(r1[0].T,r2[0].T[0],decimal=6-min_oom)
        
        min_oom = min(self.oom(r1[1].T), self.oom(r2[1].T[0]))
        np.testing.assert_array_almost_equal(r1[1].T,r2[1].T[0],decimal=6-min_oom)
		
    def test_compare_rass_and_incass(self):
        m = Sample()
        load_xsf(m, os.path.join(self._stdir,'crys.xsf'))
        mm = MM(1)
        # define three orthogonal vectors
        rp1 = np.random.rand(3)   # used for k
        r2 = np.random.rand(3)
        rp2 = np.cross(rp1,r2)
        rp3 = np.cross(rp1,rp2)
        
        # normalize a and b vectors
        rp2 /= np.linalg.norm(rp2)
        rp3 /= np.linalg.norm(rp3)
        
        #chose a random value for stag moment
        stm = np.random.ranf()*10.
        rp2 *= stm
        rp3 *= stm
        
        mm.k=np.array(rp1)
        mm.fc_set((rp2+1.j*rp3).reshape(1,3))
        m.mm = mm
        m.add_muon(np.random.rand(3))
        m.add_muon(np.random.rand(3))
        r1 = locfield(m, 's',[50,50,50],50)
        r2 = locfield(m, 'r',[50,50,50],50,nangles=1200,axis=rp1)
        r3 = locfield(m, 'i',[50,50,50],50,nangles=1200)

        r1[0].ACont = 1.0
        r2[0].ACont = 1.0
        r3[0].ACont = 1.0
        r1[1].ACont = 1.0
        r2[1].ACont = 1.0
        r3[1].ACont = 1.0
        
        min_oom = min(self.oom(r1[0].T), self.oom(r2[0].T[0]))
        np.testing.assert_array_almost_equal(r1[0].T,r2[0].T[0],decimal=7-min_oom)
        min_oom = min(self.oom(r1[1].T), self.oom(r2[1].T[0]))
        np.testing.assert_array_almost_equal(r1[1].T,r2[1].T[0],decimal=7-min_oom)
        
        r2_norms = np.apply_along_axis(np.linalg.norm, 1, r2[0].D)
        r3_norms = np.apply_along_axis(np.linalg.norm, 1, r3[0].D)
        
        min_oom=min(self.oom(np.max(r2_norms)), self.oom(np.max(r3_norms)))
        np.testing.assert_array_almost_equal(np.max(r2_norms),np.max(r3_norms),decimal=4-min_oom)
        
        min_oom=min(self.oom(np.min(r2_norms)), self.oom(np.min(r3_norms)))
        np.testing.assert_array_almost_equal(np.min(r2_norms),np.min(r3_norms),decimal=4-min_oom)

        r2_norms = np.apply_along_axis(np.linalg.norm, 1, r2[1].T)
        r3_norms = np.apply_along_axis(np.linalg.norm, 1, r3[1].T)
        
        min_oom=min(self.oom(np.max(r2_norms)), self.oom(np.max(r3_norms)))
        np.testing.assert_array_almost_equal(np.max(r2_norms),np.max(r2_norms),decimal=5-min_oom)
        
        min_oom=min(self.oom(np.min(r2_norms)), self.oom(np.min(r3_norms)))
        np.testing.assert_array_almost_equal(np.min(r3_norms),np.min(r3_norms),decimal=5-min_oom)


    def test_magmodel_input(self):
        m = Sample()
        load_mcif(m, os.path.join(self._stdir,'LiFeSO4F.mcif'))
        
        fcs_mcif = m.mm.fc
        #print(fcs_mcif)
        
        nmm = MM(m._cell.get_number_of_atoms(),m._cell.get_cell())
        
        fcs = np.zeros([m.cell.get_number_of_atoms(), 3],dtype=np.complex)
        fcs[1:3] = np.array([[0.166987,-0.248559,0.188134],
                        [-0.166987,0.248559,-0.188134]])
        fcs[16:18] = np.array([[0.166987,-0.248559,0.188134],
                        [-0.166987,0.248559,-0.188134]])
        fcs[18:20] = -np.array([[0.166987,-0.248559,0.188134],
                        [0.166987,-0.248559,0.188134]])
        fcs[20:22] = np.array([[0.166987,-0.248559,0.188134],
                        [0.166987,-0.248559,0.188134]])
                        
        
        np.testing.assert_array_almost_equal(m.mm.fcLattBMA,fcs,decimal=4)
        
        m.mm.fc_set(fcs, 1)

        np.testing.assert_array_almost_equal(m.mm.fc,fcs_mcif,decimal=4)
        
        fcs = np.zeros([m.cell.get_number_of_atoms(), 3],dtype=np.complex)
        fcs[1:3] = np.array([[1.73,-2.73,1.36],
                        [-1.73,2.73,-1.36]])
        fcs[16:18] = np.array([[1.73,-2.73,1.36],
                        [-1.73,2.73,-1.36]])
        fcs[18:20] = -np.array([[1.73,-2.73,1.36],
                        [1.73,-2.73,1.36]])
        fcs[20:22] = np.array([[1.73,-2.73,1.36],
                        [1.73,-2.73,1.36]])
        
        np.testing.assert_array_almost_equal(m.mm.fcLattBM,fcs,decimal=4)
        
        m.mm.fc_set(fcs, 2)
        
        np.testing.assert_array_almost_equal(m.mm.fc,fcs_mcif,decimal=4)
        


    def test_symbolic_fourier_components(self):
        m = Sample()
        load_xsf(m, os.path.join(self._stdir,'crys.xsf'))
        
        
        
        rand_k = np.random.rand(3)
        rand_real_fc = np.random.rand(3).reshape((1,3))+0.j # make array complex
        
        mm = MM(1)
        mm.k=np.array(rand_k)
        mm.fc_set(rand_real_fc)
        m.mm = mm
        
        
        m.add_muon([0.1,0.1,0.1])

        r = locfield(m, 's',[10,10,10],40)[0]

        if have_sympy:
            smm = SMM(1,"x,y,z")
            smm.k=np.array(rand_k)
            smm.set_symFC("[[x,y,z]]")
            smm.set_params(rand_real_fc[0].real.tolist())
    
            m.mm = smm
    
            r2 = locfield(m, 's',[10,10,10],40)[0]
    
            self.assertAlmostEqual( r.D[0], r2.D[0], places=10, msg=None, delta=None)
            self.assertAlmostEqual( r.D[1], r2.D[1], places=10, msg=None, delta=None)
            self.assertAlmostEqual( r.D[2], r2.D[2], places=7, msg=None, delta=None)           
        


    def test_mcif(self):
        m = Sample()
        load_mcif(m, os.path.join(self._stdir,'Cd2Os2O7.mcif'))
        m.add_muon([ float(x) for x in '0.125 0.125 0.125'.split() ])
        
        # distance is l(Os1-H) =  2.20122(0) Å
        r1 = locfield(m, 's',[1,1,1],2.3,4,2.3)
        r1[0].ACont = 1.
        
        # total field is 0
        np.testing.assert_array_almost_equal(r1[0].T,np.zeros(3),decimal=7)
        
        
        # distance to nnn is l(Os1-H) =  5.53962(0) Å
        r1 = locfield(m, 's',[3,3,3],5.54,4,2.3)
        r1[0].ACont = 1.
        
        # total field is 0
        np.testing.assert_array_almost_equal(r1[0].T,np.zeros(3),decimal=7)
        
        # now use a funny cell
        r1 = locfield(m, 's',[np.random.randint(3,10),np.random.randint(3,10),np.random.randint(3,10)],5.54,4,2.3)
        r1[0].ACont = 1.
        
        # total field is 0
        np.testing.assert_array_almost_equal(r1[0].T,np.zeros(3),decimal=7)


    def test_mcif2(self):
        m = Sample()
        load_mcif(m, os.path.join(self._stdir,'ScMnO3.mcif'))
        m.add_muon( [float(x) for x in '0.66666666 0.33333333 0.25'.split()])
        
        # No contact
        r1 = locfield(m, 's',[np.random.randint(10,14),np.random.randint(10,14),np.random.randint(10,14)],22,4,0.3)
        r1[0].ACont = 1.
        
        # total field is 0
        np.testing.assert_array_almost_equal(r1[0].T,np.zeros(3),decimal=7)

    def test_mcif3(self):
        m = Sample()
        load_mcif(m, os.path.join(self._stdir,'LiFeSO4F.mcif'))
        muon_set_frac(m, '0.25 0.25 0.25')
        
        # symmetry null
        r1 = locfield(m, 's',[np.random.randint(10,14),np.random.randint(10,14),np.random.randint(10,14)],22,4,10.)
        r1[0].ACont = 1.
        
        # total field is 0
        np.testing.assert_array_almost_equal(r1[0].T,np.zeros(3),decimal=7)      
        
        muon_reset(m)
        m.add_muon([float(x) for x in '0.15990  0.17820  0.14580'.split()])  # position of O3

        # one atom, distance is 2.14251 AA
        #  dipolar field should be:
        #
        #       def bfield(r,m):
        #           return 0.9274009*(3*r*np.dot(r,m)/(np.linalg.norm(r)**5)-m/(np.linalg.norm(r)**3))
        #       bfield(np.array([1.07574,1.59328,0.945826]),np.array([-1.70376,3.14961,-1.22045])) 
        #
        r1 = locfield(m, 's',[np.random.randint(1,4),np.random.randint(1,4),np.random.randint(1,4)],2.2,1,2.14251+20.*np.random.ranf())
        r1[0].ACont = 1.

        np.testing.assert_array_almost_equal(r1[0].D,np.array([0.29531083,  -0.09756855, 0.23347458]),decimal=6) 
        
        # contact field is (2/3)⋅magnetic_constant ([-1.70376 bohr_magneton ,3.14961bohr_magneton,-1.22045bohr_magneton])⋅(1angstrom^−3) =
        #   [-13.2372 T, 24.4705 T, -9.48213 T]
        
        np.testing.assert_array_almost_equal(r1[0].C,np.array([-13.2372 , 24.4705 , -9.48213]),decimal=4) 
        
        
        r2 = locfield(m, 's',[np.random.randint(3,5),np.random.randint(3,5),np.random.randint(3,5)],4.2,3,5.+20.*np.random.ranf())
        r2[0].ACont = 1.14299 # effective interaction increased since more nnn involved in acont
                                # still only 14% more
        
        # three atoms, distances are 2.14251 AA, 4.191093 AA, 4.19413AA
        #
        #      bfield(np.array([1.07574,1.59328,0.945826]),np.array([-1.70376,3.14961,-1.22045])) + \
        #      bfield(np.array([2.14396,2.77751,-2.29775]),np.array([1.70376,-3.14961,1.22045])) + \
        #      bfield(np.array([-2.28034,-2.66189,-2.29775]),np.array([1.70376,-3.14961,1.22045]))        
        
        np.testing.assert_array_almost_equal(r2[0].D,np.array([ 0.20781014, -0.07504131,  0.23329355]),decimal=6)
        
        # two atoms pointing opposite wrt the nearest neighbors...thus just a minus
        #np.testing.assert_array_almost_equal(r2[0].C,np.array([13.2372 , -24.4705 , 9.48213]),decimal=4) 
        
        
        # seven atoms, distances are 2.14251 AA, 4.191093 AA, 4.19413AA, 4.36366 AA 4.50317 AA  4.55780 AA
        #
        #      bfield(np.array([1.07574,1.59328,0.945826]),np.array([-1.70376,3.14961,-1.22045])) + \
        #      bfield(np.array([2.14396,2.77751,-2.29775]),np.array([1.70376,-3.14961,1.22045])) + \
        #      bfield(np.array([-2.28034,-2.66189,-2.29775]),np.array([1.70376,-3.14961,1.22045])) + \
        #      bfield(np.array([0.00751949,0.409054,4.1894]),np.array([1.70376,-3.14961,1.22045])) + \
        #      bfield(np.array([1.83149,-3.84612,0.945826]),np.array([1.70376,-3.14961,1.22045])) + \
        #      bfield(np.array([-4.10431,1.59328,0.945826]),np.array([1.70376,-3.14961,1.22045])) + \
        #      bfield(np.array([2.89971,-2.66189,-2.29775]),np.array([-1.70376,3.14961,-1.22045]))        
        
        #r3 = locfield(m, [np.random.randint(3,5),np.random.randint(3,5),np.random.randint(3,5)],4.6,3,5.+20.*np.random.ranf())
        #r3[0].ACont = 1.
        
        #np.testing.assert_array_almost_equal(r3[0].D,np.array([ 0.24363298, -0.0934996 ,  0.28392374]),decimal=6)

        
        
 
if __name__ == '__main__':
    unittest.main()
