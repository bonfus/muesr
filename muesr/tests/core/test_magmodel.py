#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy as np

from muesr.core.magmodel import MM, have_sympy

if have_sympy:
    from muesr.core.magmodel import SMM



class TestMM(unittest.TestCase):
 
    def setUp(self):
        self.NUMFCS = 2
        self._latt = np.random.rand(3,3)
        self._mm = MM(self.NUMFCS, self._latt)
        self._mmnolat = MM(self.NUMFCS)
    
    def test_wrong_init(self):
        with self.assertRaises(TypeError):
            mm = MM(None)
        with self.assertRaises(ValueError):
            mm = MM('a')
        with self.assertRaises(ValueError):
            mm = MM(-1)
        
    
    def test_cannot_set_anything(self):
        with self.assertRaises(TypeError):
            self._mm.ciao = 2
    
    def test_size_property(self):
        
        self.assertEqual(self._mm.size,self.NUMFCS)
        
        
    def test_lattice_params_property(self):
        
        np.testing.assert_array_equal(self._mm.lattice_params, self._latt)

    def test_k_property(self):
        
        np.testing.assert_array_equal(self._mm.k, np.zeros(3))
        
        self._mm.k = np.array([0.1,0.2,0.3])
        
        np.testing.assert_array_equal(self._mm.k, np.array([0.1,0.2,0.3]))
        
        with self.assertRaises(TypeError):
            self._mm.k = 1
        with self.assertRaises(ValueError):
            self._mm.k = np.zeros(4)
        with self.assertRaises(ValueError):
            self._mm.k = np.zeros([4,4])
        with self.assertRaises(ValueError):
            self._mm.k = np.array(['a','a','a'])
            

    def test_fc_property(self):
        
        np.testing.assert_array_equal(self._mm.fc, np.zeros([self.NUMFCS,3]))
        
        np.testing.assert_array_equal(self._mm.fcCart, np.zeros([self.NUMFCS,3]))
        
        np.testing.assert_array_equal(self._mm.fcLattBMA, np.zeros([self.NUMFCS,3]))
        
        np.testing.assert_array_equal(self._mm.fcLattBM, np.zeros([self.NUMFCS,3]))
        
        with self.assertRaises(ValueError):
            np.testing.assert_array_equal(self._mmnolat.fcLattBMA, np.zeros([self.NUMFCS,3]))
        with self.assertRaises(ValueError):    
            np.testing.assert_array_equal(self._mmnolat.fcLattBM, np.zeros([self.NUMFCS,3]))
        
        
        # check for conversions between formats
        randomfcs = np.random.rand(self.NUMFCS,3)+1.j*np.random.rand(self.NUMFCS,3)
        
        origfcs = np.copy(randomfcs)
        
        #simple way to generate lists of the form [1, 2, 0, 1] i.e. first element equal to last
        conv_list = np.mod(np.array([0,1,2,0])+np.random.randint(3),3)
        
        # convert from format in conv_list[0] to format in conv_list[1] and so on
        for i in range(3):
            self._mm.fc_set(randomfcs,int(conv_list[i]))
            randomfcs = self._mm.fc_get(int(conv_list[i+1]))
            
        np.testing.assert_array_almost_equal(origfcs, randomfcs)
        
        
        self._mm.fc = randomfcs
        np.testing.assert_array_almost_equal(randomfcs, self._mm.fc)
        np.testing.assert_array_almost_equal(randomfcs, self._mm.fcCart)
        
        self._mm.fcCart = randomfcs
        np.testing.assert_array_almost_equal(randomfcs, self._mm.fc)
        np.testing.assert_array_almost_equal(randomfcs, self._mm.fcCart)
        
        self._mm.fcLattBMA = randomfcs
        np.testing.assert_array_almost_equal(randomfcs, self._mm.fcLattBMA)
        
        self._mm.fcLattBM = randomfcs
        np.testing.assert_array_almost_equal(randomfcs, self._mm.fcLattBM)
        
        # test wrong type
        with self.assertRaises(ValueError):
            self._mm.fcLattBM = np.zeros(3)
        with self.assertRaises(ValueError):
            self._mm.fcLattBM = np.zeros(3,dtype=np.complex)
        with self.assertRaises(ValueError):
            self._mm.fcLattBM = np.zeros(3,dtype=np.complex)
            
        with self.assertRaises(TypeError):
            self._mm.fc_set(np.zeros([2,3],dtype=np.complex),'cart') #must be int
        
        
        with self.assertRaises(ValueError):
            self._mm.fc_set(np.zeros([2,3],dtype=np.complex),4) #must 0,1,2
        
        with self.assertRaises(TypeError):
            self._mm.fc_set('a',1) #must 0,1,2
        
        
        
    def test_phi_property(self):
        np.testing.assert_array_equal(self._mm.phi, np.zeros(self.NUMFCS))
        
        randomphis = np.random.random(self.NUMFCS)
        self._mm.phi = randomphis
        
        np.testing.assert_array_equal(self._mm.phi, randomphis)
        
        randomphis = np.random.random(self.NUMFCS).tolist()
        self._mm.phi = randomphis
        np.testing.assert_array_equal(self._mm.phi, randomphis)
        
        
        with self.assertRaises(TypeError):
            self._mm.phi = 1
        with self.assertRaises(ValueError):
            self._mm.phi = np.zeros([2,3])

        with self.assertRaises(ValueError):
            self._mm.phi = ['a']

        with self.assertRaises(ValueError):
            self._mm.phi = np.array(['a','a'])
            
        
        
    def test_desc_property(self):
        self._mm.desc = "ciao"
        
        #test unicode
        self._mm.desc = u"ciao"
        
        
        self.assertEqual(self._mm.desc,"ciao")
        
        with self.assertRaises(TypeError):
            self._mm.desc = 1
            
    def test_isSymbolic(self):
        
        self.assertEqual(self._mm.isSymbolic,False)


if have_sympy:
    class TestSMM(unittest.TestCase):
    
        def setUp(self):
            self.NUMFCS = 2
            self._latt = np.random.rand(3,3)
            self._smm = SMM(self.NUMFCS, 'x,y,z', self._latt)
            self._smmnolat = SMM(self.NUMFCS, 'x,y,z')
            
        def test_size_property(self):
            
            self.assertEqual(self._smm.size,self.NUMFCS)
            
        def test_lattice_params_property(self):
            
            np.testing.assert_array_equal(self._smm.lattice_params, self._latt)
    
        def test_k_property(self):
            
            np.testing.assert_array_equal(self._smm.k, np.zeros(3))
            
            self._smm.k = np.array([0.1,0.2,0.3])
            
            np.testing.assert_array_equal(self._smm.k, np.array([0.1,0.2,0.3]))
    
        def test_fc_property(self):

            # random params to test the property
            params = np.random.rand(3)+1.j*np.random.rand(3)

            
            with self.assertRaises(RuntimeError):
                self._smm.set_params(params)
                
            # test coordinate system 0
            self._smm.set_symFC('[[x,y,z],[x,y,z]]', 0)
            
            
            self._smm.set_params(params)
            
            np.testing.assert_array_almost_equal(self._smm.fc, np.array([params,params]))
            np.testing.assert_array_almost_equal(self._smm.fcCart, np.array([params,params]))
            
            paramsInCord1 = self._smm.fc_get(1)
            
            # test coordinate system 1
            self._smm.set_symFC('[[x,y,z],[x,y,z]]', 1)
                        
            self._smm.set_params(paramsInCord1[0].tolist())
                        
            np.testing.assert_array_almost_equal(self._smm.fcLattBMA, paramsInCord1)
            
            paramsInCord2 = self._smm.fc_get(2)
            
            
            # test coordinate system 2
            self._smm.set_symFC('[[x,y,z],[x,y,z]]', 2)
                        
            self._smm.set_params(paramsInCord2[0].tolist())              
            
            np.testing.assert_array_almost_equal(self._smm.fcLattBM, paramsInCord2)
            
            # compare with initial
            
            np.testing.assert_array_almost_equal(self._smm.fc_get(0), np.array([params,params]))
            
            with self.assertRaises(TypeError):
                # wrong size for symbolic fourier components
                self._smmnolat.set_symFC([1,2,3], 1)
            with self.assertRaises(ValueError):
                # wrong size for symbolic fourier components
                self._smmnolat.set_symFC('[[x,y,z]]', 1)
            with self.assertRaises(TypeError):
                # wrong coord sustem value
                self._smmnolat.set_symFC('[[x,y,z],[0,0,0]]', '!')
            
            self._smmnolat.set_symFC('[[x,y,z],[0,0,0]]', -3)
            
            # -3 is an invalid coordinate system
            with self.assertRaises(ValueError):
                self._smmnolat.set_params(params)
                #
                #np.testing.assert_array_equal(self._smmnolat.fcLattBMA, np.zeros([self.NUMFCS,3]))
                #
                #np.testing.assert_array_equal(self._smmnolat.fcLattBM, np.zeros([self.NUMFCS,3]))
            
            
        def test_phi_property(self):
            np.testing.assert_array_equal(self._smm.phi, np.zeros(self.NUMFCS))
            
            randomphis = np.random.random(self.NUMFCS)
            self._smm.phi = randomphis
            
            np.testing.assert_array_equal(self._smm.phi, randomphis)
            
        def test_desc_property(self):
            self._smm.desc = "ciao"
            
            # unicode
            self._smm.desc = u"ciao"
            
            self.assertEqual(self._smm.desc,"ciao")
            
            with self.assertRaises(TypeError):
                
                self._smm.desc = 1
                
        def test_isSymbolic(self):
            
            self.assertEqual(self._smm.isSymbolic,True)


if __name__ == '__main__':
    unittest.main()
