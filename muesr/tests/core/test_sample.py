#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import unittest
import numpy as np
import warnings 

from muesr.core.sample import Sample

from muesr.core.sampleErrors import *
from muesr.core.atoms import Atoms
from muesr.core.magmodel import MM, have_sympy
if have_sympy:
    from muesr.core.magmodel import SMM

from muesr.core.spg import Spacegroup

class TestSample(unittest.TestCase):
 
    def setUp(self):
        self._sample = Sample()
        
    def _set_a_cell(self):
        self._sample.cell = Atoms(symbols=['Co'],
                                  scaled_positions=[[0,0,0]],
                                  cell=[[3.,0,0],
                                        [0,3.,0],
                                        [0,0,3.]])
    def test_set_attribute(self):
        with self.assertRaises(TypeError):
            self._sample.blabla = 1
        
        
    def test_name_property(self):
        
        self.assertEqual(self._sample.name,"No name")
        
        self._sample.name = "test"
        
        self.assertEqual(self._sample.name,"test")
        
        #must be str or unicode
        with self.assertRaises(TypeError):
            self._sample.name = 1
            
        self.assertEqual(self._sample.name,"test")
        
        # test unicode
        self._sample.name = u"àèé"
        
    
    def test_muons_property(self):
        self._sample._reset(cell=True,muon=True,sym=True,magdefs=True)
        self._set_a_cell()
        
        with self.assertRaises(MuonError):
            self._sample.muons
            
        
        
        
    def test_add_muon_property(self):
        self._sample._reset(cell=True,muon=True,sym=True,magdefs=True)

        with self.assertRaises(CellError):
            self._sample.add_muon([0,0,0])
            
        self._set_a_cell()
        
        with self.assertRaises(TypeError):
            self._sample.add_muon('0 0 0')
            
        with self.assertRaises(ValueError):
            self._sample.add_muon([0,0,0,0])
        with self.assertRaises(ValueError):
            self._sample.add_muon(np.array([0,0,0,0]))
        
        self._sample.add_muon([0,0,0])
        np.testing.assert_array_equal(self._sample.muons[0], np.zeros(3))
        
        self._sample.add_muon([1,1,1])
        np.testing.assert_array_equal(self._sample.muons[1], np.ones(3))
        
        self._sample.add_muon([1,1,1], cartesian=True)
        np.testing.assert_array_equal(self._sample.muons[2], np.ones(3)/3.)
        
        a = 4.0  # some lattice constant
        b = a / 2
        self._sample.cell = Atoms(symbols=['Au'],
                                  positions=[0,0,0],
                                  cell=[(0, b, b), (b, 0, b), (b, b, 0)],
                                  pbc=True)
        
        self._sample.add_muon([1.,1.,1.], cartesian=True)
        np.testing.assert_array_equal(self._sample.muons[3], np.ones(3)/4.)
        
        self._sample.add_muon([.5,.5,.5], cartesian=False)
        self._sample.add_muon([2.,2.,2.], cartesian=True)
        np.testing.assert_array_equal(self._sample.muons[4], self._sample.muons[5])
        
        
    
    def test_mm_property(self):
        self._sample._reset(cell=True,magdefs=True,muon=True,sym=True)
        
        with self.assertRaises(MagDefError):
            self._sample.mm

        with self.assertRaises(CellError):
            self._sample.mm = MM(19) # randomly large number
                    
        
        self._sample._reset(magdefs=True)
        self._set_a_cell()
        with self.assertRaises(TypeError):
            self._sample.mm = 1
            
        self._sample._reset(magdefs=True, cell=True, muon=True, sym=True)
        self._sample.cell = Atoms(symbols=['C'],
                                  scaled_positions=[[0,0,0]],
                                  cell=[[3.,0,0],
                                        [0,3.,0],
                                        [0,0,3.]])
                                        
        with self.assertRaises(MagDefError):
            self._sample.mm = MM(198) # randomly large number        
            
        self._sample.mm = MM(1)
        

    def test_new_mm(self):
        self._sample._reset(magdefs=True, cell=True, muon=True, sym=True)
        with self.assertRaises(CellError):
            self._sample.new_mm()
        
        self._set_a_cell()
                                        
        self._sample.new_mm()            
        self.assertTrue(isinstance(self._sample.mm, MM))
        self.assertEqual(len(self._sample.mm.fc), 1)
        
        
    
    def test_new_smm(self):
        if have_sympy:
            self._sample._reset(magdefs=True, cell=True, muon=True, sym=True)
            with self.assertRaises(CellError):
                self._sample.new_smm("x,y,z")
            
            # TODO TESTING
        else:
            pass
            
    def test_current_mm_idx_property(self):
        self._sample._reset(magdefs=True, cell=True,muon=True, sym=True)
        self._sample.cell = Atoms(symbols=['C'],
                                  scaled_positions=[[0,0,0]],
                                  cell=[[3.,0,0],
                                        [0,3.,0],
                                        [0,0,3.]])        
        
        self._sample.new_mm()
        self._sample.mm.k=np.array([0,0,1.])
        self.assertEqual(self._sample.current_mm_idx,0)
        self._sample.new_mm()
        self._sample.mm.k=np.array([0,0,2.])
        self.assertEqual(self._sample.current_mm_idx,1)
        self._sample.new_mm()
        self._sample.mm.k=np.array([0,0,3.])
        self.assertEqual(self._sample.current_mm_idx,2)
        
        self._sample.current_mm_idx = 0
        self.assertEqual(self._sample.current_mm_idx,0)
        np.testing.assert_array_equal(self._sample.mm.k, np.array([0,0,1.]))
        
        with self.assertRaises(IndexError):
            self._sample.current_mm_idx = 3

        with self.assertRaises(IndexError):
            self._sample.current_mm_idx = -1

            
        self.assertEqual(self._sample.current_mm_idx,0)

        
    def test_sym_property(self):
        
        self._sample._reset(cell=True,muon=True,sym=True,magdefs=True)
        
        with self.assertRaises(SymmetryError):
            self._sample.sym
            
        with self.assertRaises(TypeError):
            self._sample.sym = 1
            
        self._sample.sym = Spacegroup(113)
        
        self.assertEqual(self._sample.sym.no,113)
        self.assertEqual(self._sample.sym.setting,1)
        

    def test_cell_property(self): 
        #needs better testing
        with self.assertRaises(CellError):
            self._sample.cell
            
        with self.assertRaises(TypeError):
            self._sample.cell = 1
            
        self._set_a_cell()
        current_cell = self._sample.cell
        self.assertEqual(current_cell.get_chemical_symbols(),['Co'])
        current_cell.set_chemical_symbols(['Co'])
        self.assertEqual(current_cell.get_chemical_symbols(),['Co'])

    def test_reset(self):
        
        # test cell reset
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            self._sample._reset(cell=True)
            assert len(w) == 3
            for wi in w:
                assert issubclass(wi.category, RuntimeWarning)
                
        self.assertEqual(self._sample._cell, None)
        with self.assertRaises(CellError):
            self._sample.cell
        
        # test magdef reset    
        self._sample._reset(magdefs=True)
        self.assertEqual(self._sample._magdefs, [])
        self.assertEqual(self._sample._selected_mm, -1)
        with self.assertRaises(MagDefError):
            self._sample.mm
            
        # test muon reset
        self._set_a_cell()
        self._sample.add_muon([0,1.,2])
        self._sample._reset(muon=True)
        self.assertEqual(self._sample._muon, [])
        with self.assertRaises(MuonError):
            self._sample.muons
        
        self._sample._reset(cell=True,muon=True,sym=True,magdefs=True)
            
        #test symmetry reset
        self._sample.sym = Spacegroup(113)
        self._sample._reset(sym=True)
        self.assertEqual(self._sample._sym, None)
        with self.assertRaises(SymmetryError):
            self._sample.sym
            
                       

    # TODO
    def test_check_sym(self):
        self._sample._reset(cell=True,muon=True,sym=True,magdefs=True)
        
        with self.assertRaises(SymmetryError):
            self._sample._check_sym()
        
        # vary bad from user side
        self._sample._sym = 1
        with self.assertRaises(SymmetryError):
            self._sample._check_sym()
        
        self._sample.sym = Spacegroup(113)
        self.assertTrue(self._sample._check_sym())
        
    def test_check_lattice(self):
        self._sample._reset(cell=True,muon=True,sym=True,magdefs=True)
        
        with self.assertRaises(CellError):
            self._sample._check_lattice()
        
        self._set_a_cell()
        self.assertTrue(self._sample._check_lattice())
        self._sample._cell = 1
        with self.assertRaises(CellError):
            self._sample._check_lattice()
        
            
        
    def test_check_magdefs(self):
        self._sample._reset(cell=True,muon=True,sym=True,magdefs=True)
        
        with self.assertRaises(MagDefError):
            self._sample._check_magdefs()
            
        self._set_a_cell()
        self._sample.new_mm()
        self.assertTrue(self._sample._check_magdefs())
        
        self._sample._magdefs = 1
        with self.assertRaises(MagDefError):
            self._sample._check_magdefs()        
        
        
    def test_check_muon(self):
        self._sample._reset(cell=True,muon=True,sym=True,magdefs=True)
        
        with self.assertRaises(MuonError):
            self._sample._check_muon()

        self._set_a_cell()
        self._sample.add_muon(np.zeros(3))
        self.assertTrue(self._sample._check_muon())
        
        self._sample.add_muon(np.zeros(3))
        self._sample._muon[1] = 'a'
        with self.assertRaises(MuonError):
            self._sample._check_muon()        


if __name__ == '__main__':
    unittest.main()
