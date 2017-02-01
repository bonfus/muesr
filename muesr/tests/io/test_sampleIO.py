
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    
import unittest
import numpy as np

from muesr.core.sample import Sample
from muesr.core.sampleErrors import *
from muesr.i_o.sampleIO import *

yaml_only_lattice = """
Lattice:
  Cell:
  - [2.71, 2.71, 0.0]
  - [2.71, 0.0, 2.71]
  - [0.0, 2.71, 2.71]
  ScaledPositions:
  - [0.0, 0.0, 0.0]
  - [0.25, 0.25, -0.75]
  Symbols: [S, Zn]
MagneticOrders: {}
Muon: {}
Symmetry: {}
"""

yaml_lattice_cartesian_positions = """Lattice:
  Cell:
  - [2.71, 2.71, 0.0]
  - [2.71, 0.0, 2.71]
  - [0.0, 2.71, 2.71]
  CartesianPositions:
  - [0.0, 0.0, 0.0]
  - [1.355, 1.355, 1.355 ]
  Symbols: [S, Zn]
"""

yaml_lattice_and_one_mag = """Lattice:
  Cell:
  - [2.71, 2.71, 0.0]
  - [2.71, 0.0, 2.71]
  - [0.0, 2.71, 2.71]
  ScaledPositions:
  - [0.0, 0.0, 0.0]
  - [0.25, 0.25, -0.75]
  Symbols: [S, Zn]
MagneticOrders:
  Orders:
  - fc:
    - [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    - [7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
    format: b-c
    k: [0.1, 0.2, 0.3]
    lattice:
    - [2.71, 2.71, 0.0]
    - [2.71, 0.0, 2.71]
    - [0.0, 2.71, 2.71]
    phi: [0.1, 0.2]
  Size: 2
Muon: {}
Symmetry: {}
"""

yaml_lattice_and_muon = """Lattice:
  Cell:
  - [2.71, 2.71, 0.0]
  - [2.71, 0.0, 2.71]
  - [0.0, 2.71, 2.71]
  ScaledPositions:
  - [0.0, 0.0, 0.0]
  - [0.25, 0.25, -0.75]
  Symbols: [S, Zn]
MagneticOrders: {}
Muon: 
  Positions:
  - [0.5, 0.5, 0.5]
  - [0.75, 0.25, 0.5]
Symmetry: {}
"""

yaml_only_one_mag = """
Lattice: {}
MagneticOrders:
  Orders:
  - fc:
    - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    format: b-c
    k: [0, 0, 0]
    lattice:
    - [2.71, 2.71, 0.0]
    - [2.71, 0.0, 2.71]
    - [0.0, 2.71, 2.71]
    phi: [0.0, 0.0]
  Size: 2
Muon: {}
Symmetry: {}
"""

yaml_lattice_and_two_mag = """Lattice:
  Cell:
  - [2.71, 2.71, 0.0]
  - [2.71, 0.0, 2.71]
  - [0.0, 2.71, 2.71]
  ScaledPositions:
  - [0.0, 0.0, 0.0]
  - [0.25, 0.25, -0.75]
  Symbols: [S, Zn]
MagneticOrders:
  Orders:
  - fc:
    - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    format: b-c
    k: [0, 0, 0]
    lattice:
    - [2.71, 2.71, 0.0]
    - [2.71, 0.0, 2.71]
    - [0.0, 2.71, 2.71]
    phi: [0.0, 0.0]
  - fc:
    - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    format: b-c
    k: [0, 0, 0]
    lattice:
    - [2.71, 2.71, 0.0]
    - [2.71, 0.0, 2.71]
    - [0.0, 2.71, 2.71]
    phi: [0.0, 0.0]
  Size: 2
Muon: {}
Symmetry: {}
"""

yaml_invalid = """
as;kdj
sdkfjsdlkfj
jdjfgvj
"""

yaml_latt_and_sym = """
Lattice:
  Cell:
  - [2.71, 2.71, 0.0]
  - [2.71, 0.0, 2.71]
  - [0.0, 2.71, 2.71]
  ScaledPositions:
  - [0.0, 0.0, 0.0]
  - [0.25, 0.25, -0.75]
  Symbols: [S, Zn]
MagneticOrders: {}
Muon: {}
Symmetry:
  Number: 216
  Rotations:
  - - [1, 0, 0]
    - [0, 1, 0]
    - [0, 0, 1]
  - - [0, 1, 0]
    - [0, 0, 1]
    - [-1, -1, -1]
  - - [0, 0, 1]
    - [-1, -1, -1]
    - [1, 0, 0]
  - - [-1, -1, -1]
    - [1, 0, 0]
    - [0, 1, 0]
  - - [-1, -1, -1]
    - [0, 0, 1]
    - [0, 1, 0]
  - - [1, 0, 0]
    - [-1, -1, -1]
    - [0, 0, 1]
  - - [0, 1, 0]
    - [1, 0, 0]
    - [-1, -1, -1]
  - - [0, 0, 1]
    - [0, 1, 0]
    - [1, 0, 0]
  - - [-1, -1, -1]
    - [1, 0, 0]
    - [0, 0, 1]
  - - [1, 0, 0]
    - [0, 1, 0]
    - [-1, -1, -1]
  - - [0, 1, 0]
    - [0, 0, 1]
    - [1, 0, 0]
  - - [0, 0, 1]
    - [-1, -1, -1]
    - [0, 1, 0]
  - - [1, 0, 0]
    - [-1, -1, -1]
    - [0, 1, 0]
  - - [0, 1, 0]
    - [1, 0, 0]
    - [0, 0, 1]
  - - [0, 0, 1]
    - [0, 1, 0]
    - [-1, -1, -1]
  - - [-1, -1, -1]
    - [0, 0, 1]
    - [1, 0, 0]
  - - [0, 1, 0]
    - [-1, -1, -1]
    - [0, 0, 1]
  - - [0, 0, 1]
    - [1, 0, 0]
    - [-1, -1, -1]
  - - [-1, -1, -1]
    - [0, 1, 0]
    - [1, 0, 0]
  - - [1, 0, 0]
    - [0, 0, 1]
    - [0, 1, 0]
  - - [0, 0, 1]
    - [1, 0, 0]
    - [0, 1, 0]
  - - [-1, -1, -1]
    - [0, 1, 0]
    - [0, 0, 1]
  - - [1, 0, 0]
    - [0, 0, 1]
    - [-1, -1, -1]
  - - [0, 1, 0]
    - [-1, -1, -1]
    - [1, 0, 0]
  Symbol: F -4 3 m
  Translations:
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
  - [0.0, 0.0, -0.0]
"""

class TestSampleIO(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")
    
    def test_invalid_sample(self):
        with self.assertRaises(ValueError):
            load_sample()
        with self.assertRaises(ValueError):
            load_sample("")
        with self.assertRaises(TypeError):
            save_sample(1,"")
    
    def test_load_only_latt(self):
        s = load_sample("",StringIO(yaml_only_lattice))
        np.testing.assert_array_equal(s.cell.get_cell(),
                                      np.array([[2.71, 2.71, 0.0],
                                                [2.71, 0.0, 2.71],
                                                [0.0, 2.71, 2.71]]))
        np.testing.assert_array_equal(s.cell.get_atomic_numbers(),
                                        np.array([16,30]))
                                        
                                        
        
        
    def test_load_latt_and_one_mag(self):
        s = load_sample("",StringIO(yaml_lattice_and_one_mag))
        np.testing.assert_array_equal(s.cell.get_cell(),
                                      np.array([[2.71, 2.71, 0.0],
                                                [2.71, 0.0, 2.71],
                                                [0.0, 2.71, 2.71]]))
        np.testing.assert_array_equal(s.cell.get_atomic_numbers(),
                                        np.array([16,30]))
                                        
        np.testing.assert_array_equal(s.mm.k, np.array([0.1, 0.2, 0.3]))
        np.testing.assert_array_equal(s.mm.fc, 
                                        np.array([[1.+4j, 2.+5j, 3.+6j],
                                                  [7.+10j, 8.+11j, 9.+12j]]))
        
        np.testing.assert_array_equal(s.mm.phi, np.array([0.1, 0.2]))
        self.assertEqual(s.current_mm_idx,0)
        with self.assertRaises(IndexError):
            s.current_mm_idx = 1

        
    def test_load_only_one_mag(self):
        with self.assertRaises(CellError):
            load_sample("",StringIO(yaml_only_one_mag))
    
    def test_invalid_yaml(self):
        with self.assertRaises(ValueError):
            load_sample("",StringIO(yaml_invalid))
        
    
    def test_load_lattice_and_muon(self):
        s = load_sample("",StringIO(yaml_lattice_and_muon))
        np.testing.assert_array_equal(s.muons[0], np.array([0.5,0.5,0.5]))
        np.testing.assert_array_equal(s.muons[1], np.array([0.75,0.25,0.5]))
        
    def test_load_lattice_and_sym(self):
        s = load_sample("",StringIO(yaml_latt_and_sym))
        np.testing.assert_array_equal(s.cell.get_cell(),
                                      np.array([[2.71, 2.71, 0.0],
                                                [2.71, 0.0, 2.71],
                                                [0.0, 2.71, 2.71]]))
        np.testing.assert_array_equal(s.cell.get_atomic_numbers(),
                                        np.array([16,30]))
        
        np.testing.assert_array_equal(s.sym.translations,np.zeros([24,3]))
        
        self.assertEqual(s.sym.symbol,'F -4 3 m')
        
    
    def test_store_and_load_empty_structure(self):
        l = Sample()
        myfile = StringIO()
        save_sample(l,"",myfile)
        myfile.seek(0)
        t = load_sample("",myfile)
        
    
    def test_load_cartesian_positions(self):
        s = load_sample("",StringIO(yaml_lattice_cartesian_positions))
        spos = s.cell.get_scaled_positions()
        
        np.testing.assert_array_equal(spos,
                                        np.array([[0,0,0],[0.25,0.25,-0.75]])%1)
        
        
if __name__ == '__main__':
    unittest.main()
