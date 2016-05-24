
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    
import unittest
import numpy as np

from muesr.core.sample import Sample
from muesr.core.sampleErrors import *
from muesr.io.sampleIO import *

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

yaml_invalid = """
as;kdj
sdkfjsdlkfj
jdjfgvj
"""


class TestSampleIO(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")
        
    def test_load_only_latt(self):
        
        load_sample("",StringIO(yaml_only_lattice))
        
    def test_load_latt_and_one_mag(self):
        load_sample("",StringIO(yaml_lattice_and_one_mag))
        
        
    def test_load_only_one_mag(self):
        with self.assertRaises(CellError):
            load_sample("",StringIO(yaml_only_one_mag))
    
    def test_invalid_yaml(self):
        with self.assertRaises(ValueError):
            load_sample("",StringIO(yaml_invalid))
        
        
    def test_store_and_load_empty_structure(self):
        l = Sample()
        myfile = StringIO()
        save_sample(l,"",myfile)
        myfile.seek(0)
        t = load_sample("",myfile)
        
        
        
        
if __name__ == '__main__':
    unittest.main()
