import unittest

from muesr.core.sampleErrors import CellError
from muesr.core.sample import Sample
from muesr.utilities.printer import print_cell

class TestPrinter(unittest.TestCase):
 
    def setUp(self):
        self._sample = Sample()

    def test_print_cell(self):
        # check that missing sample raises error
        with self.assertRaises(CellError):
            print_cell(self._sample)
        
        
            
if __name__ == '__main__':
    unittest.main()
