
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    
import unittest
import numpy as np
import sys

from muesr.core.sample import Sample
from muesr.core.sampleErrors import *
from muesr.i_o import load_xsf, save_xsf


class TestXsfIO(unittest.TestCase):
    def setUp(self):
        pass
    
    @unittest.skipIf(sys.version_info[0] >= 3, 'Python3 specific test')
    def test_open_invalid_file(self):
        s = Sample()
        with self.assertRaises(FileNotFoundError):
            load_xsf(s,'ciao')
            
    @unittest.skipIf(sys.version_info[0] < 3, 'Python2 specific test')
    def test_open_invalid_file(self):
        s = Sample()
        with self.assertRaises(OSError):
            load_xsf(s,'ciao')
    
        # check for python2
        with self.assertRaises(OSError):
            load_xsf(s,u'ciao')
    
    def test_save_file(self):
        s = Sample()
        with self.assertRaises(CellError):
            save_xsf(s,u'ciao')
        

        
if __name__ == '__main__':
    unittest.main()
