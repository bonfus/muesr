import unittest
import os

from muesr.settings import Settings, config


class TestSettings(unittest.TestCase):
 
    def setUp(self):
        cdir = os.path.dirname(__file__)
        self._stdir = os.path.join(cdir,'structures')
        
    def test_store(self):
        config.store()
        
    def test_class_init(self):
        sss = Settings()
        
    
