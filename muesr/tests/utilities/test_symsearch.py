try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    
import unittest
import numpy as np
from copy import deepcopy
from muesr.core.sampleErrors import CellError, MuonError
from muesr.core.sample import Sample
from muesr.utilities.symsearch import symsearch, have_spg
from muesr.i_o.cif.cif import read_cif


co_lattice = StringIO("""
#------------------------------------------------------------------------------
#$Date: 2016-02-13 19:28:24 +0000 (Sat, 13 Feb 2016) $
#$Revision: 176429 $
#$URL: file:///home/coder/svn-repositories/cod/cif/1/53/48/1534891.cif $
#------------------------------------------------------------------------------
#
# This file is available in the Crystallography Open Database (COD),
# http://www.crystallography.net/
#
# All data on this site have been placed in the public domain by the
# contributors.
#
data_1534891
loop_
_publ_author_name
'Haglund, J.'
'Fernandez Guillermet, F.'
'Grimvall, G.'
'Korling, M.'
_publ_section_title
;
 Theory of bonding in transition-metal carbides and nitrides
;
_journal_name_full
'Physical Review, Serie 3. B - Condensed Matter (18,1978-)'
_journal_page_first              11685
_journal_page_last               11691
_journal_volume                  48
_journal_year                    1993
_chemical_formula_sum            Co
_chemical_name_systematic        Co
_space_group_IT_number           225
_symmetry_space_group_name_Hall  '-F 4 2 3'
_symmetry_space_group_name_H-M   'F m -3 m'
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                90
_cell_formula_units_Z            4
_cell_length_a                   3.42
_cell_length_b                   3.42
_cell_length_c                   3.42
_cell_volume                     40.002
_citation_journal_id_ASTM        PRBMDO
_cod_data_source_file            Haglund_PRBMDO_1993_1917.cif
_cod_data_source_block           Co1
_cod_original_cell_volume        40.00169
_cod_original_formula_sum        Co1
_cod_database_code               1534891
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
-y,x,z
-x,-y,z
y,-x,z
x,-y,-z
y,x,-z
-x,y,-z
-y,-x,-z
z,x,y
-x,z,y
-z,-x,y
x,-z,y
z,-x,-y
x,z,-y
-z,x,-y
-x,-z,-y
y,z,x
y,-z,-x
z,y,-x
-y,z,-x
-z,-y,-x
-y,-z,x
z,-y,x
-z,y,x
-x,-y,-z
y,-x,-z
x,y,-z
-y,x,-z
-x,y,z
-y,-x,z
x,-y,z
y,x,z
-z,-x,-y
x,-z,-y
z,x,-y
-x,z,-y
-z,x,y
-x,-z,y
z,-x,y
x,z,y
-y,-z,-x
-y,z,x
-z,-y,x
y,-z,x
z,y,x
y,z,-x
-z,y,-x
z,-y,-x
x,y+1/2,z+1/2
-y,x+1/2,z+1/2
-x,-y+1/2,z+1/2
y,-x+1/2,z+1/2
x,-y+1/2,-z+1/2
y,x+1/2,-z+1/2
-x,y+1/2,-z+1/2
-y,-x+1/2,-z+1/2
z,x+1/2,y+1/2
-x,z+1/2,y+1/2
-z,-x+1/2,y+1/2
x,-z+1/2,y+1/2
z,-x+1/2,-y+1/2
x,z+1/2,-y+1/2
-z,x+1/2,-y+1/2
-x,-z+1/2,-y+1/2
y,z+1/2,x+1/2
y,-z+1/2,-x+1/2
z,y+1/2,-x+1/2
-y,z+1/2,-x+1/2
-z,-y+1/2,-x+1/2
-y,-z+1/2,x+1/2
z,-y+1/2,x+1/2
-z,y+1/2,x+1/2
-x,-y+1/2,-z+1/2
y,-x+1/2,-z+1/2
x,y+1/2,-z+1/2
-y,x+1/2,-z+1/2
-x,y+1/2,z+1/2
-y,-x+1/2,z+1/2
x,-y+1/2,z+1/2
y,x+1/2,z+1/2
-z,-x+1/2,-y+1/2
x,-z+1/2,-y+1/2
z,x+1/2,-y+1/2
-x,z+1/2,-y+1/2
-z,x+1/2,y+1/2
-x,-z+1/2,y+1/2
z,-x+1/2,y+1/2
x,z+1/2,y+1/2
-y,-z+1/2,-x+1/2
-y,z+1/2,x+1/2
-z,-y+1/2,x+1/2
y,-z+1/2,x+1/2
z,y+1/2,x+1/2
y,z+1/2,-x+1/2
-z,y+1/2,-x+1/2
z,-y+1/2,-x+1/2
x+1/2,y,z+1/2
-y+1/2,x,z+1/2
-x+1/2,-y,z+1/2
y+1/2,-x,z+1/2
x+1/2,-y,-z+1/2
y+1/2,x,-z+1/2
-x+1/2,y,-z+1/2
-y+1/2,-x,-z+1/2
z+1/2,x,y+1/2
-x+1/2,z,y+1/2
-z+1/2,-x,y+1/2
x+1/2,-z,y+1/2
z+1/2,-x,-y+1/2
x+1/2,z,-y+1/2
-z+1/2,x,-y+1/2
-x+1/2,-z,-y+1/2
y+1/2,z,x+1/2
y+1/2,-z,-x+1/2
z+1/2,y,-x+1/2
-y+1/2,z,-x+1/2
-z+1/2,-y,-x+1/2
-y+1/2,-z,x+1/2
z+1/2,-y,x+1/2
-z+1/2,y,x+1/2
-x+1/2,-y,-z+1/2
y+1/2,-x,-z+1/2
x+1/2,y,-z+1/2
-y+1/2,x,-z+1/2
-x+1/2,y,z+1/2
-y+1/2,-x,z+1/2
x+1/2,-y,z+1/2
y+1/2,x,z+1/2
-z+1/2,-x,-y+1/2
x+1/2,-z,-y+1/2
z+1/2,x,-y+1/2
-x+1/2,z,-y+1/2
-z+1/2,x,y+1/2
-x+1/2,-z,y+1/2
z+1/2,-x,y+1/2
x+1/2,z,y+1/2
-y+1/2,-z,-x+1/2
-y+1/2,z,x+1/2
-z+1/2,-y,x+1/2
y+1/2,-z,x+1/2
z+1/2,y,x+1/2
y+1/2,z,-x+1/2
-z+1/2,y,-x+1/2
z+1/2,-y,-x+1/2
x+1/2,y+1/2,z
-y+1/2,x+1/2,z
-x+1/2,-y+1/2,z
y+1/2,-x+1/2,z
x+1/2,-y+1/2,-z
y+1/2,x+1/2,-z
-x+1/2,y+1/2,-z
-y+1/2,-x+1/2,-z
z+1/2,x+1/2,y
-x+1/2,z+1/2,y
-z+1/2,-x+1/2,y
x+1/2,-z+1/2,y
z+1/2,-x+1/2,-y
x+1/2,z+1/2,-y
-z+1/2,x+1/2,-y
-x+1/2,-z+1/2,-y
y+1/2,z+1/2,x
y+1/2,-z+1/2,-x
z+1/2,y+1/2,-x
-y+1/2,z+1/2,-x
-z+1/2,-y+1/2,-x
-y+1/2,-z+1/2,x
z+1/2,-y+1/2,x
-z+1/2,y+1/2,x
-x+1/2,-y+1/2,-z
y+1/2,-x+1/2,-z
x+1/2,y+1/2,-z
-y+1/2,x+1/2,-z
-x+1/2,y+1/2,z
-y+1/2,-x+1/2,z
x+1/2,-y+1/2,z
y+1/2,x+1/2,z
-z+1/2,-x+1/2,-y
x+1/2,-z+1/2,-y
z+1/2,x+1/2,-y
-x+1/2,z+1/2,-y
-z+1/2,x+1/2,y
-x+1/2,-z+1/2,y
z+1/2,-x+1/2,y
x+1/2,z+1/2,y
-y+1/2,-z+1/2,-x
-y+1/2,z+1/2,x
-z+1/2,-y+1/2,x
y+1/2,-z+1/2,x
z+1/2,y+1/2,x
y+1/2,z+1/2,-x
-z+1/2,y+1/2,-x
z+1/2,-y+1/2,-x
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
Co1 Co 0 0 0 1 0.0
""")

class TestMuon(unittest.TestCase):
 
    def setUp(self):
        
        self._sample = Sample()
        
        self._sample._reset(cell=True,sym=True,magdefs=True,muon=True)
        
        co_lattice.seek(0)
        atoms, sym = read_cif(co_lattice,0) # selectd index 0
        
        if atoms:
            self._sample._reset(muon=True,sym=True)
            self._sample.cell = atoms
            self._sample.sym = sym
        else:
            raise RuntimeError
    
    def test_symsearch(self):
        
        if have_spg:
            # damn slow test!
            origsym = deepcopy(self._sample.sym)
            symsearch(self._sample)
            
            new_rotations, new_translations = self._sample.sym.get_op()
            self.assertEqual(len(new_rotations),len(origsym.get_op()[0]))
            
            for r, t in origsym.get_symop():
                ni = np.argwhere(np.all(new_rotations==r, axis=(1,2)))

                if len(ni) >= 1:
                    found = False
                    for i in ni:
                        # for some reason there is a rounding error so
                        # mod is done on 0.9999999... instead of 1
                        if np.allclose(np.mod(t,1-1e-13),np.mod(new_translations[i[0]],1-1e-13)):
                            found = True

                    if found == False:
                        raise RuntimeError
                        
                else:
                    raise RuntimeError
            
            
        else:
            self.assertFalse(symsearch(self._sample))



if __name__ == '__main__':
    unittest.main()
