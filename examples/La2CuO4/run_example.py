import numpy as np
import os
from muesr.engines.clfc import locfield
from muesr.core import Sample
from muesr.engines.clfc import find_largest_sphere
from muesr.i_o import load_cif, save_sample
from muesr.utilities import mago_add, show_structure
    
np.set_printoptions(suppress=True,precision=4)

head,tail = os.path.split(os.path.realpath(__file__))
print("working from "+head)
la2cuo4 = Sample()
load_cif(la2cuo4,os.path.join(head,"La2CuO4_Cmca_new.cif"))
la2cuo4.add_muon([-0.14, 0.1770, -0.1740]) # 1.07 A from apical O


# muon field site 1 from Budnick et al Phys Lett A 124, 103, https://www.sciencedirect.com/science/article/pii/0375960187903823,

la2cuo4.new_mm() # magnetic structure from Vaknin et al PRL 58 2802
# insert the k vector, in a pseudo-ferromagnetic cell, with a base 
la2cuo4.mm.k=np.array([0.0,0.0,0.0])
# insert the m_nu,k fourier components for each atom according to CIF 
# each list of three complex values is an atom nu, 28 atoms
# first 8 La  [0.0+0.j, 0.0+0.j, 0.0+0.j] are non magnetic
# then 4 Cu  [0.0+0.j, -+0.6+0.j, 0.0+0.j] (the red ones below) have 0.6 Bohr magnetons (2D quantum spin 1/2)
# then 16 non magnetic O
# see http://muesr.readthedocs.io/en/latest/Intro.html#description-of-magnetic-structures
la2cuo4.mm.fc= np.array([  [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j],
                                      [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j],
                                      [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j],
                                      [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j],
                                      [0.0+0.j, 0.0+0.j,  0.6+0.j], [0.0+0.j, 0.0+0.j,   0.6+0.j], # magnetic Cu
                                      [0.0+0.j, 0.0+0.j, -0.6+0.j], [0.0+0.j, 0.0+0.j, -0.6+0.j], # magnetic Cu
                                      [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j],
                                      [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j],
                                      [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j],
                                      [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j],
                                      [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j],
                                      [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j],
                                      [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j],
                                      [0.0+0.j, 0.0+0.j,  0.0+0.j], [0.0+0.j, 0.0+0.j,   0.0+0.j] ])
la2cuo4.mm.desc = 'stripe along b with k = 0'


radius=find_largest_sphere(la2cuo4,[100, 100, 100])
print radius

print "Expected answer: [ 0.0035 -0.039  -0.0131] B_dip = 0.04127550110183053 T"


show_structure(la2cuo4,visualizationTool='V')
la2cuo4.current_mm_idx=0
r=locfield(la2cuo4, 's', [100, 100, 100] ,radius)
for v in r:
   print "Calculated ",v.D, 'B_dip =',np.linalg.norm(v.D,axis=0)
   # save_sample(la2cuo4,path+"La2CuO4-stripe-as-fm.yaml")
   #------------------------------------------------------------------------------------------------------------  
   #Experimental internal field at the muon B_dip+B_contact in range of 400-430 G.
   # or 6MHz i.e 410 Gauss from https://www.sciencedirect.com/science/article/pii/0375960187903823 
   #https://link.springer.com/article/10.1007/BF02396016 and   https://link.springer.com/article/10.1007/BF02060647
   
   
