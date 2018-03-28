import numpy as np
import os
from muesr.i_o import load_sample
from muesr.engines.clfc import locfield
from muesr.core import Sample
from muesr.engines.clfc import find_largest_sphere
from muesr.i_o import load_cif 
from muesr.utilities.ms import mago_add 
import matplotlib.pyplot as P


np.set_printoptions(suppress=True,precision=4)
head,tail = os.path.split(os.path.realpath(__file__))
print("working from "+head)

def add_muon_points(smp):
     smp.add_muon([0.5,0.0,0.0]) # Agreed experimentally
     #smp.add_muon([0.5,0.0,0.25])
     #smp.add_muon([0.65,0.84,0.0])
     #smp.add_muon([0.69,0.31,0.0])
     #smp.add_muon([0.5 ,0.5 ,0.0])
     

#CoF2 a magnetic insulator\n",
cof = Sample()
load_cif(cof,os.path.join(head,"cif/CoF2.cif"))
add_muon_points(cof)
# compare your structure with CoF2.png from https://doi.org/10.1103/PhysRevB.30.186

# magnetic moment of 2.6 muB from https:doi.//org/10.1103/PhysRevB.87.121108  https://doi.org/10.1103/PhysRevB.69.014417
cof.new_mm()
cof.mm.k=np.array([0.0,0.0,1.0])
# according to CoF2.cif (setting with a,b equal, c shorter, type cif to check)
# H-M P4_2/mnm group 136, six atoms in the cell, in this order
# Co at 0.00000 0.00000 0.00000 (2b site)  
# the symmetry replica is generated at 0.5000 0.5000 0.5000
# http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-normsets?&norgens=&gnum=136 
# F at 0.30600 0.30600 0.00000  (4f site)
# the symmetry replicas are generated at 1--x 1-x 0, 0.5+x 0.5-x, 0.5, 0.5-x,0.5+x, 0.5
cof.mm.fc= np.array([[0.0+0.j, 0.0+0.j, 2.6+0.j],[0.0+0.j, 0.0+0.j, -2.6+0.j], 
                              [0.0+0.j, 0.0+0.j, 0.0+0.j],[0.0+0.j, 0.0+0.j,  0.0+0.j],
                              [0.0+0.j, 0.0+0.j, 0.0+0.j],[0.0+0.j, 0.0+0.j,  0.0+0.j] ])

   
show_structure(cof,visualizationTool='V')  # show_structure(cof,supercell=[1,1,2],visualizationTool='V')
n=18
radius=find_largest_sphere(cof,[n,n,n]) 
B_dip = np.linalg.norm(locfield(cof, 's', [n,n,n] ,radius)[0].D,axis=0)
print('R = {:.1f}, B_dip = {:.4f} T'.format(radius,B_dip))
npoints = 11
n = np.logspace(0.53,2,npoints,dtype=int)
k = -1
B_dip = np.zeros(npoints)
R = np.zeros(npoints)
for m in n:
     k += 1
     radius=find_largest_sphere(cof,[m,m,m])
     r=locfield(cof, 's', [m, m, m] ,radius) #
     R[k] = radius
     B_dip[k] =np.linalg.norm(r[0].D,axis=0)
fig,ax = P.subplots()
ax.plot(R,B_dip,'bo',label='sum')
ax.plot(R,R-R+0.265,'b--',label='exp')
ax1 = ax.twinx()
ax.set_xlabel('R')
ax.set_ylabel(r'$B_d$  [T]')
ax1.plot(R,n,'rd')
ax1.set_ylabel('m (sites per cube edge)')
ax.legend(loc=9)
P.show()
#Experimental results at site 1, Octahedral site is 2650 Gauss from (DOI:https://doi.org/10.1103/PhysRevB.30.186)

