#!/usr/bin/env python
# coding: utf-8

# In[34]:


# In this example a Dipolar field at the muon site of Superconducting Ferromagnet UCoGe is calculated
# space group H-M Pnma  4 Co, U, Ge atoms each at 4c Wyckoff positions
import numpy as np
from muesr.core.sample import Sample                   # Retains all the sample info.
from muesr.i_o.cif.cif import load_cif                 # For loading the structure from cif files
from muesr.utilities import mago_add, show_structure   # For magnetic structure description and visualization
from muesr.utilities import muon_find_equiv            # For finding and including the symmetry equivalent muon positions in the calculation
from muesr.engines.clfc import locfield, find_largest_sphere # Does the sum and returns the local field in its various contributions
np.set_printoptions(suppress=True,precision=5)


# In[35]:


# All the MuSR data of the material can be found in the 10.1103/PhysRevLett.102.167003
# More information of the structure properties can be found in https://doi.org/10.1016/0925-8388(95)02037-3
s=Sample()
load_cif(s,'UCoGe.cif')       # loading the sample crystal
#show_structure(s,[1,1,1])    # Visualise the structure with  XCrysDen.


# In[36]:


s.add_muon([0.0,0.00,0.0])   # muon site added to the sample, experimentally obtain
muon_find_equiv(s)           # find symmetry equivalent site
#show_structure(s,[1,1,1])   # Visualise the structure with  XCrysDen.


# In[37]:


#for i in s.muons:  # to see the equivalent sites
#    print(i)


# In[38]:


#Magnetic structure 

# with magnetic moment of 0.07 Bohr magneton which we assume is located at the U atom || c (10.1103/PhysRevLett.102.167003),
# The structure is ferromagnetic as specify in propagation vector  i.e 0 
# In this case we choose the moment to be || b 
# because we inter change the lattice parameter c with b as used in 10.1103/PhysRevLett.102.167003 
# check the CIF file

mu_u=0.07   # moment in b axis (10.1103/PhysRevLett.102.167003)

s.new_mm()
s.mm.k=np.array([0.0,0.0,0.0])
s.mm.fc=1.*np.array([   [0.0,0.0,0.0],
                        [0.0,0.0,0.0],
                        [0.0,0.0,0.0],
                        [0.0,0.0,0.0],
                        [0.0,mu_u,0.0],
                        [0.0,mu_u,0.0],
                        [0.0,mu_u,0.0],
                        [0.0,mu_u,0.0],
                        [0.0,0.0,0.0],
                        [0.0,0.0,0.0],
                        [0.0,0.0,0.0],
                        [0.0,0.0,0.0]], dtype=np.complex_)
#mago_add(s)    # interactive magnetic structure definition


# In[39]:


#find the largest sphere with center at the muon site(s) 
#that can be  inscribed in a nxnxn supercell to perform the sum .
# Calculate all local contributions to the field contained in r
n=100    
radius=find_largest_sphere(s,[n, n, n])
r=locfield(s, 's', [n, n, n] ,radius)


# In[40]:


cont_coup=-0.00892  #-0.0092
B_dip=np.zeros([len(s.muons),3])
B_Lor=np.zeros([len(s.muons),3])
B_Cont=np.zeros([len(s.muons),3])
B_Tot=np.zeros([len(s.muons),3])
for i in range(len(s.muons)):
    B_dip[i]=r[i].D
    B_Lor[i]=r[i].L
    r[i].ACont= cont_coup
    B_Cont[i]=r[i].C
    B_Tot[i]=r[i].T
    print('net field for site', i+1,':=', np.linalg.norm(B_Tot[i]))
print('')
print('The dipolar field components for all ' +str(len(s.muons))+ ' equivalent sites')
print(B_dip)
print('')
print('Compare the norm of calculate dipolar field = {:.5f} T with the experimental value = 0.015 T\n'.format(np.linalg.norm(B_dip[0])))
print('Lorentz Field, equal for all equivalent site')
print(' {:5.4f} {:5.4f} {:5.4f} T\n'.format(*tuple(B_Lor[0])))

print('contact field component :\n',B_Cont[0])

print('')
print(' results in agreement with 10.1103/PhysRevLett.102.167003\n')
