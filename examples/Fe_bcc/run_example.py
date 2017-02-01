# coding: utf-8

#   This is an example to calculate the muon local field at the tetrahedral site(s) in bcc-Fe.


import numpy as np
from muesr.core.sample import Sample                   # Retains all the sample info.
from muesr.i_o.cif.cif import load_cif                  # For loading the structure from cif files
from muesr.i_o.xsf.xsf import load_xsf                  # For loading the structure from xsf files
from muesr.i_o.xsf.xsf import show_supercell, show_cell # For visualisation with xcrysden (http://www.xcrysden.org/)
from muesr.utilities.ms import mago_add                # For magnetic structure description
from muesr.engines.clfc import locfield                # Does the sum and returns the local field in its diff. contributions
from muesr.engines.clfc import find_largest_sphere     # Aids in the calculation of the sphere's radius for the lattice sum.
from muesr.utilities.muon import muon_find_equiv       # For finding and including the symmetry equivalent muon positions in the calculation 


#    Declare  and load sample 
fe= Sample()                          
load_cif(fe,"./Fe.cif")

#    To add the muon position
fe.add_muon([0.50,0.25,0.0])    


#   Finds and includes the symmetry equivalent positions of the above defined muon.
#   For this example there are 12 sites "print fe.muons" will show their positions.
muon_find_equiv(fe)  


#   Description of the propagation vector k and fourier components fc with a new magnetic structure declared with fe.new_mm()
fe.new_mm()     
fe.mm.k=np.array([0.0,0.0,0.0])
fe.mm.fc= np.array([[0.0+0.j, 0.0+0.j, 2.22+0.j],
                     [0.0+0.j, 0.0+0.j, 2.22+0.j]])


"""
- This is another way to define the magnetic structure especially when using the interactive session

mago_add(fe)
0.0 0.0 0.0
Fe
0.00 0.00 2.22 
0.00 0.00 2.22 


-It should appear like this on the interactive screen 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mago_add(fe)
	 Propagation vector (w.r.t. conv. rec. cell): 0.0 0.0 0.0
	 Magnetic moments in bohr magnetons and cartesian coordinates.
	 Which atom? (enter for all)Fe
	 Lattice vectors:
	   a    2.868018200000000    0.000000000000000    0.000000000000000
	   b    0.000000000000000    2.868018200000000    0.000000000000000
	   c    0.000000000000000    0.000000000000000    2.868018200000000
	 Atomic positions (fractional):
	     1 Fe  0.00000000000000  0.00000000000000  0.00000000000000  55.845
	     2 Fe  0.50000000000000  0.50000000000000  0.50000000000000  55.845
	 FC for atom 1 Fe (3 real, [3 imag]): 0.00 0.00 2.22 
	 FC for atom 2 Fe (3 real, [3 imag]): 0.00 0.00 2.22 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

- 2.22 is the Fe experimental magnetic moment in Bohr magneton

-Both commands below can be used for visualization on interactive session with xcrysden 
show_cell(fe)                        
show_supercell(fe,[2,2,2]) 

"""



radius=find_largest_sphere(fe,[100, 100, 100])  
r=locfield(fe, 's', [100, 100, 100] ,radius) 

"""
- The find_largest_sphere func. is an experimental function.This func. must not be used, see documentation. 
-  s is for the summation method used. For help type 'help(locfield)' on the interactive python session after the locfield  must have been imported
-  [100, 100, 100] defines the supercell for the calculation
-  radius is the sphere radius
-  All the local fied contributions are contained in r, Bdipolar in r[i].D, BLorentz in r[i].L, BContact in r[i].C, Btotal in r[i].T
-  If r[i].Acont(Contact field term) is not defined r[i].C is all zero.
-  For Acont and its unit, see the documentation (http://muesr.readthedocs.io/en/latest/ContactTerm.html)  

"""


B_dip=np.zeros([len(fe.muons),3])
B_Lor=np.zeros([len(fe.muons),3])
B_Cont=np.zeros([len(fe.muons),3])
B_Tot=np.zeros([len(fe.muons),3])
for i in range(len(fe.muons)):
	B_dip[i]=r[i].D
	B_Lor[i]=r[i].L
	r[i].ACont = 0.0644
	B_Cont[i]=r[i].C
	B_Tot[i]=r[i].T


print("Dipolar Field for all the 12 tetrahedral equivalent sites")
print(B_dip) 

# This is and should be same for all the equivalent sites
print("The Lorentz field is {}".format(B_Lor[0]))   

print("The contact field is {} T".format(np.linalg.norm(B_Cont[0]))) 



"""
Quick look on the description of the muon jumps between the tetrahedral sites as 
discussed in the reference below. "The field contribution at each equivalent site
 is either parallel or antiparallel to the magnetization of the domains" such that
 B_dip(parallel)=-2B_dip(antiparallel) "the average of the dipolar field at these three sites vanishes "

Reference:  M. Schmolz et.al, Hyperfine Interactions 31 (1986) 199-204 

"""

print("Dipolar average of 1 parallel site and 2 antiparallel sites is {} T".format(np.linalg.norm(B_dip[3]+B_dip[10]+B_dip[11])))





"""
- Result on standard output should appear like this:
Using internal LFC library
Dipolar Field for all the 12 tetrahedral equivalent sites
[[  5.92276907e-18   3.82211885e-17   2.64976810e-01]
 [ -2.69848666e-17  -7.10720260e-18   2.64976810e-01]
 [ -1.86130196e-17   4.90398830e-16  -5.29953619e-01]
 [  7.49570991e-16   4.25793961e-16  -5.29953619e-01]
 [  2.85341954e-17   1.46870359e-17   2.64976810e-01]
 [ -9.54449012e-17  -3.56696222e-17   2.64976810e-01]
 [  2.32510025e-17  -9.21149791e-18   2.64976810e-01]
 [ -2.55973467e-17  -3.03263233e-17   2.64976810e-01]
 [  4.74017137e-16   4.22061112e-16  -5.29953619e-01]
 [  4.95076670e-16   4.00493505e-16  -5.29953619e-01]
 [ -2.89631243e-17   4.75531976e-19   2.64976810e-01]
 [ -1.58672778e-17  -1.93742009e-17   2.64976810e-01]]
The Lorentz field is [ 0.          0.          0.73108635]
The contact field is 1.11077214797 T
Dipolar average of 1 parallel site and 2 antiparallel sites is 9.30957573598e-14 T 

"""

