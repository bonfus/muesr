# coding: utf-8

#   This is an example to calculate the muon local field at the tetrahedral site(s) in bcc-Fe.
import numpy as np
from muesr.core.sample import Sample                   
from muesr.io.cif.cif import load_cif                                
from muesr.io.xsf.xsf import show_supercell, show_cell 
from muesr.utilities.ms import mago_add                
from muesr.engines.clfc import locfield                
from muesr.engines.clfc import find_largest_sphere     
from muesr.utilities.muon import find_equiv            


fe= Sample()                          
load_cif(fe,"./Fe.cif")
fe.add_muon([0.50,0.25,0.0])    
find_equiv(fe)  
fe.new_mm()     
fe.mm.k=np.array([0.0,0.0,0.0])
fe.mm.fc= np.array([[0.0+0.j, 0.0+0.j, 2.22+0.j],
                     [0.0+0.j, 0.0+0.j, 2.22+0.j]])

radius=find_largest_sphere(fe,[100, 100, 100])  
r=locfield(fe, 's', [100, 100, 100] ,radius) 

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
print "Dipolar Field for all the 12 tetrahedral equivalent sites"
print B_dip 
print "The Lorentz field is",B_Lor[0]   
print "The contact field is", np.linalg.norm(B_Cont[0]),"Tesla" 
print "Dipolar average of 1 parallel site and 2 antiparallel sites is",np.linalg.norm(B_dip[3]+B_dip[10]+B_dip[11]),"Tesla "
