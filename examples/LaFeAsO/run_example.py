# coding: utf-8
import glob
import numpy as np
np.set_printoptions(precision=3,suppress=True)

from muesr.core.sample import Sample
from muesr.i_o.cif.cif import load_mcif
from muesr.utilities.symsearch import symsearch
from muesr.utilities import  muon_find_equiv
from muesr.engines.clfc import locfield

#define an empty sample.
m = Sample()

# Find all mcif files.
mcifs=glob.glob("./cifs/*.mcif")

for mcif in mcifs:
# Load all maximally symmetric structures (obtained from bilbao)
#  for propagation vector (1/2,1/2,0)
#  Crystal structure is included and is parsed too.
    load_mcif(m, mcif)
# finds cell symmetries
#  these are inspected from the structure rather than parsed from the cif.
#  The reason for doing this is that structure may be read from files
#  without symmetry informations.
    symsearch(m)
# Sets muon position in fractional coordinates
    m.add_muon([0.375,  0.375, 0.635])
# Find crystalographically equivalent muon sites. 
    muon_find_equiv(m)
    res=locfield(m, 's',[100,100,50],200)
    print('Structure in file: '+mcif)
    for e in res:
        print(e.T*0.66, "Norm {: 0.4f}".format(np.linalg.norm(e.T*0.66)))

print('\n --> According to PhysRevB 80 094524, the local field at the muon site is about 1700 G')
print('\n --> So long range order is either bcs_file_22902c.mcif or ./bcs_file_22902d.mcif \n')
