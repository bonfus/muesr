# coding: utf-8

"""
This example shows how to estimate the local fields in LiFePO4 reproducing
the results of Ref. PhysRevB.84.054430. The description of the magnetic order
will be done first using the helper function which facilitate the user input
and then programmatically.
"""

#general import
import numpy as np # to calculate the norm
import os          # join path on windows.  

from muesr.core import Sample            # this object contains all the info on our sample
from muesr.i_o import load_cif            # loads lattice and symmetry from CIF file
from muesr.utilities import  mago_add    # this function provides a CLI for inserting the MAGnetic Order
from muesr.engines.clfc import locfield  # this is the function which actually performs the sum and returns
                                         #  the local fields.

# Nice printing
np.set_printoptions(precision=4, suppress=True)


print('This is an interactive example!')
print('(the automatic simulation will run next and you\'ll be able to compare the results)')

#create a new empty sample
smpl = Sample()

# load structure from cif file
load_cif(smpl, os.path.join('.','cifs','4001848.cif'))

print('To set the magnetic order see:')
print('http://pubs.acs.org/doi/pdf/10.1021/cm0300462 or http://journals.aps.org/prb/pdf/10.1103/PhysRevB.73.024410')
print('\nTo check the results see:')
print('http://journals.aps.org/prb/pdf/10.1103/PhysRevB.84.054430 (Table I in page 6) \n')

print("Tip: if you don't want to invest time in reading the above article check the comments of the python file.")

print('Now set the magnetic order:')

# helper function to add a magnetic model, this will prompt for user input.
mago_add(smpl)

# the interactive session of the above command should be like this:
#     Propagation vector (w.r.t. conv. rec. cell): 0 0 0
#     Magnetic moments in bohr magnetons and cartesian coordinates.
#     Which atom? (enter for all)Fe
#     Lattice vectors:
#       a   10.324400000000001    0.000000000000000    0.000000000000000
#       b    0.000000000000000    6.006400000000000    0.000000000000000
#       c    0.000000000000000    0.000000000000000    4.690100000000000
#     Atomic positions (fractional):
#         1 Fe  0.28220100000000  0.25000000000000  0.97474000000000  55.845
#         2 Fe  0.21779900000000  0.75000000000000  0.47474000000000  55.845
#         3 Fe  0.71779900000000  0.75000000000000  0.02526000000000  55.845
#         4 Fe  0.78220100000000  0.25000000000000  0.52526000000000  55.845
#     FC for atom 1 Fe (3 real, [3 imag]): 0 4.19 0
#     FC for atom 2 Fe (3 real, [3 imag]): 0 -4.19 0
#     FC for atom 3 Fe (3 real, [3 imag]): 0 -4.19 0
#     FC for atom 4 Fe (3 real, [3 imag]): 0 4.19 0




# set muon positions in lattice coordinates
smpl.add_muon([0.1225,0.3772,0.8679])
smpl.add_muon([0.0416,0.2500,0.9172])
smpl.add_muon([0.3901,0.2500,0.3599])
smpl.add_muon([0.8146,0.0404,0.8914])


print('Local fields at the muon sites (Tesla/Bohr Magneton):')
# calculate magnetic fields, the output is a list of LocalFields objects
#  (one for each muon site) that contain information on Total, Dipolar, 
#  Lorentz and Fermi Contact fields (for the Contact term please read the
#  documentation carefully!)
#  
#  The parameters of the following function are:
#            sample   type of calculation     supercell    Lorentz radius
r = locfield( smpl,          's',           [100,100,100],       40)

# we report the results in T/mu_B (same as the paper cited above)
print(r[0].T/4.19, "Norm {: 0.4f}".format(np.linalg.norm(r[0].T/4.19)))
print(r[1].T/4.19, "Norm {: 0.4f}".format(np.linalg.norm(r[1].T/4.19)))
print(r[2].T/4.19, "Norm {: 0.4f}".format(np.linalg.norm(r[2].T/4.19)))
print(r[3].T/4.19, "Norm {: 0.4f}".format(np.linalg.norm(r[3].T/4.19)))


# ==== Now Automatically ====

print('Now automatic run...\n')


# create a new magnetic model for the sample
smpl.new_mm()

# The new magnetic order is automatically selected.
#  one can get the current magnetic model with the
#  property "mm". 

# Set a description for the magnetic order
smpl.mm.desc = "Ferromagnetic"

# the k property can provides and set the propagation vector,
#  always defined in reciprocal lattice units (r.l.u.).
smpl.mm.k=np.array([0.,0.,0.])

# Define Fourier components, in CARTESIAN coordinates and bhor magnetons!
#  if this is done in the script, the components for all the atoms must be set.
FCs = np.array([[ 0.00+0.j,  4.19+0.j,  0.00+0.j],
       [ 0.00+0.j, -4.19+0.j,  0.00+0.j],
       [ 0.00+0.j, -4.19+0.j,  0.00+0.j],
       [ 0.00+0.j,  4.19+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j]])

# set the Fourier components, by default in Cartesian coordinates.
smpl.mm.fc_set(FCs)


print('Local fields at the muon sites (Tesla/Bohr Magneton):')
r = locfield(smpl, 's', [100,100,100],40)

print(r[0].T/4.19, "Norm {: 0.4f}".format(np.linalg.norm(r[0].T/4.19)))
print(r[1].T/4.19, "Norm {: 0.4f}".format(np.linalg.norm(r[1].T/4.19)))
print(r[2].T/4.19, "Norm {: 0.4f}".format(np.linalg.norm(r[2].T/4.19)))
print(r[3].T/4.19, "Norm {: 0.4f}".format(np.linalg.norm(r[3].T/4.19)))
