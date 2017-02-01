# coding: utf-8
import numpy as np
from muesr.core.sample import Sample
from muesr.i_o import load_cif
from muesr.utilities import print_cell, muon_find_equiv
from muesr.engines.clfc import locfield
from matplotlib import pyplot as plt

"""

Crystal structure: space group P213, No. 198; Mn-ion at the position 
(0.138,0.138,0.138), and Si-ion at the position (0.845,0.845,0.845); 
lattice constant 4.558 Å.

The magnetic structure of MnSi is characterized by spins forming a 
left-handed incommensurate helix with a propagation vector k≃0.036 Å^−1
in the [111] direction [5–7]. The static Mn moments ( ∼0.4μB for T→0 K)
point in a plane perpendicular to the propagation vector. 

"""


s = Sample()
load_cif(s, 'MnSi.cif')

print_cell(s)

scaled_pos = s.cell.get_scaled_positions()
pos_Mn1 = scaled_pos[0]
pos_Mn2 = scaled_pos[1]
pos_Mn3 = scaled_pos[2]
pos_Mn4 = scaled_pos[3]

# let's do some simple math

a = 4.558 # Ang
a_star = 2*np.pi/a #Ang^-1

norm_k = 0.036 # Å −1

k_astar = k_bstar = k_cstar = (1/np.sqrt(3))*norm_k

k_rlu = np.array([k_astar,k_bstar,k_cstar])/a_star  # this only works for a cubic lattic


# rotation plane must be perpendicular to k // [111]
#  let's chose one direction to be in the [1,-1,0] direction and find the other

def uv(vec):
    return vec/np.linalg.norm(vec)

k_u = uv(k_rlu)     # unitary vector in k direction
a = uv([1,-1,0])    # one of the two directions which define the plane of 
                    #  the rotation

b = np.cross(k_u,a) # the other vector (perp. to a) defining the plane 
                    #  where the moments lie.
                    # There are two choices here: left handed or right 
                    # handed spiral. We will do both.



# assuming that, at x=0, a Mn moment is // a, the spiral can be obtained as
RH_FCs = 0.4 * np.array([(a+1j*b)*np.exp(-2*np.pi*1j*np.dot(k_rlu,pos_Mn1)),
                         (a+1j*b)*np.exp(-2*np.pi*1j*np.dot(k_rlu,pos_Mn2)),
                         (a+1j*b)*np.exp(-2*np.pi*1j*np.dot(k_rlu,pos_Mn3)),
                         (a+1j*b)*np.exp(-2*np.pi*1j*np.dot(k_rlu,pos_Mn4)),
                         [0,0,0],
                         [0,0,0],
                         [0,0,0],
                         [0,0,0]
                        ])
# Notice the minus sign.
LH_FCs = 0.4 * np.array([(a-1j*b)*np.exp(-2*np.pi*1j*np.dot(k_rlu,pos_Mn1)),
                         (a-1j*b)*np.exp(-2*np.pi*1j*np.dot(k_rlu,pos_Mn2)),
                         (a-1j*b)*np.exp(-2*np.pi*1j*np.dot(k_rlu,pos_Mn3)),
                         (a-1j*b)*np.exp(-2*np.pi*1j*np.dot(k_rlu,pos_Mn4)),
                         [0,0,0],
                         [0,0,0],
                         [0,0,0],
                         [0,0,0]
                        ])





s.new_mm()
s.mm.desc = "Right handed spiral"
s.mm.k = k_rlu
s.mm.fc = RH_FCs

s.new_mm()
s.mm.desc = "Left handed spiral"
s.mm.k = k_rlu
s.mm.fc = LH_FCs



s.add_muon([0.532, 0.532, 0.532])
muon_find_equiv(s)

# For Reft Handed
s.current_mm_idx = 0;
r_RH = locfield(s, 'i',[100,100,100],200,nnn=3,nangles=360)

s.current_mm_idx = 1;
r_LH = locfield(s, 'i',[100,100,100],200,nnn=3,nangles=360)


for i in range(4):
    r_RH[i].ACont = -0.066
    r_LH[i].ACont = -0.066


# Two subplots, unpack the axes array immediately
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

ax1.set_title('right-handed')
ax2.set_title('left-handed')


for i in range(4):
    ax1.plot(np.linalg.norm(r_RH[i].T, axis=1))
    ax2.plot(np.linalg.norm(r_LH[i].T, axis=1))

ax1.set_ylabel('Total field [T]')
ax1.set_xlabel('angles [deg]')
ax2.set_xlabel('angles [deg]')

plt.show()
