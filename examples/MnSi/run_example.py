# coding: utf-8
import numpy as np
from muesr.core.sample import Sample
from muesr.i_o import load_cif
from muesr.utilities import print_cell, muon_find_equiv
from muesr.utilities.muon import muon_reset, muon_find_equiv
from muesr.engines.clfc import locfield, dipten
from matplotlib import pyplot as plt

np.set_printoptions(suppress=True)

"""

Crystal structure: space group P213, No. 198; Mn-ion at the position 
(0.138,0.138,0.138), and Si-ion at the position (0.845,0.845,0.845); 
lattice constant 4.558 Å.

The magnetic structure of MnSi is characterized by spins forming a 
left-handed incommensurate helix with a propagation vector k≃0.036 Å^−1
in the [111] direction [5–7]. The static Mn moments ( ∼0.4μB for T→0 K)
point in a plane perpendicular to the propagation vector. 

"""

print("Create sample...", end='')
s = Sample()
print("done!")

print("Load CIF file", end='')
load_cif(s, 'MnSi.cif')
print("done!")

print("Calculate dipolar tensor for equivalent sites...\n")

# this is a general position along the 111, 
# just to identify the form of the dipolar tensor for the sites 
# along the 111
s.add_muon([0.5,0.5,0.5])

# we find the remainig eq muon sites
muon_find_equiv(s)

# apply an arbitrary small field to select magnetic atoms.
# the abslute value is not used, but it must be different from 0.
APP_FCs = 0.001*np.array([[0,0,1],
                          [0,0,1],
                          [0,0,1],
                          [0,0,1],
                          [0,0,0],
                          [0,0,0],
                          [0,0,0],
                          [0,0,0]], dtype=np.complex)


s.new_mm()
s.mm.desc = "Applied field"
s.mm.k = np.array([0,0,0])
s.mm.fc = APP_FCs

# Calculate the dipolar tensor. Result is in Ang^-3
dts = dipten(s, [30,30,30],50) # supercell size set to 30 unit cells, 
                                  # a 50 Ang sphere si certainly contained.

# Print the muon site for all the 4 equivalent site.
for i, pos in enumerate(s.muons):
    print(("\nFrac. muon position: {:2.3f} {:2.3f} {:2.3f}\n" + \
            "Dipolar Tensor: {:2.3f} {:2.3f} {:2.3f}\n" + \
            "                {:2.3f} {:2.3f} {:2.3f}\n" + \
            "                {:2.3f} {:2.3f} {:2.3f}\n").format( \
              *(pos.tolist() + (dts[i] * 6.022E24/1E24/4).flatten().tolist()) \
            )
          )
      

# We will now identify the position of the muon site using the TF data
# remove all defined positions for the search
muon_reset(s)      

# let's find the muon sites with the dipolar tensor values

# move along 111
positions = np.linspace(0,1,100)
for pos in positions:
    s.add_muon([pos,pos,pos])


dts = dipten(s, [30,30,30],50) # supercell size set to 30 unit cells, 
                                  # a 50 Ang sphere si certainly contained.

values_for_plot = np.zeros_like(positions)                              
for i, dt in enumerate(dts):
    # to reproduce the plot of Amato et. al, we convert to emu/mol
    values_for_plot[i] = dt[0,1] * 6.022E24/1E24/4 # (cm^3/angstrom^3) = 1E24

# remove all defined positions for the search
muon_reset(s)     

fig = plt.figure()
ax = fig.gca()
ax.axhline(y=0, xmin=0., xmax=1., color='b')
ax.axhline(y=-0.2044, xmin=0.4, xmax=0.9, color='r', ls='-',lw=2)
ax.axvline(x=0.532, ymin=-0.0, ymax=1, color='g', ls=':',lw=2)
ax.axvline(x=0.712, ymin=-0.0, ymax=1., color='g', ls=':',lw=2)  # typo in the article?

ax.plot(positions, values_for_plot)
ax.set_ylim([-0.7,0.7])
ax.set_xlabel("Position along the 111 (frac. unit)")
ax.set_ylabel("Dip. Tens. elements [emu/mol]")
plt.savefig("DipolarTensor.png")




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
s.current_mm_idx = 1;      # N.B.: indexes start from 0 but idx=0 is the transverse field!
r_RH = locfield(s, 'i',[50,50,50],100,nnn=3,nangles=360)

s.current_mm_idx = 2;
r_LH = locfield(s, 'i',[50,50,50],100,nnn=3,nangles=360)


for i in range(4):
    r_RH[i].ACont = -0.066
    r_LH[i].ACont = -0.066


# Two subplots, unpack the axes array immediately
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,5))

ax1.set_title('right-handed')
ax2.set_title('left-handed')


for i in range(4):
    ax1.plot(np.linspace(0,360,360),np.linalg.norm(r_RH[i].T, axis=1))
    ax2.plot(np.linspace(0,360,360),np.linalg.norm(r_LH[i].T, axis=1))

ax1.set_ylabel('Total field [T]')
ax1.set_xlabel('angles [deg]')
ax2.set_xlabel('angles [deg]')

plt.savefig("TotalFields.png")


# Repeat with much more angles
s.current_mm_idx = 1;      # N.B.: indexes start from 0 but idx=0 is the transverse field!
r_RH = locfield(s, 'i',[50,50,50],100,nnn=3,nangles=36000)

s.current_mm_idx = 2;
r_LH = locfield(s, 'i',[50,50,50],100,nnn=3,nangles=36000)


for i in range(4):
    r_RH[i].ACont = -0.066
    r_LH[i].ACont = -0.066

# Two subplots, unpack the axes array immediately
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,5))

ax1.set_title('right-handed')
ax2.set_title('left-handed')

N_BINS=1000

LH_Hist=np.zeros(N_BINS)
RH_Hist=np.zeros(N_BINS)
bin_range=np.zeros(N_BINS+1)
for i in range(4):
    hist, bin_range = np.histogram(np.linalg.norm(r_RH[i].T, axis=1), bins=N_BINS, range=(0.08,0.24))
    RH_Hist += hist
    hist, bin_range = np.histogram(np.linalg.norm(r_LH[i].T, axis=1), bins=N_BINS, range=(0.08,0.24))
    LH_Hist += hist

mid_of_bin = bin_range[0:-1]+0.5*np.diff(bin_range)

ax1.plot(mid_of_bin, RH_Hist)
ax2.plot(mid_of_bin, LH_Hist)

ax1.set_ylim([0,1600])
ax2.set_ylim([0,1600])

ax1.set_ylabel('P(B)')
ax1.set_xlabel('B[T]')
ax2.set_xlabel('B[T]')

plt.savefig("Histogram.png")
