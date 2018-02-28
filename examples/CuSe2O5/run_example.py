# coding: utf-8
#
# This example reproduces the results of PhysRevB 87 104413 
#
from muesr.i_o import load_cif
from muesr.engines.clfc import locfield
from muesr.core import Sample

from muesr.utilities import mago_add
from muesr.utilities import show_structure

# Generate a rotation matrix
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    naxis = axis/np.linalg.norm(axis)
    a = np.cos(theta/2.0)
    b, c, d = -naxis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

s = Sample()
load_cif(s,'./cif/CuSe2O5.cif')


show_structure(s)
#mago_add(s,coordinates='b-l')
#0.13 0.5 0
#0.13 -0.5 0
#-0.13 -0.5 0
#-0.13 0.5 0

prbfcs = np.array([[ 0.13+0.j,  0.50+0.j,  0.00+0.j],
       [ 0.13+0.j, -0.50+0.j,  0.00+0.j],
       [-0.13+0.j, -0.50+0.j,  0.00+0.j],
       [-0.13+0.j,  0.50+0.j,  0.00+0.j],
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
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j],
       [ 0.00+0.j,  0.00+0.j,  0.00+0.j]])

mago_add(s,coordinates='b-l', fcs=prbfcs, kvalue=np.array([0.,0.,0.]))

s.add_muon([0.19,0.01,0.23])
s.add_muon([0.33,0.4,0.06])
s.add_muon([0.32,0.44,0.02])
s.add_muon([0.35,0.49,0.32])

results = locfield(s,'s',[50,50,50],100)

# In order to obtain the results in the a*bc coordinate system
# a rotation on -20.7 (110.7-90) degrees must be applied.
for i, f in enumerate(results):
    # rotate in a*bc coordinate system
    tmp=np.dot(f.D,rotation_matrix([0,1,0],-np.pi*(20.7/180)))
    print("Field for site ",i," is ", tmp, " in a*bc coord. system")
