# Copyright (C) 2011 Atsushi Togo
# All rights reserved.
#
# This file is part of phonopy.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# * Neither the name of the phonopy project nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import numpy as np
from muesr.core.atoms import Atoms
from muesr.core.sampleErrors import MagDefError


def get_simple_supercell(sample,  multi):
    """
    This function creates a simple supercell by expanding the unit cell
    in a, b and c directions.
    
    
    """
    # Scaled positions within the frame, i.e., create a supercell that
    # is made simply to multiply the input cell.
    
    unitcell = sample.cell
    
    if type(multi) is np.ndarray:
        multi = multi.tolist()
    
    if type(multi) is list:
        if len(multi) == 3:
            try:
                multi = [int(x) for x in multi]
            except:
                raise TypeError('Cannot convert multi to int')
            if multi[0] <= 0 or multi[1] <= 0 or multi[2] <= 0:
                raise ValueError('Supercell values must be strictly positive.')
            else:
                # everything is fine!
                pass
        else:
            raise ValueError('multi must be a 3D vector!')
    else:
        raise TypeError('multi must me list or numpy array of integers' +
                        ' (automatically converted)')
    
    have_mag_structure = True
    
    FC=None
    K=None
    PHI=None
    
    try:
        FC = sample.mm.fc
        K = sample.mm.k 
        PHI  = sample.mm.phi
    except MagDefError:
        have_mag_structure = False
        
    
    positions = unitcell.get_scaled_positions()
    numbers = unitcell.get_atomic_numbers()
    masses = unitcell.get_masses()
    lattice = unitcell.get_cell()
    
    # these lists are used to strore new values
    positions_multi = []
    numbers_multi = []
    masses_multi = []
    
    if not have_mag_structure:
        magmoms_multi = None
    else:
        magmoms_multi = []

    for l, pos in enumerate(positions):
        if numbers[l] == 0:    #  Check again if muon in there!
            raise RuntimeError #  This shuld never happen!
            continue
        for k in range(multi[2]):
            for j in range(multi[1]):
                for i in range(multi[0]):
                    positions_multi.append([ (pos[0] + i) / multi[0],
                                             (pos[1] + j) / multi[1],
                                             (pos[2] + k) / multi[2] ])
                    numbers_multi.append(numbers[l])
                    masses_multi.append(masses[l])
                    
                    
                    if have_mag_structure:

                        
                        c = np.cos ( 2.0*np.pi * (np.dot(K,[float(i),float(j),float(k)]) + PHI[l]))
                        s = np.sin ( 2.0*np.pi * (np.dot(K,[float(i),float(j),float(k)]) + PHI[l]));
                        
                        sk = np.real(FC[l])
                        isk = np.imag(FC[l])
                        
                        m = np.zeros(3)
                        m = c*sk + s*isk
                        
                        #print "Norma: " , np.linalg.norm(m)
                        magmoms_multi.append(m)                               



    return Atoms(numbers = numbers_multi,
                 masses = masses_multi,
                 magmoms = magmoms_multi,
                 scaled_positions = positions_multi,
                 cell = np.dot( np.diag( multi ), lattice ),
                 pbc=True)




def print_cell(cell, mapping=None):
    """
    Print lattice structure.
    
    cell : Aroms instance, the atomic structure to be printed
    mapping : list or None, can be used to show relations between atoms.
    
    """
    symbols = cell.get_chemical_symbols()
    masses = cell.get_masses()
    magmoms = cell.get_magnetic_moments()
    lattice = cell.get_cell()
    print ("Lattice vectors:")
    print ("  a %20.15f %20.15f %20.15f" % tuple( lattice[0] ))
    print ("  b %20.15f %20.15f %20.15f" % tuple( lattice[1] ))
    print ("  c %20.15f %20.15f %20.15f" % tuple( lattice[2] ))
    print ("Atomic positions (fractional):")
    for i, v in enumerate(cell.get_scaled_positions()):
        if magmoms == None:
            print ("%5d %-2s%18.14f%18.14f%18.14f %7.3f" % \
                (i+1, symbols[i], v[0], v[1], v[2], masses[i]))
        else:
            print ("%5d %-2s%18.14f%18.14f%18.14f %7.3f  %s" % \
                (i+1, symbols[i], v[0], v[1], v[2], masses[i], ', '.join('%.3f' % val for val in magmoms[i])))
        # print 
        if mapping == None:
            print()
        else:
            print (">", mapping[i]+1)



#
# Get distance between a pair of atoms
#
def get_distance(atoms, a0, a1, tolerance=1e-5):
    """
    Return the shortest distance between a pair of atoms in PBC
    """
    reduced_bases = get_reduced_bases( atoms.get_cell(), tolerance)
    scaled_pos = np.dot( atoms.get_positions(), np.linalg.inv(reduced_bases) )
    # move scaled atomic positions into -0.5 < r <= 0.5
    for pos in scaled_pos:
        pos -= pos.round()

    # Look for the shortest one in surrounded 3x3x3 cells 
    distances = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                distances.append( np.linalg.norm(
                        np.dot(scaled_pos[a0] - scaled_pos[a1] + np.array([i,j,k]),
                               reduced_bases) ) )
    return min(distances)

def get_distance_with_center( atoms, center, atom_num, tolerance=1e-5 ):
    """
    Return the shortest distance to atom from specified position
    """
    reduced_bases = get_reduced_bases( atoms.get_cell(), tolerance)
    scaled_pos = np.dot( atoms.get_positions(), np.linalg.inv(reduced_bases) )
    # move scaled atomic positions into -0.5 < r <= 0.5
    #for pos in scaled_pos:
    #    pos -= pos.round()
    #print pos
    #print scaled_pos
    # Look for the shortest one in surrounded 3x3x3 cells 
    distances = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                distances.append( np.linalg.norm(
                        np.dot(scaled_pos[atom_num] - center + np.array([i,j,k]),
                               reduced_bases) ) )
    return min(distances)

def get_atom_distance( scaled_pos, center, cell, tolerance=1e-5 ):
    """
    Return the shortest distance to atom from specified position
    """
    reduced_bases = get_reduced_bases(cell, tolerance)
    # move scaled atomic positions into -0.5 < r <= 0.5
    #for pos in scaled_pos:
    #    pos -= pos.round()
    #print pos
    #print scaled_pos
    # Look for the shortest one in surrounded 3x3x3 cells 
    distances = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                distances.append( np.linalg.norm(
                        np.dot(scaled_pos - center + np.array([i,j,k]),
                               reduced_bases) ) )
    return min(distances)    

def get_atom_vec( scaled_pos, center, cell, tolerance=1e-5 ):
    """
    Return the shortest distance to atom from specified position
    """
    reduced_bases = get_reduced_bases(cell, tolerance)
    # move scaled atomic positions into -0.5 < r <= 0.5
    #for pos in scaled_pos:
    #    pos -= pos.round()
    #print pos
    #print scaled_pos
    # Look for the shortest one in surrounded 3x3x3 cells 
    distances = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                a = np.dot(scaled_pos - center + np.array([i,j,k]),
                               reduced_bases)
                distances.append(
                        np.dot(scaled_pos - center + np.array([i,j,k]),
                               reduced_bases) )
                print (np.linalg.norm(a))
                print ("i j k: %i %i %i" % (i,j,k))
    #print distances
    #print [np.linalg.norm(item) for item in distances]
    return min(distances, key = lambda item: (np.linalg.norm(item)))     


#
# Delaunay reduction
#    
def get_reduced_bases(cell, tolerance=1e-5):
    """
    This is an implementation of Delaunay reduction.
    Some information is found in International table.
    """
    return get_Delaunay_reduction(cell, tolerance)

def get_Delaunay_reduction(lattice, tolerance):
    extended_bases = np.zeros((4,3), dtype=float)
    extended_bases[:3,:] = lattice
    extended_bases[3] = -np.sum(lattice, axis=0)

    for i in range(100):
        if reduce_bases(extended_bases, tolerance):
            break
    if i == 99:
        print ("Delaunary reduction is failed.")

    shortest = get_shortest_bases_from_extented_bases(extended_bases, tolerance)

    return shortest

def reduce_bases(extended_bases, tolerance):
    metric = np.dot(extended_bases, extended_bases.transpose())
    for i in range(4):
        for j in range(i+1, 4):
            if metric[i][j] > tolerance:
                for k in range(4):
                    if (not k == i) and (not k == j):
                        extended_bases[k] += extended_bases[i]
                extended_bases[i] = -extended_bases[i]
                extended_bases[j] = extended_bases[j]
                return False

    # Reduction is completed.
    # All non diagonal elements of metric tensor is negative.
    return True

def get_shortest_bases_from_extented_bases(extended_bases, tolerance):

    def mycmp(x, y):
        return cmp(np.vdot(x,x), np.vdot(y,y))

    basis = np.zeros((7,3), dtype=float)
    basis[:4] = extended_bases
    basis[4]  = extended_bases[0] + extended_bases[1]
    basis[5]  = extended_bases[1] + extended_bases[2]
    basis[6]  = extended_bases[2] + extended_bases[0]
    # Sort bases by the lengthes (shorter is earlier)
    basis = sorted(basis, cmp=mycmp)
    
    # Choose shortest and linearly independent three bases
    # This algorithm may not be perfect.
    for i in range(7):
        for j in range(i+1, 7):
            for k in range(j+1, 7):
                if abs(np.linalg.det([basis[i],basis[j],basis[k]])) > tolerance:
                    return np.array([basis[i],basis[j],basis[k]])

    print ("Delaunary reduction is failed.")
    return basis[:3]

#
# Other tiny tools
#    
def get_angles( lattice ):
    """
    Get alpha, beta and gamma angles from lattice vectors.
    
    >>> get_angles( np.diag([1,2,3]) )
    (90.0, 90.0, 90.0)
    """
    a, b, c = get_cell_parameters( lattice )
    alpha = np.arccos(np.vdot(lattice[1], lattice[2]) / b / c ) / np.pi * 180
    beta  = np.arccos(np.vdot(lattice[2], lattice[0]) / c / a ) / np.pi * 180
    gamma = np.arccos(np.vdot(lattice[0], lattice[1]) / a / b ) / np.pi * 180
    return alpha, beta, gamma

def get_cell_parameters( lattice ):
    return np.sqrt( np.dot ( lattice, lattice.transpose() ).diagonal() )

def get_cell_matrix( a, b, c, alpha, beta, gamma ):
    # These follow 'matrix_lattice_init' in matrix.c of GDIS
    alpha *= np.pi / 180
    beta *= np.pi / 180
    gamma *= np.pi / 180
    a1 = a
    a2 = 0.0
    a3 = 0.0
    b1 = np.cos( gamma )
    b2 = np.sin( gamma )
    b3 = 0.0
    c1 = np.cos( beta )
    c2 = ( 2 * np.cos( alpha ) + b1**2 + b2**2 - 2 * b1 * c1 - 1 ) / ( 2 * b2 )
    c3 = np.sqrt( 1 - c1**2 - c2**2 )
    lattice = np.zeros( ( 3, 3 ), dtype=float )
    lattice[ 0, 0 ] = a
    lattice[ 1 ] = np.array( [ b1, b2, b3 ] ) * b
    lattice[ 2 ] = np.array( [ c1, c2, c3 ] ) * c
    return lattice

def get_reciprocal_lattice( lattice ):
    #return reciprocal lattice vectors
    
    #Construct the metric tensor
    g = gtensor(lattice)

    # inverse of metric tensor
    gi = 2 * np.pi * np.linalg.inv(g)
    r_lattice = np.dot(gi,lattice)
    #V = np.sqrt(np.linalg.det(g))
    #print V
    return r_lattice


def gtensor(lattice):
    "calculates the metric tensor of a lattice"
    g=np.zeros((3, 3), 'd')
    #print 'shape ', g.shape
    a, b, c = get_cell_parameters(lattice)
    alpha, beta, gamma = get_angles(lattice)

    alpha *= np.pi / 180
    beta *= np.pi / 180
    gamma *= np.pi / 180
    
    g[0, 0]=a**2;
    g[0, 1]=a*b*np.cos(gamma)
    g[0, 2]=a*c*np.cos(beta)

    g[1, 0]=g[0, 1]
    g[1, 1]=b**2
    g[1, 2]=c*b*np.cos(alpha)

    g[2, 0]=g[0, 2]
    g[2, 1]=g[1, 2]
    g[2, 2]=c**2
    return g

    
def reciprocate(x, y, z, cell):
    "calculate miller indexes of a vector defined by its fractional cell coords"
    g = gtensor(cell)

    h=g[0, 0]*x+g[1, 0]*y+g[2, 0]*z;
    k=g[0, 1]*x+g[1, 1]*y+g[2, 1]*z;
    l=g[0, 2]*x+g[1, 2]*y+g[2, 2]*z;
    return h, k, l

def S2R(qx, qy, qz, cell):
    "Given cartesian coordinates of a vector in the S System, calculate its Miller indexes."
    #definire x y e z
    H=qx*x[0, :]+qy*y[0, :]+qz*z[0, :];
    K=qx*x[1, :]+qy*y[1, :]+qz*z[1, :];
    L=qx*x[2, :]+qy*y[2, :]+qz*z[2, :];
    q=np.sqrt(qx**2+qy**2+qz**2);
    return H, K, L, q

def R2S(x, y, z, H, K, L):
    "Given reciprocal-space coordinates of a vecotre, calculate its coordinates in the Cartesian space."
    qx=self.scalar(H, K, L, x[0, :], x[1, :], x[2, :], 'latticestar');
    qy=self.scalar(H, K, L, y[0, :], y[1, :], y[2, :], 'latticestar');
    qz=self.scalar(H, K, L, z[0, :], z[1, :], z[2, :], 'latticestar');
    q=self.modvec(H, K, L, 'latticestar');
    return qx, qy, qz, q


if __name__ == "__main__":
    import doctest
    doctest.testmod()
