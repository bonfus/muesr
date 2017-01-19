from muesr.core.parsers import *
from muesr.core.cells import get_reduced_bases
from muesr.settings import config

import numpy as np

#@profile
def build_uniform_grid(sample, size, min_distance_from_atoms=1.0):
    """
    Generates a grid of symmetry inequivalent interstitial
    positions with a specified minimum distance from the atoms of the
    sample. Especially intended for DFT simulations.
    
    :param sample: A sample object.
    :param int size: The number of steps in the three lattice directions. 
                     Only equispaced grids are supported at the moment.
    :param float min_distance_from_atoms: Minimum distance between a 
                                          interstitial position and the
                                          atoms of the lattice.
                                          Units are Angstrom.
    :returns: A list of symmetry inequivalent positions.
    :rtype: list 
    """
    
    tolerance = config.FCRD
    
    #build uniform grid
    npoints = size**3
    x_ = np.linspace(0., 1., size, endpoint=False)
    y_ = np.linspace(0., 1., size, endpoint=False)
    z_ = np.linspace(0., 1., size, endpoint=False)

    uniform_grid = np.meshgrid(x_, y_, z_, indexing='ij')
    x,y,z = uniform_grid
    
    
    equiv=np.ones_like(x)*(npoints)
    
    
    alt = True
    nb_cells = np.array([[-1., -1., -1.],
                         [-1., -1.,  0.],
                         [-1., -1.,  1.],
                         [-1.,  0., -1.],
                         [-1.,  0.,  0.],
                         [-1.,  0.,  1.],
                         [-1.,  1., -1.],
                         [-1.,  1.,  0.],
                         [-1.,  1.,  1.],
                         [ 0., -1., -1.],
                         [ 0., -1.,  0.],
                         [ 0., -1.,  1.],
                         [ 0.,  0., -1.],
                         [ 0.,  0.,  0.],
                         [ 0.,  0.,  1.],
                         [ 0.,  1., -1.],
                         [ 0.,  1.,  0.],
                         [ 0.,  1.,  1.],
                         [ 1., -1., -1.],
                         [ 1., -1.,  0.],
                         [ 1., -1.,  1.],
                         [ 1.,  0., -1.],
                         [ 1.,  0.,  0.],
                         [ 1.,  0.,  1.],
                         [ 1.,  1., -1.],
                         [ 1.,  1.,  0.],
                         [ 1.,  1.,  1.]])

    
    for i in range(size):
        for j in range(size):
            for k in range(size):
                
                if equiv[i,j,k] < npoints:
                    #this point is equivalent to someone else!
                    continue
                
                for r,t in sample.sym.get_symop():
                    # new position for the muon
                    n = np.zeros(3)
                    # apply symmetry and bring back to unit cell
                    n = np.round(np.dot(r,[x[i,j,k],y[i,j,k],z[i,j,k]])+t,decimals=config.FCRD)%1
                    if (np.abs(n*size - np.rint(n*size)) < 10**-(config.FCRD)).all():
                        
                        #get index of point
                        ii,jj,kk=np.rint(n*size).astype(int)
                        if (ii*(size**2)+jj*size+kk > i*(size**2)+j*size+k):
                            equiv[ii,jj,kk] -= 1
                            equiv[i,j,k] += 1

    
    
    reduced_bases = sample.cell.get_cell()
    scaled_pos = sample.cell.get_scaled_positions()
                        
    positions = []
    for i in range(size):
        for j in range(size):
            for k in range(size):
                if equiv[i,j,k] >= npoints:
                    #saves distances with all atoms
                    distances = []
                    dists = np.zeros(27*len(sample._cell))
                    center = [x[i,j,k],y[i,j,k],z[i,j,k]]

                    #check distances form atoms (also in neighbouring cells)
                    

                    for a in range(len(sample._cell)):
                        dists[a*27:(a+1)*27] = np.linalg.norm(
                                                  np.dot(scaled_pos[a] - center + nb_cells,
                                                         reduced_bases), axis=1 )
                    #
                    # old method 20 times slower!
                    #
                    #for a in range(len(sample._cell)):
                    #    for ii in (-1, 0, 1):
                    #        for jj in (-1, 0, 1):
                    #            for kk in (-1, 0, 1):
                    #                distances.append( np.linalg.norm(
                    #                        np.dot(scaled_pos[a] - center + np.array([ii,jj,kk]),
                    #                            reduced_bases) ) )
                    #
                    #print(np.allclose(dists,distances))
                    #if min(distances) > min_distance_from_atoms:
                    
                    
                    if dists.min() > min_distance_from_atoms:
                        positions.append([x[i,j,k],y[i,j,k],z[i,j,k]])
                        
    return positions


if __name__ == "__main__":
    unittest.main()
        
