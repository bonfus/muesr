import numpy as np
import re

from muesr.core.atoms import *
from muesr.core.nprint import nprint
from muesr.core.ninput import ninput_mt
from muesr.core.parsers import mybool



def write_xsf(fileobj, images, data=None):
    
    # this should work with unicode too.
    if not hasattr(fileobj, 'write'):
        fileobj = open(fileobj, 'w')
        
    if not isinstance(images, (list, tuple)):
        images = [images]

    #fileobj.write('ANIMSTEPS %d\n' % len(images))

    numbers = images[0].get_atomic_numbers()
    
    #pbc = images[0].get_pbc()
    pbc = (True,True,True)
    
    if pbc[2]:
        fileobj.write('CRYSTAL\n')
    elif pbc[1]:
        fileobj.write('SLAB\n')
    elif pbc[0]:
        fileobj.write('POLYMER\n')
    else:
        fileobj.write('MOLECULE\n')

    for n, atoms in enumerate(images):
        if True: #there should be a check concerning pbc
            fileobj.write('PRIMVEC %d\n' % (n + 1))
            cell = atoms.get_cell()
            for i in range(3):
                fileobj.write(' %.14f %.14f %.14f\n' % tuple(cell[i]))

        fileobj.write('PRIMCOORD %d\n' % (n + 1))

        # Write magnetic moments as forces for xcrysden

        forces = atoms.get_magnetic_moments()
        
        if not forces is None:
            # reduce absolute value for nice plotting
            forces *= 0.01

        pos = atoms.get_positions()

        fileobj.write(' %d 1\n' % len(pos))
        for a in range(len(pos)):
            fileobj.write(' %2d' % numbers[a])
            fileobj.write(' %20.14f %20.14f %20.14f' % tuple(pos[a]))
            if forces is None:
                fileobj.write('\n')
            else:
                fileobj.write(' %20.14f %20.14f %20.14f\n' % tuple(forces[a]))
            
    if data is None:
        fileobj.close()
        return

    fileobj.write('BEGIN_BLOCK_DATAGRID_3D\n')
    fileobj.write(' data\n')
    fileobj.write(' BEGIN_DATAGRID_3Dgrid#1\n')

    data = np.asarray(data)
    if data.dtype == complex:
        data = np.abs(data)

    shape = data.shape
    fileobj.write('  %d %d %d\n' % shape)

    cell = atoms.get_cell()
    origin = np.zeros(3)
    for i in range(3):
        if not pbc[i]:
            origin += cell[i] / shape[i]
    fileobj.write('  %f %f %f\n' % tuple(origin))

    for i in range(3):
        fileobj.write('  %f %f %f\n' %
                      tuple(cell[i] * (shape[i] + 1) / shape[i]))

    for x in range(shape[2]):
        for y in range(shape[1]):
            fileobj.write('   ')
            fileobj.write(' '.join(['%f' % d for d in data[x, y]]))
            fileobj.write('\n')
        fileobj.write('\n')

    fileobj.write(' END_DATAGRID_3D\n')
    fileobj.write('END_BLOCK_DATAGRID_3D\n')
    fileobj.close()
    return
    
def read_xsf_data(f):
    
    found_data = False
    
    def search_data_block(f):
        while f.readline().find("BEGIN_BLOCK_DATAGRID_3D") == -1:
            continue
        return True
        
    while search_data_block(f):
        nprint ("Found data block:")
        nprint("Desc" + f.readline())
        nprint("type" + f.readline())
        if ninput_mt("Use this? ", mybool):
            found_data = True
            break 
        
    if not found_data:
        nprint ("No data found, exiting...", 'warn')
        return None
    
    grid = map(int, f.readline().split())
    
    # Lattice vectors along lines
    orig = np.array(map(float, f.readline().split()), dtype=object)
    v1 = np.array(map(float, f.readline().split()), dtype=object)
    v2 = np.array(map(float, f.readline().split()), dtype=object)
    v3 = np.array(map(float, f.readline().split()), dtype=object)
    
    v1n = v1/grid[0]
    v2n = v2/grid[1]
    v3n = v3/grid[2]
    
    sdata = np.zeros(np.prod(grid))
    
    data_buffer = []
    
    nprint("Grid has " + str(np.prod(grid)) + " points.")
    
    #grid_steps = [int(x/10) for x in grid]
    #nprint("Setting steps throught grid to %i %i %i." % tuple(grid_steps) )

    
    nprint ("Parsing data...")
    
    for k in range(0, grid[2]):
        for j in range(0,grid[1]):
            for i in range(0, grid[0]):
                if data_buffer == []:
                    data_buffer = [float(x) for x in f.readline().split()] # a random number of data may be on the line

                sdata[i+grid[0]*j+grid[1]*grid[0]*k] = data_buffer[0] 
                data_buffer.pop(0) #empty the buffer
    return [sdata, [grid, orig, [v1n,v2n,v3n]]]


def read_xsf(fileobj, index=-1, read_data=False):
    if not hasattr(fileobj, 'read'):
        fileobj = open(fileobj,'r')

    readline = fileobj.readline
    while True:
        line = readline().strip()
        if line != '':
            if line[0] != '#':
                break
    
    if 'INFO' in line:
        while True:
            line = readline().strip()
            if line != '':
                if line[0] != '#':
                    if 'END_INFO' in line:
                        break    

        line = readline().strip()

    if 'ANIMSTEPS' in line:
        nimages = int(line.split()[1])
        line = readline().strip()
    else:
        nimages = 1

    if 'CRYSTAL' in line:
        pbc = True
    elif 'SLAB' in line:
        pbc = (True, True, False)
    elif 'POLYMER' in line:
        pbc = (True, False, False)
    elif 'DIM-GROUP' in line:
        line = readline().strip()
        if int(line.split(' ')[0]) == 3:
            pbc = True
        else:
            nprint ("WARNING: I'm missing something in your XCrysDen file.",'warn')
            pbc = False            
    else:
        nprint ("WARNING: I'm missing something in your XCrysDen file.",'warn')
        pbc = False

    images = []
    for n in range(nimages):
        cell = None
        if pbc:
            while True:
                line = readline().strip()
                if line != '':
                    if line[0] != '#':
                        break
            assert 'PRIMVEC' in line
            cell = []
            for i in range(3):
                cell.append([float(x) for x in readline().split()])
        else:
            nprint ("This program cannot work without PBC...create a fake cell.",'warn')
            return None

        while not 'PRIMCOORD' in readline().strip():
            continue
        natoms = int(readline().split()[0])
        symbols = []
        numbers = []
        positions = []
        for a in range(natoms):
            line = readline().split()
            try:
                numbers.append(int(line[0]))
            except:
                atm_name = re.sub("\d","",str(line[0]))
                symbols.append(atm_name)
            positions.append([float(x) for x in line[1:]])

        positions = np.array(positions)
        if len(positions[0]) == 3:
            forces = None
        else:
            positions = positions[:, :3]
            forces = positions[:, 3:]
        if len(symbols) != 0:
            image = Atoms(symbols, positions, cell=cell, pbc=pbc)
        elif len(numbers) != 0:
            image = Atoms(numbers = numbers, positions = positions, cell=cell, pbc=pbc)

        images.append(image)

    if read_data:
        data = read_xsf_data(fileobj)
        #line = readline()
        #assert 'BEGIN_BLOCK_DATAGRID_3D' in line
        #line = readline()
        #assert 'BEGIN_DATAGRID_3D' in line
        #
        #shape = [int(x) for x in readline().split()]
        #start = [float(x) for x in readline().split()]
        #
        #for i in range(3):
        #    readline()
        #    
        #n_data = shape[0]*shape[1]*shape[2]
        #data = np.array([float(readline())
        #                 for s in range(n_data)]).reshape(shape[::-1])
        #data = np.swapaxes(data, 0, 2)
        #
        fileobj.close()
        return images[index], data
        
    fileobj.close()
    return images[index]

