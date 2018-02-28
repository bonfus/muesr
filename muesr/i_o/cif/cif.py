"""
This file was adapted from ASE.
Module to read and write atoms in cif file format.

See http://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for a
description of the file format.  STAR extensions as save frames,
global blocks, nested loops and multi-data values are not supported.
"""

import re, os
import shlex
import warnings

import numpy as np

from muesr.core.magmodel import MM
from muesr.core.spg import Spacegroup
from muesr.core.nprint import nprint, nprintmsg
from muesr.core.isstr import isstr

from muesr.i_o.cif.crystal import crystal
from muesr.core.spg import spacegroup_from_data
from muesr.i_o.cif.cell import cellpar_to_cell

# Old conventions:
old_spacegroup_names = {'Abm2': 'Aem2',
                        'Aba2': 'Aea2',
                        'Cmca': 'Cmce',
                        'Cmma': 'Cmme',
                        'Ccca': 'Ccc1'}

        
def load_cif(sample, filename, reset_muon=True, reset_sym=True):
    """

    Loads the structural information from a Cif file.
    
    
    :param muesr.core.sample.Sample sample: the sample object.
    :param str filename:   the cif file path (path + filename).
    :param bool reset_muon: if true the muon positions is reinitialized.
    :param bool reset_sym:  if true the symmetry is reinitialized. ALWAYS TRUE
    :returns: True if succesfull, false otherwise
    :rtype: bool
    """

    filename = str(filename)
    
    f = open(os.path.expanduser(filename),'r')
    
        #nprint ("Problems with file name...file not found!", 'warn')
        #return False
    
    #print(read_cif(f,0))
    atoms, sym = read_cif(f,0) # selectd index 0
    f.close()
    if atoms:
        sample._reset(cell=True,sym=True,magdefs=True,muon=reset_muon)
        sample.cell = atoms
        sample.sym = sym
        return True
    else:
        nprint ("Atoms not loaded!", 'warn')
        return False

def load_mcif(sample, filename, reset_muon=True, reset_sym=True):
    """
    Loads both the crystalline structure and the magnetic order from a
    mcif file.
    N.B.: This function is EXPERIMENTAL.
    
    .. note::
       Only the first lattice and magnetic structure is loaded.
       mCif files hosting multiple structures are not supported.
    
    :param muesr.core.sample.Sample sample: the sample object.
    :param str filename:   the mcif file path (path + filename).
    :param bool reset_muon: if true the muon positions is reinitialized.
    :param bool reset_sym:  if true the symmetry is reinitialized.
    :returns: True if succesfull, false otherwise
    :rtype: bool
    """
    
    
    # DEFINITION OF UNITS AND SETTINGS: http://cmswiki.byu.edu/wiki/Magnetic_Coordinates
    #   new link http://magcryst.org/resources/magnetic-coordinates/
    
    f = open(os.path.expanduser(filename),'r')
    
    data = parse_cif(f)
    
    f.close()
    
    #print('Parsing data from: ' + data[0][0])
    
    tags = data[0][1]
    
    # Load symmetry
    #sym = tags2spacegroup(tags)
    
    # load cell info
    a = tags['_cell_length_a']
    b = tags['_cell_length_b']
    c = tags['_cell_length_c']
    alpha = tags['_cell_angle_alpha']
    beta = tags['_cell_angle_beta']
    gamma = tags['_cell_angle_gamma']

    
    # Find magnetic atoms
    mag_atoms_labels = tags['_atom_site_moment_label']
    # load mag moments
    mag_atoms_moments = np.zeros([len(mag_atoms_labels),3],dtype=np.complex)

    scaled_positions = np.array([tags['_atom_site_fract_x'], 
                                tags['_atom_site_fract_y'], 
                                tags['_atom_site_fract_z']]).T
    
    scaled_positions = np.mod(scaled_positions, 1.)
    
    all_scaled_pos = np.copy(scaled_positions)
    

    
    symbols = []
    if '_atom_site_type_symbol' in tags:
        labels = tags['_atom_site_type_symbol']
    else:
        labels = tags['_atom_site_label']
    for s in labels:
        # Strip off additional labeling on chemical symbols
        m = re.search(r'([A-Z][a-z]?)', s)  
        symbol = m.group(0)
        symbols.append(symbol)
    
    
    
    fcs = np.zeros_like(all_scaled_pos,dtype=np.complex)
    
    for i, al in enumerate(tags['_atom_site_label']):
        if al in tags['_atom_site_moment_label']:
            mi = tags['_atom_site_moment_label'].index(al)
            fcs[i] = [complex(tags['_atom_site_moment_crystalaxis_x'][mi]),
                        complex(tags['_atom_site_moment_crystalaxis_y'][mi]),  
                        complex(tags['_atom_site_moment_crystalaxis_z'][mi])
                        ]
    # THESE ARE IN CRYSTAL AXIS COORDINATE SYSTEM!!!
    # bohr magneton units are used
    # the magnetic metric tensor is M = L.G.L^(-1), which is unitless. 
    # NOW GO TO REDUCED LATTICE COORDINATE SYSTEM TO DO THE SYMMETRY
    L = np.diag([1./a,1./b,1./c])
    fcs = np.dot(fcs,L)
    
    # we copy the fourier components that were present in the mcif.
    # the others will be obtained from symmetry operations
    all_fcs = np.copy(fcs)
    
    for j, m_a_p in enumerate(scaled_positions):
        for cent in tags['_space_group_symop.magn_centering_xyz']: # magn_centering_xyz
            rc,tc,trc = parse_magn_operation_xyz_string(cent)
            # apply centering and go back to the unit cell
            cm_a_p = (np.dot(rc,m_a_p)+tc)%1.
            
            #print ('cmap is :' + str(cm_a_p))
            for s in tags['_space_group_symop.magn_operation_xyz']: # magn_operation_xyz
                
                r,t,tr = parse_magn_operation_xyz_string(s)
                
                symp = (np.dot(r,cm_a_p)+t)%1.
                
                #print('Symp is: '+ str(symp))
                # check if this position was already present
                for l, pp in enumerate(all_scaled_pos):
                    if np.allclose(symp,pp,rtol=1e-3):
                        break
                else:
                    # Append position just found with symmetry
                    symbols.append(symbols[j])
                    all_scaled_pos = np.append(all_scaled_pos,[symp],axis=0)
                    
                    # go to crystal units
                    crysfc = trc*np.linalg.det(rc)*tr*np.linalg.det(r)*np.dot(r, np.dot(rc,fcs[j]))                   
                    
                    all_fcs = np.append(all_fcs ,[crysfc],axis=0)
    


    # 
    mag_crys2car = cellpar_to_cell([a, b, c, alpha, beta, gamma], (0,0,1), None)
    # go to cartesian coordinates
    cartfc = np.dot(all_fcs,mag_crys2car)

    
    sample._reset(muon=reset_muon,sym=reset_sym)
    # symmetry is already introduced when parsing magnetism
    sample.cell, _ = crystal(symbols=symbols, basis=all_scaled_pos, cellpar=[a, b, c, alpha, beta, gamma])
    
    # initialization needs the number of atoms in the unit cell
    nmm=MM(len(all_scaled_pos),sample._cell.get_cell())
    # magnetic moments are specified in cartesian coordinates.
    # position in crystal coordinates...not nice but simple!
    nmm.fc_set(cartfc)
    # propagation vector is 0 since the cell "contains" the magnetic 
    #   structure, no incommensurate structures, sorry! 
    nmm.k=np.array([0.,0.,0.])
    sample.mm=nmm
    
    return True

        

def convert_value(value):
    """Convert CIF value string to corresponding python type."""
    value = value.strip()
    if re.match('(".*")|(\'.*\')$', value):
        return value[1:-1]
    elif re.match(r'[+-]?\d+$', value):
        return int(value)
    elif re.match(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$', value):
        return float(value)
    elif re.match(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\(\d+\)$',
                  value):
        return float(value[:value.index('(')])  # strip off uncertainties
    elif re.match(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\(\d+$',
                  value):
        warnings.warn('Badly formed number: "{0}"'.format(value))
        return float(value[:value.index('(')])  # strip off uncertainties
    else:
        return value


def parse_multiline_string(lines, line):
    """Parse semicolon-enclosed multiline string and return it."""
    assert line[0] == ';'
    strings = [line[1:].lstrip()]
    while True:
        line = lines.pop().strip()
        if line[:1] == ';':
            break
        strings.append(line)
    return '\n'.join(strings).strip()


def parse_singletag(lines, line):
    """Parse a CIF tag (entries starting with underscore). Returns
    a key-value pair."""
    kv = line.split(None, 1)
    if len(kv) == 1:
        key = line
        line = lines.pop().strip()
        while not line or line[0] == '#':
            line = lines.pop().strip()
        if line[0] == ';':
            value = parse_multiline_string(lines, line)
        else:
            value = line
    else:
        key, value = kv
    return key, convert_value(value)


def parse_loop(lines):
    """Parse a CIF loop. Returns a dict with column tag names as keys
    and a lists of the column content as values."""
    header = []
    line = lines.pop().strip()
    while line.startswith('_'):
        tokens = line.split()
        header.append(tokens[0].lower())
        if len(tokens) == 1:
            line = lines.pop().strip()
        else:
            line = ' '.join(tokens[1:])
            break
    columns = dict([(h, []) for h in header])
    if len(columns) != len(header):
        seen = set()
        dublicates = [h for h in header if h in seen or seen.add(h)]
        warnings.warn('Duplicated loop tags: {0}'.format(dublicates))

    tokens = []
    while True:
        lowerline = line.lower()
        if (not line or
            line.startswith('_') or
            lowerline.startswith('data_') or
            lowerline.startswith('loop_')):
            break
        if line.startswith('#'):
            line = lines.pop().strip()
            continue
        if line.startswith(';'):
            t = [parse_multiline_string(lines, line)]
        else:
            if len(header) == 1:
                t = [line]
            else:
                t = shlex.split(line, posix=False)

        line = lines.pop().strip()

        tokens.extend(t)
        if len(tokens) < len(columns):
            continue
        if len(tokens) == len(header):
            for h, t in zip(header, tokens):
                columns[h].append(convert_value(t))
        else:
            warnings.warn('Wrong number of tokens: {0}'.format(tokens))
        tokens = []
    if line:
        lines.append(line)
    return columns


def parse_items(lines, line):
    """Parse a CIF data items and return a dict with all tags."""
    tags = {}
    while True:
        if not lines:
            break
        line = lines.pop()
        if not line:
            break
        line = line.strip()
        lowerline = line.lower()
        if not line or line.startswith('#'):
            continue
        elif line.startswith('_'):
            key, value = parse_singletag(lines, line)
            tags[key.lower()] = value
        elif lowerline.startswith('loop_'):
            tags.update(parse_loop(lines))
        elif lowerline.startswith('data_'):
            if line:
                lines.append(line)
            break
        elif line.startswith(';'):
            parse_multiline_string(lines, line)
        else:
            raise ValueError('Unexpected CIF file entry: "{0}"'.format(line))
    return tags


def parse_block(lines, line):
    """Parse a CIF data block and return a tuple with the block name
    and a dict with all tags."""
    assert line.lower().startswith('data_')
    blockname = line.split('_', 1)[1].rstrip()
    tags = parse_items(lines, line)
    return blockname, tags


def parse_cif(fileobj):
    """Parse a CIF file. Returns a list of blockname and tag
    pairs. All tag names are converted to lower case."""
    if isstr(fileobj):
        fileobj = open(fileobj)
    lines = [''] + fileobj.readlines()[::-1]  # all lines (reversed)
    blocks = []
    while True:
        if not lines:
            break
        line = lines.pop()
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        blocks.append(parse_block(lines, line))
    return blocks


def tags2atoms(tags, store_tags=False, primitive_cell=False,
               subtrans_included=True):
    """Returns an Atoms object from a cif tags dictionary.  See read_cif()
    for a description of the arguments."""
    if primitive_cell and subtrans_included:
        raise RuntimeError(
            'Primitive cell cannot be determined when sublattice translations '
            'are included in the symmetry operations listed in the CIF file, '
            'i.e. when `subtrans_included` is True.')

    a = tags['_cell_length_a']
    b = tags['_cell_length_b']
    c = tags['_cell_length_c']
    alpha = tags['_cell_angle_alpha']
    beta = tags['_cell_angle_beta']
    gamma = tags['_cell_angle_gamma']

    scaled_positions = np.array([tags['_atom_site_fract_x'],
                                 tags['_atom_site_fract_y'],
                                 tags['_atom_site_fract_z']]).T

    symbols = []
    if '_atom_site_type_symbol' in tags:
        labels = tags['_atom_site_type_symbol']
    else:
        labels = tags['_atom_site_label']
    for s in labels:
        # Strip off additional labeling on chemical symbols
        m = re.search(r'([A-Z][a-z]?)', s)
        symbol = m.group(0)
        symbols.append(symbol)

    # Symmetry specification, see
    # http://www.iucr.org/resources/cif/dictionaries/cif_sym for a
    # complete list of official keys.  In addition we also try to
    # support some commonly used depricated notations
    no = None
    if '_space_group.it_number' in tags:
        no = tags['_space_group.it_number']
    elif '_space_group_it_number' in tags:
        no = tags['_space_group_it_number']
    elif '_symmetry_int_tables_number' in tags:
        no = tags['_symmetry_int_tables_number']

    symbolHM = None
    if '_space_group.Patterson_name_h-m' in tags:
        symbolHM = tags['_space_group.patterson_name_h-m']
    elif '_symmetry_space_group_name_h-m' in tags:
        symbolHM = tags['_symmetry_space_group_name_h-m']
    elif '_space_group_name_h-m_alt' in tags:
        symbolHM = tags['_space_group_name_h-m_alt']

    if symbolHM is not None:
        symbolHM = old_spacegroup_names.get(symbolHM.strip(), symbolHM)

    for name in ['_space_group_symop_operation_xyz',
                 '_space_group_symop.operation_xyz',
                 '_symmetry_equiv_pos_as_xyz']:
        if name in tags:
            sitesym = tags[name]
            break
    else:
        sitesym = None

    spacegroup = 1
    if sitesym is not None:
        subtrans = [(0.0, 0.0, 0.0)] if subtrans_included else None
        spacegroup = spacegroup_from_data(
            no=no, symbol=symbolHM, sitesym=sitesym, subtrans=subtrans)
    elif no is not None:
        spacegroup = no
    elif symbolHM is not None:
        spacegroup = symbolHM
    else:
        spacegroup = 1

    if store_tags:
        kwargs = {'info': tags.copy()}
    else:
        kwargs = {}

    if 'D' in symbols:
        deuterium = [symbol == 'D' for symbol in symbols]
        symbols = [symbol if symbol != 'D' else 'H' for symbol in symbols]
    else:
        deuterium = False

    atoms, spacegroup = crystal(symbols, basis=scaled_positions,
                                cellpar=[a, b, c, alpha, beta, gamma],
                                spacegroup=spacegroup, primitive_cell=primitive_cell,
                                **kwargs)
    if deuterium:
        masses = atoms.get_masses()
        masses[atoms.numbers == 1] = 1.00783
        masses[deuterium] = 2.01355
        atoms.set_masses(masses)

    return atoms, spacegroup


def read_cif(fileobj, index, store_tags=False, primitive_cell=False,
             subtrans_included=True):
    """Read Atoms object from CIF file. *index* specifies the data
    block number or name (if string) to return.

    If *index* is None or a slice object, a list of atoms objects will
    be returned. In the case of *index* is *None* or *slice(None)*,
    only blocks with valid crystal data will be included.

    If *store_tags* is true, the *info* attribute of the returned
    Atoms object will be populated with all tags in the corresponding
    cif data block.

    If *primitive_cell* is true, the primitive cell will be built instead
    of the conventional cell.

    If *subtrans_included* is true, sublattice translations are
    assumed to be included among the symmetry operations listed in the
    CIF file (seems to be the common behaviour of CIF files).
    Otherwise the sublattice translations are determined from setting
    1 of the extracted space group.  A result of setting this flag to
    true, is that it will not be possible to determine the primitive
    cell.
    """
    blocks = parse_cif(fileobj)
    # Find all CIF blocks with valid crystal data
    images = []
    for name, tags in blocks:
        try:
            atoms, spg = tags2atoms(tags, store_tags, primitive_cell,
                               subtrans_included)
            images.append([atoms,spg])
        except KeyError:
            pass
    for data in images[index]:
        yield data


def split_chem_form(comp_name):
    """Returns e.g. AB2  as ['A', '1', 'B', '2']"""
    split_form = re.findall(r'[A-Z][a-z]*|\d+',
                            re.sub('[A-Z][a-z]*(?![\da-z])',
                                   r'\g<0>1', comp_name))
    return split_form


def write_cif(fileobj, images, format='default'):
    """Write *images* to CIF file."""
    if isstr(fileobj):
        fileobj = paropen(fileobj, 'w')

    if hasattr(images, 'get_positions'):
        images = [images]

    for i, atoms in enumerate(images):
        fileobj.write('data_image%d\n' % i)

        a, b, c, alpha, beta, gamma = atoms.get_cell_lengths_and_angles()

        if format == 'mp':

            comp_name = atoms.get_chemical_formula(mode='reduce')
            sf = split_chem_form(comp_name)
            formula_sum = ''
            ii = 0
            while ii < len(sf):
                formula_sum = formula_sum + ' ' + sf[ii] + sf[ii + 1]
                ii = ii + 2

            formula_sum = str(formula_sum)
            fileobj.write('_chemical_formula_structural       %s\n' %
                          atoms.get_chemical_formula(mode='reduce'))
            fileobj.write('_chemical_formula_sum      "%s"\n' % formula_sum)

        fileobj.write('_cell_length_a       %g\n' % a)
        fileobj.write('_cell_length_b       %g\n' % b)
        fileobj.write('_cell_length_c       %g\n' % c)
        fileobj.write('_cell_angle_alpha    %g\n' % alpha)
        fileobj.write('_cell_angle_beta     %g\n' % beta)
        fileobj.write('_cell_angle_gamma    %g\n' % gamma)
        fileobj.write('\n')

        if atoms.pbc.all():
            fileobj.write('_symmetry_space_group_name_H-M    %s\n' % '"P 1"')
            fileobj.write('_symmetry_int_tables_number       %d\n' % 1)
            fileobj.write('\n')

            fileobj.write('loop_\n')
            fileobj.write('  _symmetry_equiv_pos_as_xyz\n')
            fileobj.write("  'x, y, z'\n")
            fileobj.write('\n')

        fileobj.write('loop_\n')

        if format == 'mp':
            fileobj.write('  _atom_site_type_symbol\n')
            fileobj.write('  _atom_site_label\n')
            fileobj.write('   _atom_site_symmetry_multiplicity\n')
            fileobj.write('  _atom_site_fract_x\n')
            fileobj.write('  _atom_site_fract_y\n')
            fileobj.write('  _atom_site_fract_z\n')
            fileobj.write('  _atom_site_occupancy\n')
        else:
            fileobj.write('  _atom_site_label\n')
            fileobj.write('  _atom_site_occupancy\n')
            fileobj.write('  _atom_site_fract_x\n')
            fileobj.write('  _atom_site_fract_y\n')
            fileobj.write('  _atom_site_fract_z\n')
            fileobj.write('  _atom_site_thermal_displace_type\n')
            fileobj.write('  _atom_site_B_iso_or_equiv\n')
            fileobj.write('  _atom_site_type_symbol\n')

        scaled = atoms.get_scaled_positions()
        no = {}
        for i, atom in enumerate(atoms):
            symbol = atom.symbol
            if symbol in no:
                no[symbol] += 1
            else:
                no[symbol] = 1
            if format == 'mp':
                fileobj.write(
                    '  %-2s  %4s  %4s  %7.5f  %7.5f  %7.5f  %6.1f\n' %
                    (symbol, symbol + str(no[symbol]), 1,
                     scaled[i][0], scaled[i][1], scaled[i][2], 1.0))
            else:
                fileobj.write(
                    '  %-8s %6.4f %7.5f  %7.5f  %7.5f  %4s  %6.3f  %s\n' %
                    ('%s%d' % (symbol, no[symbol]),
                     1.0,
                     scaled[i][0],
                     scaled[i][1],
                     scaled[i][2],
                     'Biso',
                     1.0,
                     symbol))

#### Addition for magnetic CIF files ####
def parse_magn_operation_xyz_string(tag):
    """ Parse the symmetry operations of the magnetic part of mcif """
    # this could probably be used also for cif
    r1,r2,r3,p=tag.split(',')
    r = np.zeros([3,3]) # rotations
    t = np.zeros([3])   # translations
    
    # compile traslation matrix
    t[0] = convert_to_float(r1[-4:]) if '/' in r1 else 0.
    t[1] = convert_to_float(r2[-4:]) if '/' in r2 else 0.
    t[2] = convert_to_float(r3[-4:]) if '/' in r3 else 0.
    
    # compile rotation matrix
    for i, line in enumerate([r1,r2,r3]):
        x=y=z=0.
        m = re.search(r'-?\d?x',line)
        if m:
            xs=m.group()
            x = 1. if ('x' in xs and not '-' in xs) else -1.
            m = re.search(r'\d',xs)
            if m:
                x *= float(m.group())

        m = re.search(r'-?\d?y',line)
        if m:
            ys=m.group()
            y = 1. if ('y' in ys and not '-' in ys) else -1.
            m = re.search(r'\d',ys)
            if m:
                y *= float(m.group())
                
        m = re.search(r'-?\d?z',line)
        if m:
            zs=m.group()
            z = 1. if ('z' in zs and not '-' in zs) else -1.
            m = re.search(r'\d',zs)
            if m:
                z *= float(m.group())                            

            
        r[i,0] = x
        r[i,1] = y
        r[i,2] = z
        
    return (r,t,float(p))
        
        
        


def convert_to_float(frac_str):
    "This function converts fractions to float, ex. -1/3 = -0.33333..."
    try:
        return float(frac_str)
    except ValueError:
        try:
            num, denom = frac_str.split('/')
        except ValueError:
            return None
        try:
            leading, num = num.split(' ')
        except ValueError:
            return float(num) / float(denom)        
        if float(leading) < 0:
            sign_mult = -1
        else:
            sign_mult = 1
        return float(leading) + sign_mult * (float(num) / float(denom))
