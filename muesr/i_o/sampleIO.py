import os, warnings
import numpy as np
from muesr.core.parsers  import parse_bool
from muesr.core.ninput   import *
from muesr.core.atoms    import Atoms
from muesr.core.sample   import Sample
from muesr.core.sampleErrors   import *
from muesr.core.spg      import spacegroup_from_data
from muesr.core.magmodel import MM


have_yaml = True


try:
    from yaml import load, dump
    try:
        from yaml import CLoader as Loader, CDumper as Dumper
    except ImportError:
        from yaml import Loader, Dumper
except:
    have_yaml = False



def save_sample(sample, filename="", fileobj=None, overwrite=False):
    """
    This function saves the sample provided in YAML format.
    
    :param sample: the sample object
    :param str filename: the filename used to store data.
    :param file fileobj: a file object used in place of filename.
    :param overwrite bool: if selected file should be overwritten.
    :return: None
    :rtype: None
    :raises: TypeError, ValueError, IsADirectoryError
    """        
    
    
    if not have_yaml:
        warnings.warn("Warning, YAML python package not present!\n" +
                      "You cannot use this function")
        return False

    if not isinstance(sample, Sample):
        raise TypeError('Sample argument must be a Sample object.')

    outdict = {}
    
    outdict['Name'] = sample.name
    
    lattdict = {}
    try:
        sample._check_lattice()
        lattdict['Cell'] = sample._cell.get_cell().tolist()
        lattdict['Symbols'] = sample._cell.get_chemical_symbols()
        lattdict['ScaledPositions'] = sample._cell.get_scaled_positions().tolist()
    except CellError:
        pass
        
        
    outdict['Lattice'] = lattdict
    
    muondict = {}
    try:
        sample._check_muon()
        muondict['Positions'] = [x.tolist() for x in sample.muons]
    except MuonError:
        pass
    
    outdict['Muon'] = muondict
    
    magdict = {}
    try:
        sample._check_magdefs()
        magdict['Size'] = sample.mm.size
        maglist = []
        for md in sample._magdefs:  #not nice, should use properties
            maglist.append(
                {'format':'b-c',
                    'k': md.k.tolist(),
                    'fc': np.hstack([md.fc.real,md.fc.imag]).tolist(),
                    'phi': md.phi.tolist(),
                    'desc': md.desc
                }
            )
            if not (md.lattice_params is None):
                maglist[-1]['lattice'] = md.lattice_params.tolist()
                
        magdict['Orders'] = maglist
    except MagDefError:
        pass
    
    outdict['MagneticOrders'] = magdict
    
    
    symdict = {}
    try:
        sample._check_sym()
        symdict['Number'] = sample.sym.no
        symdict['Symbol'] = sample.sym.symbol
        symdict['Rotations'] = sample.sym.rotations.tolist()
        symdict['Translations'] = sample.sym.translations.tolist()
    except SymmetryError:
        pass
        
        
    outdict['Symmetry'] = symdict

    output = dump(outdict, Dumper=Dumper)
    
    if fileobj is None:
        if filename == "":
            raise ValueError("Specify filename or provide a file object")
            
        if os.path.isfile(filename) and (not overwrite):
            warnings.warn('File not (over)written.', RuntimeWarning)
            return False
        
        
        with open(filename,'w') as f:
            f.write(output)
        return True
        
    else:
        fileobj.write(output)
        return True

    
def load_sample(filename="", fileobj=None):

    """
    This function load a sample from a file in YAML format.
        
    :param str filename: the filename used to store data.
    :param file fileobj: an optional file object. If specified, this
                         supersede the filename input.
    :return: a sample object
    :rtype: :py:class:`~Sample` object or None
    :raises: ValueError, FileNotFoundError, IsADirectoryError
    
    """      
    
    # fail if YAML is not available
    if not have_yaml:
        warnings.warn("Warning, YAML python package not present!")
        return
    
    sample = Sample()
    
    data = {}
    if fileobj is None:
        if filename == "":
            raise ValueError("Specify filename or File object")
        with open(filename,'r') as f:
            data = load(f, Loader=Loader)
    else:
        data = load(fileobj, Loader=Loader)

    if not(type(data) is dict):
        raise ValueError('Invalid data file. (problems with YAML?)')

    
    if data is None:
        raise ValueError('Invalid/empty data file. (problems with YAML?)')
    
    if 'Name' in data.keys():
        sample.name = str(data['Name'])
    
    if 'Lattice' in data.keys():
        l = data['Lattice']
        
        spos = None
        cpos = None
        if 'ScaledPositions' in l.keys():
            spos = np.array(l['ScaledPositions'])
        elif 'CartesianPositions' in l.keys():
            cpos = np.array(l['CartesianPositions'])
        
        cell = None
        if 'Cell' in l.keys():
            cell = np.array(l['Cell'])
            
        symbols = None
        if 'Symbols' in l.keys():
            symbols = l['Symbols']
        
        if (cell is None) or (symbols is None):
            warnings.warn('Cell not loaded!', RuntimeWarning)
        else:
            if not spos is None:
                sample.cell = Atoms(symbols = symbols, scaled_positions = spos, cell=cell, pbc=True)
            elif not cpos is None:
                sample.cell = Atoms(symbols = symbols, positions = cpos, cell=cell, pbc=True)
            else:
                warnings.warn('Cell not loaded!', RuntimeWarning)

    else:
        warnings.warn('Cell not loaded!', RuntimeWarning)

            
    if 'Muon' in data.keys():
        m = data['Muon']
        if 'Positions' in m:
            for p in m['Positions']:
                sample.add_muon(p)
        else:
            warnings.warn('Muon positions not loaded!', RuntimeWarning)
    else:
        warnings.warn('Muon positions not loaded!', RuntimeWarning)
                
    if 'Symmetry' in data.keys():
        s = data['Symmetry']
        
        if ('Number' in s.keys()) and \
            ('Symbol' in s.keys()) and \
            ('Rotations' in s.keys()) and \
            ('Translations' in s.keys()):
            sample.sym = spacegroup_from_data(s['Number'],s['Symbol'],
                                            rotations=np.array(s['Rotations']),
                                            translations=np.array(s['Translations'])) 
        else:
            warnings.warn('Symmetry not loaded.', RuntimeWarning)
    else:
        warnings.warn('Symmetry not loaded!', RuntimeWarning)
            
    if 'MagneticOrders' in data.keys():
        m = data['MagneticOrders']
        
        if len(m) > 0:
        
            msize = int(m['Size'])
    
            for mo in m['Orders']:
                
                if 'lattice' in mo.keys():                
                    n = MM(msize, \
                            np.array(mo['lattice']))
                else:
                    n = MM(msize)
                
                sample.mm = n
                sample.mm.k = np.array(mo['k'])
                sample.mm.phi = np.array(mo['phi'])
                
                if 'desc' in mo.keys():
                    sample.mm.desc = str(mo['desc'])
                
                rfcs, ifcs = np.hsplit(np.array(mo['fc']),2)
                
            
                if mo['format'].lower() in ['bohr-cartesian', 'b-c']:
                    sample.mm.fcCart=(rfcs + 1.j*ifcs)
                    
                elif mo['format'].lower() in ['bohr/angstrom-lattic', 'b/a-l']:
                    sample.mm.fcLattBMA=(rfcs + 1.j*ifcs)
                    
                elif mo['format'].lower() in ['bohr-lattice','b-l']:
                    sample.mm.fcLattBM=(rfcs + 1.j*ifcs)
                    
                else:
                    raise ValueError('Invalid Fourier Components format specifier in YAML file.')
        else:
            warnings.warn('Magnetic definitions not loaded!', RuntimeWarning)
    else:
        warnings.warn('Magnetic definitions not loaded!', RuntimeWarning)
                
    return sample
        
        
if __name__ == '__main__':
    unittest.main()
        
