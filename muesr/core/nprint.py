from muesr.core.prettytable import PrettyTable

HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'

bcolors = {'ok': OKGREEN , 'end' : ENDC, 'warn' : WARNING, 'error' : FAIL}

messages = {'ecrystal': [ 'No lattice structure in workspace!','warn'], \
            'esym'    : [ 'No symmetries defined!','warn'], \
            'esupcell': [ 'No supercell structure in workspace!','warn'], \
            'emuon'   : [ 'Muon position not defined!','warn'], \
            'emag'    : [ 'No magnetic structures defined!','warn'], \
            'ema'     : [ 'Command argument missing/invalid!','warn'], \
            'efile'   : [ 'Error opening file!','warn'], \
            'NS'      : [ 'Nothing set!','warn'] 
            }


def cstring(stri,c):
    "returns a coloured string"
    if c in bcolors.keys():
        return bcolors[c] + stri + bcolors['end']
    else:
        raise RuntimeWarning("Color not defined")
        

def nprint(arg, color = None):
    if type(arg) != str:
        print("\t " + bcolors['ok'] + str(arg) + bcolors['end'])
    else:
        if color:
            print ('\t ' + bcolors[color] + arg + bcolors['end'])
        else:
            print ("\t %s" % arg)

def nprinttab(arg, header):
    x = PrettyTable(header, border=True)
    for line in arg:
        x.add_row(line)
    print(x)
    
    

def nprintmsg(arg):
    if arg in messages.keys():
        nprint(messages[arg][0],messages[arg][1])
    else:
        nprint ("Standard message NOT DEFINED!",'error')


def print_cell(cell, atom_type = None):
    # cell: cell to print
    # atomtype: function to select which atom to print
    
    symbols = cell.get_chemical_symbols()
    masses = cell.get_masses()
    magmoms = cell.get_magnetic_moments()
    lattice = cell.get_cell()
    nprint ("Lattice vectors:")
    nprint ("  a %20.15f %20.15f %20.15f" % tuple( lattice[0] ),'ok')
    nprint ("  b %20.15f %20.15f %20.15f" % tuple( lattice[1] ),'ok')
    nprint ("  c %20.15f %20.15f %20.15f" % tuple( lattice[2] ),'ok')
    nprint ("Atomic positions (fractional):")
    if atom_type == None:
        atom_type = lambda x: True
        
    for i, v in enumerate(cell.get_scaled_positions()):
        if atom_type(symbols[i]):
            # we should print this atom
            if magmoms == None:
                nprint ("%5d %-2s%18.14f%18.14f%18.14f %7.3f" % \
                    (i+1, symbols[i], v[0], v[1], v[2], masses[i]),'ok')
            else:
                nprint ("%5d %-2s%18.14f%18.14f%18.14f %7.3f  %s" % \
                    (i+1, symbols[i], v[0], v[1], v[2], masses[i], ', '.join('%.3f' % val.real for val in magmoms[i])),'ok')
