from muesr.core.nprint  import nprint
from re import split as resplit
import numpy as np

mybool = lambda x: x in ['true', '1', 't', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh']


def parse_plotstr(arg):
    plt_res_commands = ['list', 'set', 'show', 'ion', 'ioff', 'xlabel', 'ylabel', 'pylab', 'legend']
    if arg not in plt_res_commands:
        return arg.strip().replace(' ','_')
    else:
        raise ValueError

def parse_int(arg):
    'Convert a series of zero or more numbers to an argument tuple'
    if type(arg) is str:
        try:
            return tuple(map(int, arg.split() ))
        except:
            raise ValueError
    elif type(arg) is tuple:
        try:
            return tuple([int(x) for x in arg])
        except:
            raise ValueError            
    elif type(arg) is list:
        try:
            return tuple([int(x) for x in arg])
        except:
            raise ValueError
    else:
        raise ValueError

def parse_float(arg):
    'Convert a series of zero or more numbers to an argument tuple'
    if type(arg) is str:
        try:
            return tuple(map(float, arg.split() ))
        except:
            raise ValueError
    elif type(arg) is tuple:
        try:
            return tuple([float(x) for x in arg])
        except:
            raise ValueError            
    elif type(arg) is list:
        try:
            return tuple([float(x) for x in arg])
        except:
            raise ValueError
    else:
        raise ValueError

def parse_vector(arg, dimension = 3):

    choice = parse_float(arg)
        
    if len(choice) == 1 and choice[0] == 0:
        return ([0]*dimension)

    if len(choice) == dimension:
        return choice
    else:
        raise ValueError
        
def parse_complex_vector(arg, dimension = 6):

    choice = parse_float(arg)
        
    if len(choice) == 1 and choice[0] == 0:
        return ([0]*dimension)

    elif len(choice) == 3:
        return (choice+(0.,)*int(dimension/2))

    elif len(choice) == dimension:
        return choice
    else:
        raise ValueError

def parse_bool(arg):
    return mybool(arg.lower())


