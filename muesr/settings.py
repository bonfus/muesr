## GLOBAL SETTINGS ##
import os
import tempfile
try:
    from appdirs import *
except:
    from muesr.core.appdirs import *  # Use own version if system version
    # is missing
try:
    import configparser as CP
except:
    import ConfigParser as CP  # Python 2


class Settings(object):
    def __init__(self):
        self._need_config = False

        d = user_config_dir('muesr')
        # create config dir if missing
        if not os.path.exists(d):
            os.makedirs(d)

        config_path = os.path.join(d, 'muesr.cfg')
        if not os.path.isfile(config_path):
            open(config_path, 'a').close()

        self._cfg = CP.ConfigParser()
        self._cfg.read(config_path)
        if not ('Directories' in self._cfg.sections()):
            self._need_config = True
            self._cfg.add_section('Directories')
            self._cfg.set('Directories', 'xcrysden_tmp', tempfile.gettempdir())

        if not ('Executables' in self._cfg.sections()):
            self._need_config = True
            self._cfg.add_section('Executables')
            self._cfg.set('Executables', 'xcrysden_exec', 'xcrysden')

        if not ('Numerical' in self._cfg.sections()):
            self._need_config = True
            self._cfg.add_section('Numerical')
            self._cfg.set('Numerical', 'RoundingFractionalCoordinates', '7')

        #Fractiona coordinate rounding_decimals
        try:
            self._FCRD = int(self._cfg.get('Numerical',
                                           'RoundingFractionalCoordinates'))
        except:
            raise ValueError('Cannot set value from config.')

        #Name of xcrysden executable file
        self._XCRSEXEC = self._cfg.get('Executables', 'xcrysden_exec')
        self._XCRSEXEC = self._which(self._XCRSEXEC)

        #Path for xcrysden temp files (with ending /!)
        self._XCRSTMP = self._cfg.get('Directories', 'xcrysden_tmp')
        #check directory exists
        if not os.path.exists(self._XCRSTMP):
            raise ValueError('Temp dir for XCrysDen does not exists.')
        #check is writable

    def store(self):
        d = user_config_dir('muesr')
        # create config dir if missing
        if not os.path.exists(d):
            os.makedirs(d)

        config_path = os.path.join(d, 'muesr.cfg')
        try:
            with open(config_path, 'w') as f:
                self._cfg.write(f)
        except EnvironmentError:  # parent of IOError, OSError *and* WindowsError where available
            print('Cannot save configs')

    @property
    def XCrysExec(self):
        return self._XCRSEXEC

    @property
    def XCrysTmp(self):
        return self._XCRSTMP

    @property
    def FCRD(self):
        return self._FCRD

    # from http://stackoverflow.com/a/377028
    def _which(self, program):
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file

        return None


config = Settings()
