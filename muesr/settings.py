## GLOBAL SETTINGS ##
import os, warnings
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
            self._cfg.set('Executables', 'vesta_exec', 'VESTA')

        if not ('Numerical' in self._cfg.sections()):
            self._need_config = True
            self._cfg.add_section('Numerical')
            self._cfg.set('Numerical', 'RoundingFractionalCoordinates', '7')

        if not ('Visualization' in self._cfg.sections()):
            self._need_config = True
            self._cfg.add_section('Visualization')
            self._cfg.set('Visualization', 'DefaultApplication', 'xcrysden')
            self._cfg.set('Visualization', 'xcrysden', self._cfg.get('Executables', 'xcrysden_exec'))
            self._cfg.set('Visualization', 'vesta', self._cfg.get('Executables', 'vesta_exec'))


        #Fractiona coordinate rounding_decimals
        try:
            self._FCRD = self._cfg.getint('Numerical',
                                           'RoundingFractionalCoordinates')
        except:
            raise ValueError('Cannot set value from config.')

        #Name of xcrysden executable file
        self._XCRSEXEC = ""
        try:
            self._XCRSEXEC = self._cfg.get('Visualization', 'xcrysden')
        except CP.NoOptionError:
            self._XCRSEXEC = self._cfg.get('Executables', 'xcrysden_exec')
            self._cfg.set('Visualization', 'xcrysden', self._XCRSEXEC)
            
        self._XCRSEXEC = self._which(self._XCRSEXEC)
        #Name of VESTA executable file.
        # Set it if not found to keep compatibility with old log files.
        try:
            self._VESTAEXEC = self._cfg.get('Visualization', 'vesta')
        except CP.NoOptionError:
            self._cfg.set('Executables', 'vesta_exec', 'VESTA')
            self._cfg.set('Visualization', 'vesta', 'VESTA')
            self._VESTAEXEC = self._cfg.get('Visualization', 'vesta')
            
        self._VESTAEXEC = self._which(self._VESTAEXEC)

        #Path for xcrysden temp files (with ending /!)
        self._XCRSTMP = self._cfg.get('Directories', 'xcrysden_tmp')
        #check directory exists
        if not os.path.exists(self._XCRSTMP):
            warnings.warn('Temp dir for XCrysDen files "' + self._XCRSTMP \
                          + '" does not exists. Changing to ' + tempfile.gettempdir())
            self._XCRSTMP = tempfile.gettempdir()
        #check is writable
        test_file = os.path.join(self._XCRSTMP, 'tmp.test')
        try:
            open(test_file, 'a').close()
        except:
            raise ValueError("Cannot write into '{}'. Change TMP directory.".format(self._XCRSTMP))
        
        # Set default editor
        try:
            self._DEFAULTVISAPP = self._cfg.get('Visualization', 'DefaultApplication')
        except CP.NoSectionError:
            self._cfg.add_section('Visualization')
            self._cfg.set('Visualization', 'DefaultApplication', 'xcrysden')
            self._cfg.set('Visualization', 'vesta', self._cfg.set('Executables', 'vesta_exec'))
            self._cfg.set('Visualization', 'xcrysden', self._cfg.set('Executables', 'xcrysden_exec'))
        except CP.NoOptionError:
            self._cfg.set('Visualization', 'DefaultApplication', 'xcrysden')
        
        self._DEFAULTVISAPP = self._cfg.get('Visualization', 'DefaultApplication')

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
    
    @XCrysExec.setter
    def XCrysExec(self, value):
        abs_path = self._which(value)
        if abs_path:
            self._XCRSEXEC = abs_path
            # Supporting both old and new config. To be dropped soon.
            self._cfg.set('Executables', 'xcrysden_exec', value)
            self._cfg.set('Visualization', 'xcrysden', value)
            self.store()
        else:
            raise ValueError("File '{}' not found or not in PATH.".format(value))
    
    @property
    def VESTAExec(self):
        return self._VESTAEXEC
    
    @VESTAExec.setter
    def VESTAExec(self, value):
        abs_path = self._which(value)
        if abs_path:
            self._VESTAEXEC = abs_path
            # Supporting both old and new config. To be dropped soon.
            self._cfg.set('Executables', 'vesta_exec', value)
            self._cfg.set('Visualization', 'vesta', value)
            self.store()
        else:
            raise ValueError("File '{}' not found or not in PATH.".format(value))
    
    @property
    def AllVisExecs(self):
        apps = []
        for k, v in self._cfg.items('Visualization'):
            if k.lower() == 'defaultapplication':
                continue
            app_path = self._which(v)
            if app_path:
                if k.lower() == self._DEFAULTVISAPP.lower():
                    apps.insert(0, app_path)
                else:
                    apps.append(app_path)
        return apps
            
    
    @property
    def XCrysTmp(self):
        """
        This is the directory where temporary structure files in XCrysDen
        format are stored.
        """

        return self._XCRSTMP
        
    @XCrysTmp.setter
    def XCrysTmp(self, value):
        """
        This is the directory where temporary structure files in XCrysDen
        format are stored.
        """

        if not os.path.exists(value):
            raise ValueError("Directory '{}' does not exists".format(value))

        if not os.path.isdir(value):
            raise ValueError("Path '{}' is not a directory".format(value))

        test_file = os.path.join(value, 'tmp.test')
        try:
            open(test_file, 'a').close()
        except:
            raise ValueError("Cannot write into '{}'".format(value))
        
        self._XCRSTMP = value
        self._cfg.set('Directories', 'xcrysden_tmp', value)
        self.store()

    @property
    def DefaultVisualizationApp(self):
        """
        Default application for visualization
        """
        return self._DEFAULTVISAPP

    @DefaultVisualizationApp.setter
    def DefaultVisualizationApp(self, value):
        """
        Default application for visualization
        """
        try:
            self._DEFAULTVISAPP = str(value).lower()
        except:
            raise TypeError("Cannot convert value to string")
        self.store()

    @property
    def FCRD(self):
        """
        Rounding used in parsing fractional coordinates.
        """
        return self._FCRD
    
    @FCRD.setter
    def FCRD(self, value):
        """
        Rounding used in parsing fractional coordinates.
        """
        # avoid problems with long type in python2
        try:
            v = value + 1
            value = int(value)
        except:
            raise TypeError("Value must be int")
        self._FCRD = value
        self._cfg.set('Numerical', 'RoundingFractionalCoordinates', str(value))
        self.store()

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
