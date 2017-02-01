from muesr.settings import config
import subprocess, warnings

def run_xcrysden(fname):
    if config.XCrysExec == None:
        warnings.warn("XCrysDen executable not found. Check configs.")
        return
    p = subprocess.Popen([config.XCrysExec, "--xsf", fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

