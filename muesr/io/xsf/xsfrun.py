from muesr.settings import config
import subprocess

def run_xcrysden(fname):
    p = subprocess.Popen([config.XCrysExec, "--xsf", fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

