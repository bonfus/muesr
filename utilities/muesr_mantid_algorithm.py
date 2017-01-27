

from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import (LoadMuonNexus, DeleteWorkspace,
                              LoadDetectorsGroupingFile, CreateWorkspace,
                              GroupWorkspaces, CreateEmptyTableWorkspace)

# there seem to be different names depending on the mantid version
try:
    from mantidplot import proxies
except:
    from pymantidplot import proxies

#python2 or python3?
try:
    from queue import Queue, Empty
except:
    from Queue import Queue, Empty

#this is useful if you don't want to install muesr
import sys, os
sys.path.insert(0,os.path.dirname(os.path.realpath(__file__)))


#from PyQt4 import QtCore, QtGui
from muesr.gui.qt.muesr_gui import MainWindow
from muesr.core.sample import Sample as Smpl
from muesr.core.atoms import Atoms

#from muesrgui import MainWindow


#grab data from Mantid

class MuESRTool(PythonAlgorithm):
 
    def category(self):
        return 'Muon'
    
    def summary(self):
        return "Run MUESR from within MantidPlot"

    def PyInit(self):
        self.declareProperty(WorkspaceProperty(name="Input_Workspace",
                                               defaultValue="",
                                               direction=Direction.Input,
                                               optional=PropertyMode.Optional),
                             doc="Definition of sample structure")    
    def PyExec(self):
        ws = self.getProperty('Input_Workspace').value
        start_gui(ws)


AlgorithmFactory.subscribe(MuESRTool)




def queueSampleDefFromWS(in_q, ws):
    
    # nothing to do
    if not ws:
        return
    
    # workspaces may not have a crystal structure defined!
    try:
        sample = ws.sample()
    except: #what sort of error is raised?!
        return
    
    try:
        # parse lattice data and create a Ase atom
        crystalStruct = sample.getCrystalStructure()
        unitCell = crystalStruct.getUnitCell()
        atoms = crystalStruct.getScatterers()
        sg = crystalStruct.getSpaceGroup()
    except:
        #this is an error! Maybe Mantid apis changed?
        print("Could not parse lattice structure from Workspace.")
        return
        
        
    symbols = []
    scaled_positions = []
    
    
    
    for atm in atoms:
        data = atm.split()
        symbols.append(data[0])
        scaled_positions.append([float(x) for x in data[1:4]])
    
    cell = unitCell.getBinv()
                
    
    all_sym = []
    all_pos = []
    for i, p in enumerate(scaled_positions):
        eqps = sg.getEquivalentPositions(p)
        
        for eqp in eqps:
            all_sym.append(symbols[i])
            all_pos.append(eqp)
        
    atoms = Atoms(symbols=all_sym, 
                  scaled_positions=all_pos, 
                  cell=cell)
                  
    s = Smpl()
    s.cell = atoms 
    in_q.put(s)



def get_sample_definition_from_workspaces(in_q, ws=None):
    
    # if no workspace is given, all of them are loaded...
    if ws == None:
        # load info from all warkspaces
        for name in mtd.getObjectNames():
            if name == "":
                continue
                
            ws = None
            try:
                ws = mtd[name]
            except KeyError: #this is not the correct error! Should never happen anyway!
                print("something very odd")
                continue
            queueSampleDefFromWS(in_q, ws)
    else:
        queueSampleDefFromWS(in_q, ws)


def start_gui(ws = None):
    mtd.importAll()
    
    mantid_input_queue = Queue()
    mantid_output_queue = Queue()
    
    get_sample_definition_from_workspaces(mantid_input_queue, ws)
    
    muesr_gui = proxies.threadsafe_call(MainWindow,None,mantid_input_queue,mantid_output_queue)
    
    print("Starting MUESR GUI")
    proxies.threadsafe_call(muesr_gui.show)
    

    while True:
        
        try:
            data = mantid_output_queue.get(True, 1)
        except Empty:
            if not muesr_gui.isHidden():
                continue
            else:

                proxies.threadsafe_call(muesr_gui.deleteLater)
                muesr_gui = None

                break
        
        #not needed
        if not data:
            proxies.threadsafe_call(muesr_gui.deleteLater)
            muesr_gui = None
            break
        
        
        wpname = '_'.join(data[0].split()) #other restrictions on the names?
        tableWS = CreateEmptyTableWorkspace()
        
        
        tableWS.setTitle(data[0]) # Create a new table with 5 rows and 3 columns

        if len(data[1])>0:
            if data[1][0].T.ndim == 1:
                tableWS.addColumn(type="int",name="Muon site", plottype=1)
                tableWS.addColumn(type="double",name="Lorentz X", plottype=2)
                tableWS.addColumn(type="double",name="Lorentz Y", plottype=2)
                tableWS.addColumn(type="double",name= "Lorentz Z", plottype=2)
                tableWS.addColumn(type="double",name= "Dipolar X", plottype=2)
                tableWS.addColumn(type="double",name= "Dipolar Y", plottype=2)
                tableWS.addColumn(type="double",name= "Dipolar Z", plottype=2)
                tableWS.addColumn(type="double",name="Contact X", plottype=2)
                tableWS.addColumn(type="double",name="Contact Y", plottype=2)
                tableWS.addColumn(type="double",name= "Contact Z", plottype=2)
                tableWS.addColumn(type="double",name= "Total X", plottype=2)
                tableWS.addColumn(type="double",name= "Total Y", plottype=2)
                tableWS.addColumn(type="double",name= "Total Z", plottype=2)
                
                for i, locField in enumerate(data[1]):
                    nextRow = { 'Muon site': i +1,
                        'Lorentz X': locField.L[0],
                        'Lorentz Y': locField.L[1],
                        'Lorentz Z': locField.L[2],
                        'Dipolar X': locField.D[0],
                        'Dipolar Y': locField.D[1],
                        'Dipolar Z': locField.D[2],
                        'Contact X': locField.C[0],
                        'Contact Y': locField.C[1],
                        'Contact Z': locField.C[2],
                        'Total X': locField.T[0],
                        'Total Y': locField.T[1],
                        'Total Z': locField.T[2],
                         }
                    tableWS.addRow ( nextRow )
            else:
                print("Not supported yet")
    
