import sys, os, threading
try:
    from queue import Queue, Empty
except ImportError:
    from Queue import Queue, Empty
from copy import copy
import logging as log
import numpy as np

# import muesr stuff
from muesr.core.sample import Sample
from muesr.core.sampleErrors import CellError, MuonError, MagDefError
from muesr.engines.clfc import locfield, find_largest_sphere
from muesr.io.xsf.xsf import load_xsf, show_cell, show_supercell
from muesr.io.cif.cif import load_cif
from muesr.io.sampleIO import load_sample, save_sample
from muesr.io.exportFPS import export_fpstudio

# tool to check is the interface is started as a Mantid Algorithm
from muesr.gui.mantidtools import running_inside_mantid

# Provides plots if matplotlib is available
from muesr.gui.qt.muesr_gui_matplotlib_interface import DataPlotDialog

# icons...just for fun
from muesr.gui.icons import get_icon

# mantid rins QT4 as of the time of writing...maybe update this in the 
# future
if running_inside_mantid:
    os.environ['QT_PREFERRED_BINDING'] = os.pathsep.join(['PyQt4'])


# import all needed QT stuff
from .Qt import QtCore, QtGui, QtWidgets
from .Qt.QtCore import QAbstractItemModel, QModelIndex

# constants affecting visulization

CONST_MAX_MUON_SITES_IN_EDITOR=12*12*12
CONST_MAX_SC_CELL_FOR_XCRYSDEN=5


class MagnModel(QAbstractItemModel):
    """
    An Item Model for reading and setting Fourier Components.
    Validates input and reacts on Coordinate system changes.
    """
    def __init__(self, parent=None, sample=None, coord = 0):
        QAbstractItemModel.__init__(self, parent)

        self.__headers = ["Atom", "Position", "Re X", "Re Y", "Re Z", "Im X", "Im Y", "Im Z"]

        if (coord >= 0 and coord < 3):
            self.__coord = coord
        else:
            raise ValueError
        
        # saves Chemical symbols and positions to help input
        # (sample variable is lost after initialization)
        self.__chem_symbols = sample.cell.get_chemical_symbols()
        self.__pos = sample.cell.get_scaled_positions()
        
        # a copy of the magnetic order is stored to simplify the 
        # handling of the coordinate systems.
        # 
        self.__mm = copy(sample.mm) #should this be a deepcopy?
        
        # get current structure in selected coordinates
        self.__mmdata = self.__mm.fc_get(coord)
        
        # save a backup copy to reset when user tries to store invalid
        # values
        self.__mmbkp = (self.__mm.fc_get(coord),coord)
        
        
        # needed to highlight the changed (bold font is displayed)
        self.__changed = [[False,]*6 for _ in range(len(self.__chem_symbols))]


    def setSample(self, sample):
        """
        Reset the sample definition. This function should be never
        called.
        """
        self.beginResetModel()
        self.__chem_symbols = sample.cell.get_chemical_symbols()
        self.__pos = sample.cell.get_scaled_positions()
        self.__mm = copy(sample.mm)
        self.__mmdata = self.__mm.fc_get(self.__coord)
        self.__mmbkp = (copy(self.__mm.fc_get(coord)),copy(coord))
        self.__changed = [[False,]*6 for _ in range(len(self.__chem_symbols))]
        self.endResetModel()
        
    def changeCoord(self, coord):
        """
        Change coordinate system for Fourier Components input.
         
        """
        if coord == self.__coord:
            return
        else:
            if (coord >=0 and coord < 3):
                self.beginResetModel()
                self.__mm.fc_set(self.__mmdata,self.__coord)
                self.__coord = coord
                self.__mmdata = self.__mm.fc_get(self.__coord)
                self.endResetModel()

    def reset(self):
        self.beginResetModel()
        self.__mm.fc_set(self.__mmbkp[0],self.__mmbkp[1])
        self.__mmdata = self.__mm.fc_get(self.__coord)
        self.__changed = [[False,]*6 for _ in range(len(self.__chem_symbols))]
        self.endResetModel()
    
    def getData(self):
        """
        Return the current (updated) magnetic model object and the 
        currently selected coordinate system.
        
        """
        return (self.__mmdata, self.__coord)

    def rowCount(self, parent=QModelIndex()):
        if parent.isValid():
            return 0
        else:
            return len(self.__chem_symbols)

    def columnCount(self, parent=QModelIndex()):
        if parent.isValid():
            return 0
        else:
            return len(self.__headers)
            

    def parent(self, index):
        return QModelIndex()

    def index(self, row, column=0, parent=QModelIndex()):
        if parent.isValid() or \
                column < 0 or column >= self.columnCount() or \
                row < 0 or row >= self.rowCount():
            return QModelIndex()

        return self.createIndex(row, column, row)

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):

        if section >= 0 and section < 8 and orientation == QtCore.Qt.Horizontal:
            if role == QtCore.Qt.DisplayRole:
                return self.__headers[section]

        return QAbstractItemModel.headerData(self, section, orientation, role)

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if self._valid(index):
            #key = self._keyFromIndex(index)
            #print(index.row(),index.column())
            column = index.column()
            row = index.row()
            if role == QtCore.Qt.DisplayRole:
                if column == 0:
                    return self.__chem_symbols[row]
                elif column == 1:
                    p = self.__pos[row].tolist()
                    return "{:3.3f} {:3.3f} {:3.3f}".format(*p)
                elif column >= 2:
                    #template string used to format Fourier Components in 
                    # the table. 4 digits for floats should be enough
                    templ_str = "{:3.4f}"
                    
                    # first put real values (for cols 2,3,4) then Imag.
                    if column <= 4:
                        return templ_str.format(self.__mmdata[row][column-2].real)
                    else:
                        return templ_str.format(self.__mmdata[row][column-5].imag)
                return self
            elif role == QtCore.Qt.FontRole and self.__changed[row][column-2]:
                font= QtGui.QFont()
                font.setBold(True)
                return font

        return None

    def flags(self, index):
        """
        Only set Fourier components fields to editable. All the rest is
        untouchable.
        """
        if self._valid(index):
            flags = QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable
            if index.column() >= 2:
                return QtCore.Qt.ItemIsEditable | flags
            else:
                return flags
        return QtCore.Qt.NoItemFlags


    def setData(self, index, value, role=QtCore.Qt.EditRole):
        if self._valid(index) and index.column() >= 2:
            
            r = index.row()
            c = index.column()-2
            
            try:
                if (c>2): #imaginary part
                    self.__mmdata[r,c-3] = self.__mmdata[r,c-3].real + 1j*float(value)
                else:     #real part
                    self.__mmdata[r,c] = float(value) + 1j*self.__mmdata[r,c].imag
                
                self.__changed[r][c] = True
                
                
                self.__mm.fc_set(self.__mmdata,self.__coord)
                
            except (TypeError, ValueError) as ex:
                log.error("Failed to set value (%r) for element (%d, %d)", value, r,c,
                          exc_info=True)
            else:
                self.dataChanged.emit(index, index)
                return True

        return False

    def _valid(self, index):
        row = index.row()
        return row >= 0 and row < self.rowCount()


class MuonsTabEntry(QtWidgets.QWidget):
    """
    This is the widget that provides the muon input functions
    """
    
    # signal that is emitted every time a new list of muons is applied
    muonAdd = QtCore.pyqtSignal([list])
    
    def __init__(self, parent):
        #self.parent=parent
        
        super(MuonsTabEntry,self).__init__(parent)
        
        self.verticalLayout = QtWidgets.QVBoxLayout(self)

        self.horizontalLayout = QtWidgets.QHBoxLayout()

        self.TxtEditMuonSites = QtWidgets.QPlainTextEdit(self)

        self.horizontalLayout.addWidget(self.TxtEditMuonSites)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()

        self.BtnFormat = QtWidgets.QPushButton(self)
        self.verticalLayout_2.addWidget(self.BtnFormat)
        
        self.BtnValidate = QtWidgets.QPushButton(self)
        self.verticalLayout_2.addWidget(self.BtnValidate)
        
        self.BtnApply = QtWidgets.QPushButton(self)
        self.verticalLayout_2.addWidget(self.BtnApply)
        
        self.BtnShow = QtWidgets.QPushButton(self)
        self.verticalLayout_2.addWidget(self.BtnShow)

        self.BtnEq = QtWidgets.QPushButton(self)
        self.verticalLayout_2.addWidget(self.BtnEq)
        
        
        
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.frame = QtWidgets.QFrame(self)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)

        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.frame)

        self.label = QtWidgets.QLabel(self.frame)

        self.horizontalLayout_2.addWidget(self.label)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.SBMuonGenGridX = QtWidgets.QSpinBox(self.frame)
        self.SBMuonGenGridX.setMinimum(1)
        self.SBMuonGenGridX.setMaximum(10000)

        self.horizontalLayout_2.addWidget(self.SBMuonGenGridX)
        self.SBMuonGenGridY = QtWidgets.QSpinBox(self.frame)
        self.SBMuonGenGridY.setMinimum(1)
        self.SBMuonGenGridY.setMaximum(10000)

        self.horizontalLayout_2.addWidget(self.SBMuonGenGridY)
        self.SBMuonGenGridZ = QtWidgets.QSpinBox(self.frame)
        self.SBMuonGenGridZ.setMinimum(1)
        self.SBMuonGenGridZ.setMaximum(10000)

        self.horizontalLayout_2.addWidget(self.SBMuonGenGridZ)
        self.BtnMuonGenGrid = QtWidgets.QPushButton(self.frame)

        self.horizontalLayout_2.addWidget(self.BtnMuonGenGrid)
        self.verticalLayout.addWidget(self.frame)

        self.BtnFormat.setText(("Format"))
        self.BtnValidate.setText(("Validate"))
        self.BtnApply.setText(("Apply"))
        self.BtnShow.setText(("Show"))
        self.BtnEq.setText(("Find Equivalent"))
        self.label.setText(("Generate Grid (insert number of points along a, b, c)"))
        self.BtnMuonGenGrid.setText(("Generate"))
        

        
        self.BtnFormat.clicked.connect(self._format)
        self.BtnValidate.clicked.connect(self._validate)
        self.BtnApply.clicked.connect(self._apply)
        self.BtnApply.setEnabled(False)
        self.BtnShow.setEnabled(False)
        self.BtnEq.setEnabled(False)
        
        self.BtnMuonGenGrid.clicked.connect(self._generate_grid)
        
        self.txt = ""
        self.validation = False
        
    
    def _generate_grid(self):
        """
        Generate a list of muon positions and possibly show them.
        """
        
        nx = float(self.SBMuonGenGridX.value())
        ny = float(self.SBMuonGenGridY.value())
        nz = float(self.SBMuonGenGridZ.value())
        
        tot = int(nx*ny*nz)
        
        
        mlst = []
        self.txt = ""
        if tot < CONST_MAX_MUON_SITES_IN_EDITOR:
            for i in range(int(nx)):
                for j in range(int(ny)):
                    for k in range(int(nz)):
                        p=[i/nx,j/ny,k/nz]
                        self.txt += '{}\t{}\t{}\n'.format(*p)
                        mlst.append(p)
            
        else:
            self.txt = "{} muon sites defined. Way too many to be shown.".format(tot)
            for i in range(int(nx)):
                  for j in range(int(ny)):
                      for k in range(int(nz)):
                          mlst.append([i/nx,j/ny,k/nz])
                        
        self.TxtEditMuonSites.setPlainText(self.txt)
        #self.parent().add_muons(mlst)
        self.muonAdd.emit(mlst)
        
    
    def _format(self):
        """
        Format a sort of table as txt tabulators. Useful to check input.
        """
        
        txt = self.TxtEditMuonSites.toPlainText()
        self.TxtEditMuonSites.clear()
        # this is rather impossible to read. Here's what it does:
        # 1) removes white lines
        # 2) spaces (" ") becomes Tabs "\t" if line is not a comment

        self.txt = "\n".join(['\t'.join(ll.strip().split()) \
                          if ll.lstrip()[0] != '#' else \
                             ll.rstrip() \
                          for ll in txt.splitlines() if ll.strip()])
            
        self.TxtEditMuonSites.setPlainText(self.txt)
        
    def _validate(self):
        """
        Checks if input data is valid.
        """
        if self.validation == False:
            self._format()
            self.BtnFormat.setEnabled(False)
            self.BtnValidate.setText("Edit")
            self.TxtEditMuonSites.setReadOnly(True)
            self.validation = True
        else:
            self.TxtEditMuonSites.clear()
            self.TxtEditMuonSites.setPlainText(self.txt)
            self.BtnFormat.setEnabled(True)
            self.BtnApply.setEnabled(False)
            self.BtnShow.setEnabled(False)
            self.TxtEditMuonSites.setReadOnly(False)
            self.BtnValidate.setText("Validate")
            self.validation = False
            return
      

        
        self.TxtEditMuonSites.clear()
        
        n_errors = 0
        
        for ll in self.txt.splitlines():
            if ll.lstrip()[0] == '#':
                self.TxtEditMuonSites.appendHtml( \
                        '<font color="orange"><pre>'+ll.strip()+'</pre></font>')
                continue
                
            try:
                if len([float(x) for x in ll.split()]) != 3:
                    self.TxtEditMuonSites.appendHtml( \
                        '<font color="red"><pre>'+ll.strip()+'</pre> #invalid number of entries</font>')
                    n_errors += 1
                    continue
                
            except ValueError:
                self.TxtEditMuonSites.appendHtml( \
                        '<font color="red"><pre>'+ll.strip()+'</pre> #invalid values</font>')
                n_errors += 1
                continue
                
            self.TxtEditMuonSites.appendHtml( \
                        '<font color="green"><pre>'+ll.strip()+'</pre></font>')
                
        if n_errors == 0:
            self.BtnApply.setEnabled(True)
                
    def _apply(self):
        """
        This function converts the (validated) text to a 2D list of
        floats that is then sent to the main Sample Object with the 
        emit command. 
        """
        v = []
        for ll in self.txt.splitlines():
            if ll.lstrip()[0] == '#':
                continue
                
            try:
                vv = [float(x) for x in ll.split()]
                if len(vv) != 3:
                    log.error("Problems parsing the muon sites")
                    break
                else:
                    v.append(vv)
                
            except:
                log.error("Problems parsing the muon sites")
                break
        
        #self.parent().add_muons(v)
        self.muonAdd.emit(v)
        
    def update_data(self, sample):
        """
        This function updates the interface when new positions are available
        (for example loaded from a file).
        """
        
        self.TxtEditMuonSites.clear()
        
        #check len before adding too many lines
        try:
            m = sample.muons
        except MuonError:
            self.TxtEditMuonSites.clear()
            self.txt = ""
            self._setEditableState()
            return
        
        if len(m) == 0:
            self.txt = ""
            return
        elif len(m) >= CONST_MAX_MUON_SITES_IN_EDITOR:
            self.txt = "{} muon sites defined. Way too many to be shown.".format(len(m))
            self.TxtEditMuonSites.setPlainText(self.txt)
            self._setEditableState()
            return
        
        self.txt = ""
        for p in sample.muons:
            self.txt += "{}\t{}\t{}\n".format(*p.tolist())
        
        self._setEditableState()
        
    def _setEditableState(self):
        self.TxtEditMuonSites.setPlainText(self.txt)
        
        self.BtnFormat.setEnabled(False)
        self.BtnValidate.setText("Edit")
        self.BtnValidate.setEnabled(True)
        self.BtnApply.setEnabled(False)
        self.BtnShow.setEnabled(True)
        
        self.TxtEditMuonSites.setReadOnly(True)
        
        self.validation = True
        
        

class CrystalTabEntry(QtWidgets.QWidget):
    """
    This is the interface that deals with crystal structures.
    
    """
    
    # signal emitted when user presses the Load Structure button
    crysLoadCell = QtCore.pyqtSignal()
    
    # signal emitted when user presses the Show Structure button
    crysShowCell = QtCore.pyqtSignal()
    
    def __init__(self, parent, *args):
        super(CrystalTabEntry,self).__init__(parent)

        
        self.horizontalLayout = QtWidgets.QHBoxLayout(self)
        
        self.InfoLabel = QtWidgets.QPlainTextEdit(self)
        self.InfoLabel.setReadOnly(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.InfoLabel.sizePolicy().hasHeightForWidth())
        self.InfoLabel.setSizePolicy(sizePolicy)
        
        self.horizontalLayout.addWidget(self.InfoLabel)
        
        self.verticalLayout = QtWidgets.QVBoxLayout()
        
        self.BtnLoadStructure = QtWidgets.QPushButton(self)

        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.BtnLoadStructure.sizePolicy().hasHeightForWidth())
        
        self.BtnLoadStructure.setSizePolicy(sizePolicy)
        
        self.verticalLayout.addWidget(self.BtnLoadStructure)
        
        self.BtnShowStructure = QtWidgets.QPushButton(self)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.BtnShowStructure.sizePolicy().hasHeightForWidth())
        self.BtnShowStructure.setSizePolicy(sizePolicy)
        
        self.verticalLayout.addWidget(self.BtnShowStructure)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.horizontalLayout.addLayout(self.verticalLayout)

        self.InfoLabel.setPlainText(("Structure not defined"))
        self.BtnLoadStructure.setText(("Load (xsf or CIF)"))
        self.BtnShowStructure.setText(("Show"))
        
        
        self.BtnLoadStructure.clicked.connect(self._load_file)
        self.BtnShowStructure.clicked.connect(self._show_cell)

    def _load_file(self):
        """
        Send request to load a new cell through by emitting the 
        crysLoadCell signal.
        """
        self.crysLoadCell.emit()
        
    def _show_cell(self):
        self.crysShowCell.emit()
    
    def update_data(self,sample):
        try:
            sample._check_lattice()
        except:
            return
        
        if (sample.cell != None):
            cell = sample.cell
            symbols = cell.get_chemical_symbols()
            lattice = cell.get_cell()
            buf = "Lattice vectors:\n"
            buf += ("  a %20.15f %20.15f %20.15f\n" % tuple( lattice[0] ))
            buf += ("  b %20.15f %20.15f %20.15f\n" % tuple( lattice[1] ))
            buf += ("  c %20.15f %20.15f %20.15f\n" % tuple( lattice[2] ))
            buf += '\n'
            buf += ("Atomic positions (fractional):\n")
        
            for i, v in enumerate(cell.get_scaled_positions()):
                buf += ("%5d %-2s%18.14f%18.14f%18.14f\n" % \
                            (i+1, symbols[i], v[0], v[1], v[2]))
        
              
            self.InfoLabel.setPlainText(buf)


class MagnetismTabEntry(QtWidgets.QWidget):
    magDescChanged = QtCore.pyqtSignal([int,str])
    magKChanged = QtCore.pyqtSignal([int,list])
    magFCChanged = QtCore.pyqtSignal([int,tuple])
    magSelectOrder = QtCore.pyqtSignal([int])
    magAddOrder = QtCore.pyqtSignal()
    
    magShowMagStructure = QtCore.pyqtSignal()
    magExportMagStructure = QtCore.pyqtSignal()
    
    def __init__(self, parent):

        self.model = None
        
        super(QtWidgets.QWidget, self).__init__(parent=parent)


        self.verticalLayout = QtWidgets.QVBoxLayout(self)

        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()

        self.MagnList = QtWidgets.QComboBox(self)

        self.horizontalLayout_3.addWidget(self.MagnList)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.BtnExportFP = QtWidgets.QPushButton(self)

        self.horizontalLayout_3.addWidget(self.BtnExportFP)
        self.BtnShow = QtWidgets.QPushButton(self)

        self.horizontalLayout_3.addWidget(self.BtnShow)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.FrameMM = QtWidgets.QFrame(self)
        self.FrameMM.setEnabled(False)
        self.FrameMM.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.FrameMM.setFrameShadow(QtWidgets.QFrame.Raised)

        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.FrameMM)

        self.horizontalLayout = QtWidgets.QHBoxLayout()

        self.label_3 = QtWidgets.QLabel(self.FrameMM)

        self.horizontalLayout.addWidget(self.label_3)
        self.TxtDesc = QtWidgets.QLineEdit(self.FrameMM)

        self.horizontalLayout.addWidget(self.TxtDesc)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.label = QtWidgets.QLabel(self.FrameMM)

        self.horizontalLayout.addWidget(self.label)
        self.TxtKx = QtWidgets.QLineEdit(self.FrameMM)

        self.horizontalLayout.addWidget(self.TxtKx)
        self.TxtKy = QtWidgets.QLineEdit(self.FrameMM)

        self.horizontalLayout.addWidget(self.TxtKy)
        self.TxtKz = QtWidgets.QLineEdit(self.FrameMM)

        self.horizontalLayout.addWidget(self.TxtKz)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.line = QtWidgets.QFrame(self.FrameMM)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)

        self.verticalLayout_2.addWidget(self.line)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()

        self.label_2 = QtWidgets.QLabel(self.FrameMM)

        self.horizontalLayout_2.addWidget(self.label_2)
        self.CmbCoordSystem = QtWidgets.QComboBox(self.FrameMM)

        self.horizontalLayout_2.addWidget(self.CmbCoordSystem)
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem2)
        self.label_4 = QtWidgets.QLabel(self.FrameMM)

        self.horizontalLayout_2.addWidget(self.label_4)
        self.BtnResetFC = QtWidgets.QPushButton(self.FrameMM)

        self.horizontalLayout_2.addWidget(self.BtnResetFC)
        self.BtnApplyFC = QtWidgets.QPushButton(self.FrameMM)

        self.horizontalLayout_2.addWidget(self.BtnApplyFC)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.FCTableView = QtWidgets.QTableView(self.FrameMM)

        self.verticalLayout_2.addWidget(self.FCTableView)
        self.verticalLayout.addWidget(self.FrameMM)

        self.BtnExportFP.setText(( "Export FP Studio"))
        self.BtnShow.setText(( "Show"))
        self.label_3.setText(( "Description"))
        self.label.setText(( "Propagation Vector (r.l.u.)"))
        self.label_2.setText(( "Coordinate system:"))
        self.label_4.setText(( "<font color=\'red\'> Remember to apply changes to Fourier Components!</font>"))
        self.BtnResetFC.setText(( "Reset"))
        self.BtnApplyFC.setText(( "Apply"))










        
        
        
        self.MagnList.insertItems(0,['Add new'])
        self.MagnList.activated.connect(self._select_mag_order)
        
        self.BtnShow.clicked.connect(self._show_mag_order)
        self.BtnExportFP.clicked.connect(self._export_mag_order)
        
        self.CmbCoordSystem.insertItems(0,['Cartesian'])
        self.CmbCoordSystem.insertItems(1,['uB/Ang (x||a, y||b, z||c)'])
        self.CmbCoordSystem.insertItems(2,['uB (x||a, y||b, z||c)'])
        
        self.CmbCoordSystem.activated.connect(self._update_coordinates)
        
        self.TxtDesc.textChanged.connect(self._update_description)
        
        
        validator = QtGui.QRegExpValidator(QtCore.QRegExp("[-+]?[0-9]*\.?[0-9]+"))
        
        self.TxtKx.setValidator(validator)
        self.TxtKx.setText("0")
        
        self.TxtKy.setValidator(validator)
        self.TxtKy.setText("0")
        
        self.TxtKz.setValidator(validator)
        self.TxtKz.setText("0")

        self.TxtKx.editingFinished.connect(self._apply_K)
        self.TxtKy.editingFinished.connect(self._apply_K)
        self.TxtKz.editingFinished.connect(self._apply_K)
        
        
        self.BtnResetFC.clicked.connect(self._reset_FC)
        self.BtnApplyFC.clicked.connect(self._apply_FC)
    
    def _show_mag_order(self):
        idx =  self._get_mm_idx_from_combo()
        
        if idx >= 0:
            #self.parent().show_magnetic_supercell()
            self.magShowMagStructure.emit()
        else:
            log.error("No magnetic order selected")
            
    def _export_mag_order(self):
        idx =  self._get_mm_idx_from_combo()
        
        if idx >= 0:
            #self.parent().export_magneitc_order()
            self.magExportMagStructure.emit()
        else:
            log.error("No magnetic order selected")
    
    def _get_mm_idx_from_combo(self):
        return self.MagnList.currentIndex()-1
        
    def _update_description(self, text):
        self.MagnList.setItemText(self.MagnList.currentIndex(),text)
        if self.model:
            self.magDescChanged.emit(self._get_mm_idx_from_combo(),text)
        #self.parentWidget().update_mm_desc(self._get_mm_idx_from_combo(),text)
        
    def _apply_K(self):
        """
        Parse propagation vector and set it by emitting the corresponding
        signal.
        """
        
        v = [0.,0.,0.]
        
        if self.TxtKx.text() == "":
            self.TxtKx.setText("0")
        if self.TxtKy.text() == "":
            self.TxtKy.setText("0")
        if self.TxtKz.text() == "":
            self.TxtKz.setText("0")
        
        try:
            v[0] = float(self.TxtKx.text())
            v[1] = float(self.TxtKy.text())
            v[2] = float(self.TxtKz.text())
        except:
            log.error("Invalid value for propagation vector.")
        else:
            #self.parent().update_mm_K(self._get_mm_idx_from_combo(),v)
            self.magKChanged.emit(self._get_mm_idx_from_combo(),v)
    
    def _reset_FC(self):
        if self.model:
            self.model.reset()
            
    def _apply_FC(self):
        if self.model:
            #self.parent().update_mm_FC(self._get_mm_idx_from_combo(), self.model.getData())
            self.magFCChanged.emit(self._get_mm_idx_from_combo(),self.model.getData())
    
    def _select_mag_order(self, index):
        if (index == 0):
            #self.parent().add_mag_order()
            self.magAddOrder.emit()
        else:
            #self.parent().select_mag_order(index-1)
            if self.model:
                self.magSelectOrder.emit(index-1)
            
    def _update_coordinates(self, index):
        if (self.model):
            self.model.changeCoord(index)
        
    def update_data(self,sample):
        i = 0
        selected = -1
        if (sample.mm_count > 0):
            selected = sample.current_mm_idx
        else:
            # reset everything
            self.model = None
            self.FCTableView.setModel(QtGui.QStandardItemModel())
            self.MagnList.clear()
            self.MagnList.insertItems(0,['Add new'])
            self.TxtDesc.setText("")
            
            
        
        for i in range(sample.mm_count):
            sample.current_mm_idx=i
            if (self.MagnList.count() <= i+1):
                self.MagnList.insertItems(sample.current_mm_idx+1,[sample.mm.desc])
        
        if selected >= 0:
            self.FrameMM.setEnabled(True)
            sample.current_mm_idx=selected
            
            self.MagnList.setCurrentIndex(sample.current_mm_idx+1)
            self.TxtDesc.setText(sample.mm.desc)
        
            coord = self.CmbCoordSystem.currentIndex()
            self.model = MagnModel(sample = sample, coord=coord)
        
            self.FCTableView.setModel(self.model)
        
            self.TxtKx.setText(str(sample.mm.k[0]))
        
            self.TxtKy.setText(str(sample.mm.k[1]))
            
            self.TxtKz.setText(str(sample.mm.k[2]))
        
class CalcTabEntry(QtWidgets.QWidget):
    """
    This widget defines all the details of the calculation.
    """
    
    
    calcRequestLorSphere = QtCore.pyqtSignal([list])
    calcQueueSum = QtCore.pyqtSignal([dict])
    
    def __init__(self, parent, *args):
        
        super(CalcTabEntry,self).__init__(parent, *args)
        
        self.verticalLayout = QtWidgets.QVBoxLayout(self)

        self.groupBox = QtWidgets.QGroupBox(self)

        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox)

        self.horizontalLayout = QtWidgets.QHBoxLayout()

        self.label = QtWidgets.QLabel(self.groupBox)

        self.horizontalLayout.addWidget(self.label)
        
        spacerItem = QtWidgets.QSpacerItem(497, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        
        self.SCa = QtWidgets.QSpinBox(self.groupBox)
        self.SCa.setMinimum(1)
        self.SCa.setMaximum(1000)
        self.SCa.setSingleStep(2)

        self.horizontalLayout.addWidget(self.SCa)
        
        self.SCb = QtWidgets.QSpinBox(self.groupBox)
        self.SCb.setMinimum(1)
        self.SCb.setMaximum(1000)
        self.SCb.setSingleStep(2)

        self.horizontalLayout.addWidget(self.SCb)
        
        self.SCc = QtWidgets.QSpinBox(self.groupBox)
        self.SCc.setMinimum(1)
        self.SCc.setMaximum(1000)
        self.SCc.setSingleStep(2)

        self.horizontalLayout.addWidget(self.SCc)
        
        self.verticalLayout_3.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.label_2 = QtWidgets.QLabel(self.groupBox)

        self.horizontalLayout_2.addWidget(self.label_2)
        
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        
        self.LELorentzSphere = QtWidgets.QLineEdit(self.groupBox)

        self.horizontalLayout_2.addWidget(self.LELorentzSphere)
        self.BtnFindLorSphere = QtWidgets.QPushButton(self.groupBox)

        self.horizontalLayout_2.addWidget(self.BtnFindLorSphere)
        self.verticalLayout_3.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()

        self.label_6 = QtWidgets.QLabel(self.groupBox)

        self.horizontalLayout_7.addWidget(self.label_6)
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_7.addItem(spacerItem2)
        self.SBContNNN = QtWidgets.QSpinBox(self.groupBox)
        self.SBContNNN.setMinimum(1)
        self.SBContNNN.setMaximum(20)
        self.SBContNNN.setSingleStep(1)

        self.horizontalLayout_7.addWidget(self.SBContNNN)
        self.verticalLayout_3.addLayout(self.horizontalLayout_7)
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout()

        self.label_7 = QtWidgets.QLabel(self.groupBox)

        self.horizontalLayout_8.addWidget(self.label_7)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem3)
        self.SBContRad = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.SBContRad.setMinimum(0.0)
        self.SBContRad.setSingleStep(0.5)
        self.SBContRad.setProperty("value", 3.0)

        self.horizontalLayout_8.addWidget(self.SBContRad)
        self.verticalLayout_3.addLayout(self.horizontalLayout_8)
        self.verticalLayout.addWidget(self.groupBox)
        self.groupBox_2 = QtWidgets.QGroupBox(self)

        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_2)

        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()

        self.label_3 = QtWidgets.QLabel(self.groupBox_2)

        self.horizontalLayout_3.addWidget(self.label_3)
        spacerItem4 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem4)
        self.BtnCalcSimple = QtWidgets.QRadioButton(self.groupBox_2)
        self.BtnCalcSimple.setChecked(True)

        self.BtnGCalcType = QtWidgets.QButtonGroup(self)

        self.BtnGCalcType.addButton(self.BtnCalcSimple)
        self.horizontalLayout_3.addWidget(self.BtnCalcSimple)
        self.BtnCalcRotate = QtWidgets.QRadioButton(self.groupBox_2)

        self.BtnGCalcType.addButton(self.BtnCalcRotate)
        self.horizontalLayout_3.addWidget(self.BtnCalcRotate)
        self.BtnCalcIncom = QtWidgets.QRadioButton(self.groupBox_2)

        self.BtnGCalcType.addButton(self.BtnCalcIncom)
        self.horizontalLayout_3.addWidget(self.BtnCalcIncom)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()

        self.label_4 = QtWidgets.QLabel(self.groupBox_2)

        self.horizontalLayout_4.addWidget(self.label_4)
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem5)
        self.SBnAngles = QtWidgets.QSpinBox(self.groupBox_2)
        self.SBnAngles.setMinimum(1)
        self.SBnAngles.setMaximum(10000)

        self.horizontalLayout_4.addWidget(self.SBnAngles)
        self.verticalLayout_2.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()

        self.label_5 = QtWidgets.QLabel(self.groupBox_2)

        self.horizontalLayout_5.addWidget(self.label_5)
        spacerItem6 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_5.addItem(spacerItem6)
        self.LEAxisX = QtWidgets.QLineEdit(self.groupBox_2)
        self.LEAxisX.setMaximumSize(QtCore.QSize(100, 16777215))

        self.horizontalLayout_5.addWidget(self.LEAxisX)
        self.LEAxisY = QtWidgets.QLineEdit(self.groupBox_2)
        self.LEAxisY.setMaximumSize(QtCore.QSize(100, 16777215))

        self.horizontalLayout_5.addWidget(self.LEAxisY)
        self.LEAxisZ = QtWidgets.QLineEdit(self.groupBox_2)
        self.LEAxisZ.setMaximumSize(QtCore.QSize(100, 16777215))

        self.horizontalLayout_5.addWidget(self.LEAxisZ)
        self.verticalLayout_2.addLayout(self.horizontalLayout_5)
        self.verticalLayout.addWidget(self.groupBox_2)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()

        spacerItem7 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_6.addItem(spacerItem7)
        self.BtnQueue = QtWidgets.QPushButton(self)

        self.horizontalLayout_6.addWidget(self.BtnQueue)
        self.verticalLayout.addLayout(self.horizontalLayout_6)
        
        
        
        
        # set text
        
        self.groupBox.setTitle(("Supercell definition"))
        self.label.setText(("Supercell size (repetitions along a,b,c)"))
        self.label_2.setText(("Lorentz Radius (Ang.)"))
        self.BtnFindLorSphere.setText(("Suggest Optimal"))
        self.label_6.setText(("Number of Nearest Neighbors for Contact Term"))
        self.label_7.setText(("Maximum muon-atom distance for contact term (Ang.)"))
        self.groupBox_2.setTitle(("Calculation"))
        self.label_3.setText(("Type"))
        self.BtnCalcSimple.setText(("Simple"))
        self.BtnCalcRotate.setText(("Rotation"))
        self.BtnCalcIncom.setText(("Incommensurate"))
        self.label_4.setText(("Angles"))
        self.label_5.setText(("Axis"))
        self.BtnQueue.setText(("Queue Sum"))
        
        
        # connect signals
        parent.newLorentzSphere.connect(self._updateLorentzSphere)
        
        
        
        #self.parent = parent
        
        self.BtnQueue.clicked.connect(self._queueCalc)
        self.BtnFindLorSphere.clicked.connect(self._findLargestSphere)
        
        validator = QtGui.QRegExpValidator(QtCore.QRegExp("[-+]?[0-9]*\.?[0-9]+"),self.LELorentzSphere)
        self.LELorentzSphere.setValidator(validator)
        
        self.LELorentzSphere.setText("0")
        
        self.LEAxisX.setText("0")
        self.LEAxisX.setValidator(validator)
        self.LEAxisY.setText("0")
        self.LEAxisY.setValidator(validator)
        self.LEAxisZ.setText("0")
        self.LEAxisZ.setValidator(validator)
        
        
        self.BtnGCalcType.setId(self.BtnCalcSimple,0)
        self.BtnGCalcType.setId(self.BtnCalcRotate,1)
        self.BtnGCalcType.setId(self.BtnCalcIncom, 2)
    
    def _findLargestSphere(self):
        
        sc_a = int(self.SCa.value())
        sc_b = int(self.SCb.value())
        sc_c = int(self.SCc.value())
      
        #v = self.parent().getLargestLorentzSphere([sc_a,sc_b,sc_c])
        self.calcRequestLorSphere.emit([sc_a,sc_b,sc_c])
        
    def _updateLorentzSphere(self, v):
        self.LELorentzSphere.setText("{:03.4f}".format(v))
    
    def _queueCalc(self):
        sc_a = int(self.SCa.value())
        sc_b = int(self.SCb.value())
        sc_c = int(self.SCc.value())
        
        c_type_str = ''
        c_type_int = self.BtnGCalcType.checkedId()
        
        if c_type_int == 0:
            c_type_str = 's'
        elif c_type_int == 1:
            c_type_str = 'r'
        elif c_type_int == 2:
            c_type_str = 'i'
        else:
            log.error("Invalid calculation type")
            return
        
        r_float = 0.
        try:
            r_float = float(self.LELorentzSphere.text())
        except:
            log.error("Invalid Lorentz Sphere value")
            return

        i_nangles = None
        try:
            i_nangles = int(self.SBnAngles.value())
            assert(i_nangles>=1)
        except:
            log.error("Invalid numer of angles")
            return
        
        i_nnn = 2
        try:
            i_nnn = int(self.SBContNNN.value())
            assert(i_nnn >= 0)
        except:
            log.error("Invalid numer of nearest neighbors for Contact term")
            return
            
        r_cont_float = 10.0
        try:
            r_cont_float = float(self.SBContRad.value())
        except:
            log.error("Invalid value for distance between muon and atoms for Contact term")
            return
            
        if c_type_int == 0:
            data = dict(ctype=c_type_str,supercellsize=[sc_a,sc_b,sc_c],
                    radius=r_float, nnn=i_nnn, rcont=r_cont_float)
        elif c_type_int == 1:
            lst_axis = [0,0,0]
            try:
                lst_axis[0] = float(self.LEAxisX.text())
                lst_axis[1] = float(self.LEAxisY.text())
                lst_axis[2] = float(self.LEAxisZ.text())
                assert(lst_axis[0]**2+lst_axis[1]**2+lst_axis[2]**2 > 0)
            except (AssertionError, ValueError):
                log.error("Invalid rotation axis")
                return
            data = dict(ctype=c_type_str,supercellsize=[sc_a,sc_b,sc_c],
                    radius=r_float, nnn=i_nnn, rcont=r_cont_float,
                    nangles=i_nangles,axis=lst_axis)
            
        elif c_type_int == 2:
            data = dict(ctype=c_type_str,supercellsize=[sc_a,sc_b,sc_c],
                    nnn=i_nnn, rcont=r_cont_float,
                    radius=r_float, nangles=i_nangles)
        
        # emit signal to request calculation
        self.calcQueueSum.emit(data)

class CalcThread(threading.Thread):
    """
    Class defining the htread that suns the simulations.
    """
  
    def __init__(self, in_q, out_q):
        super(CalcThread, self).__init__()
        #get sample
        self.in_q = in_q
        self.out_q = out_q
        
    def run(self):
        """
        This function actually runs the simulations. Each sample has its
        own thread to run the simulations. This means that all simulations
        for a given sample enter into a FIFO queue. This is done to avoid
        scrambling the order of the simulations.
        Calculations can be done in parallel if they belong to different
        samples.
        
        The function trys to get a task 3 times with a timeout of 1s then 
        exit.
        """

        for i in range(3): # timeout is 1s so wait 3s and leave function if there is nothing to do
            try:
                sample, calc_data = self.in_q.get(True,1)
            except Empty:
                continue
            
            try:
                r = locfield(sample,**calc_data)
                self.out_q.put(r)
            except:
                print("something very bad happened")



class TotalFieldListItem(QtGui.QStandardItem):
    """
    Simple class overloading QtGui.QStandardItem, just for style.
    """
    def __init__(self, LFobj, index = -1):
        
        self._lf = LFobj
        self._idx = index
        
        if (self._idx <0):
            super(TotalFieldListItem,self).__init__('{:03.4f}, {:03.4f}, {:03.4f}'.format(*self._lf.T.tolist()))
        else:
            super(TotalFieldListItem,self).__init__('{:03.4f}, {:03.4f}, {:03.4f}'.format(*self._lf.T[self._idx].tolist()))
        
        self.setEditable(False)
        
class LorentzFieldListItem(QtGui.QStandardItem):
    """
    Simple class overloading QtGui.QStandardItem, just for style.
    """
    def __init__(self, LFobj, index = -1):
        
        self._lf = LFobj
        self._idx = index
        
        if (self._idx <0):
            super(LorentzFieldListItem,self).__init__('{:03.4f}, {:03.4f}, {:03.4f}'.format(*self._lf.L.tolist()))
        else:
            super(LorentzFieldListItem,self).__init__('{:03.4f}, {:03.4f}, {:03.4f}'.format(*self._lf.L[self._idx].tolist()))
        
        self.setEditable(False)
        
class DipolarFieldListItem(QtGui.QStandardItem):
    """
    Simple class overloading QtGui.QStandardItem, just for style.
    """
    def __init__(self, LFobj, index = -1):
        
        self._lf = LFobj
        self._idx = index
        
        if (self._idx <0):
            super(DipolarFieldListItem,self).__init__('{:03.4f}, {:03.4f}, {:03.4f}'.format(*self._lf.D.tolist()))
        else:
            super(DipolarFieldListItem,self).__init__('{:03.4f}, {:03.4f}, {:03.4f}'.format(*self._lf.D[self._idx].tolist()))
        
        self.setEditable(False)
        
        
class ContactFieldListItem(QtGui.QStandardItem):
    """
    This custom QStandardItem allows to update the value of the contact term
    every time a new contact coupling is inserted.
    """
    def __init__(self, LFobj, index = -1):
        
        self._lf = LFobj
        self._idx = index
        
        if (self._idx <0):
            super(ContactFieldListItem,self).__init__('{:03.4f}, {:03.4f}, {:03.4f}'.format(*self._lf.C.tolist()))
        else:
            super(ContactFieldListItem,self).__init__('{:03.4f}, {:03.4f}, {:03.4f}'.format(*self._lf.C[self._idx].tolist()))
        
        self.setEditable(False)
        
    def data(self, role):
        if role == QtCore.Qt.DisplayRole:
            if (self._idx < 0):
                return '{:03.4f}, {:03.4f}, {:03.4f}'.format(*self._lf.C.tolist())
            else:
                return '{:03.4f}, {:03.4f}, {:03.4f}'.format(*self._lf.C[self._idx].tolist())
        else:
            return super(ContactFieldListItem,self).data(role)



# overloaded items for list
class ContactCouplingListItem(QtGui.QStandardItem):
    """
    This class allows to validate the input of the contact coupling term.
    """
    def __init__(self, *args):
        super(ContactCouplingListItem,self).__init__(*args)
        
    def setData(self,value,role):

        if role == QtCore.Qt.EditRole:
            try:
                v = float(value)
            except:
                log.error("Invalid value for Contact coupling term.")
                return
        #if everything fine finally set
        return super(ContactCouplingListItem,self).setData(value,role)



class LFObjListItem(QtGui.QStandardItem):
    """
    This class defines the item which is inserted in the results tree view
    after a succesful simulation.
    It behaves differently for 1D (no angles involved) and 2D (many angles
    for the same muon site) results.
    
    """
    def __init__(self, LFObjs, *args):
        super(LFObjListItem,self).__init__(*args)
        self._LFObjs = LFObjs
        
        
        for i, lf in enumerate(self._LFObjs):
            muidx = QtGui.QStandardItem('mu site {}'.format(i))
            muidx.setEditable(False)
            nlf = QtGui.QStandardItem('{}'.format(1 if lf.T.ndim == 1 else lf.T.shape[0]))
            nlf.setEditable(False)
            
            acont = ContactCouplingListItem('{}'.format(lf.ACont))
            
            if (lf.T.ndim == 1): # 1D array
                L = LorentzFieldListItem(lf)
                C = ContactFieldListItem(lf)
                D = DipolarFieldListItem(lf)
                T = TotalFieldListItem(lf)
            elif (lf.T.ndim == 2): #2D array
                if lf.T.shape[0] < 20:
                    L = QtGui.QStandardItem("List")
                    C = QtGui.QStandardItem("List")
                    D = QtGui.QStandardItem("List")
                    T = QtGui.QStandardItem("List")
                    L.setEditable(False)
                    C.setEditable(False)
                    D.setEditable(False)
                    T.setEditable(False)
                    
                    for l in range(lf.T.shape[0]):
                        
                        items = [QtGui.QStandardItem(),
                                 QtGui.QStandardItem(),
                                 QtGui.QStandardItem(),
                                 LorentzFieldListItem(lf,l),
                                 DipolarFieldListItem(lf,l),
                                 ContactFieldListItem(lf,l),
                                 TotalFieldListItem(lf,l)
                                 ]
                        for item in items:
                            item.setEditable(False)
                            
                        muidx.appendRow(items)
                        
                                          
                else:
                    L = QtGui.QStandardItem("Too many items")
                    C = QtGui.QStandardItem("Too many items")
                    D = QtGui.QStandardItem("Too many items")
                    T = QtGui.QStandardItem("Too many items")
                    
                    
            else:
                log.debug("Something very odd!")
            
            
            self.appendRow([muidx, nlf, acont, L, D, C, T])
    
    
    def getLFObj(self):
        return self._LFObjs
        
    def setACont(self, muon_idx, new_value):
        if muon_idx < len(self._LFObjs):
            self._LFObjs[muon_idx].ACont = new_value
            self.emitDataChanged()
            
    
class ResultsTreeView(QtWidgets.QTreeView):
    """
    This class overloads QTreeView to provide a export data signal
    for mantid (all this is not very clean...).
    
    """
    mantidExportResult = QtCore.pyqtSignal([tuple])
    
    def __init__(self, parent, *args):
        super(ResultsTreeView,self).__init__(parent=parent, *args)
        
        #self.setSelectionBehavior(QtCore.QAbstractItemView.SelectRows)
        self._model = QtGui.QStandardItemModel()
        self._model.setHorizontalHeaderLabels(['Muon Site', 'Num loc fields', \
                                              'ACont','Lorentz (x, y, z)', \
                                              'Dipolar (x, y, z)', \
                                              'Contact (x, y, z)', \
                                              'Total (x, y, z)'])
                                              
                                              
        self._model.dataChanged.connect(self._update_LFCs)
        
        self.setModel(self._model)
        
        self.setUniformRowHeights(True) 
        
        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._openMenu)
        
        
        self.res_counter = 0 #unique counter
    
    def _update_LFCs(self, idx1, idx2):
        #this is not so nice...better logic would be great!

        item = self._model.itemFromIndex(idx1)
        
        if idx1.column() != 2:
            return
        
        obj = None
        index = idx1
        while True:
            obj = self._model.itemFromIndex(index)
            if isinstance(obj, LFObjListItem):
                break
            if index.parent().isValid():
                index = index.parent()
            else:
                break


        obj.setACont(idx1.row(),float(item.text()))
        
    def _openMenu(self, position):
        
        indexes = self.selectedIndexes()
        
        level = 0

        if len(indexes) > 0:
            index = indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1
        else:
            level = -1
        
        menu = QtWidgets.QMenu()
        if level == 0:
            if running_inside_mantid:
                a = menu.addAction(("Export to Mantid"))
                a.triggered.connect(self._export_mantid)
                
            a = menu.addAction(("Export data"))
            a.triggered.connect(self._export_data)
            a = menu.addAction(("View data"))
            a.triggered.connect(self._view_data)
        elif level >= 1:
            a = menu.addAction(("View data"))
            a.triggered.connect(self._view_data)
        else:
            a = menu.addAction(("Invalid selection"))
            a.setEnabled(False)
            
        menu.exec_(self.viewport().mapToGlobal(position))

    def _view_data(self):
        cur_sel = self.selectedIndexes()
        if len(cur_sel) > 0:
            
            index = cur_sel[0]
            LFObjElem = None

            while True:
                obj = self._model.itemFromIndex(index)
                if isinstance(obj, LFObjListItem):
                    LFObjElem = obj
                    break
                if index.parent().isValid():
                    index = index.parent()
                else:
                    break
                
            
            if LFObjElem==None:
                log.debug("Cannot get results object")
                return
            else:
                Obj = LFObjElem.getLFObj()
                self._open_view_dialog(Obj)
                
    def _export_mantid(self):
        if not running_inside_mantid:
            log.debug("Callong export_mantid outside mantid or with no queue!")
            return
      
        # get selected simulation
        index = self.selectedIndexes()
        
        try:
            obj = self._model.itemFromIndex(index[0])
            if isinstance(obj, LFObjListItem):
                self.mantidExportResult.emit((obj.text(), obj.getLFObj()))
                log.info("Results exported to mantid")
        except:
            log.debug("Something bad happened in _export_mantid")
              
    def _export_data(self):
      
        
        index = self.selectedIndexes()
        
        
        try:
            obj = self._model.itemFromIndex(index[0])
            if isinstance(obj, LFObjListItem):
                filename = select_file_dialog(self, name = "Export data",action='save', ffilter="*.dat",initdir=obj.text()+".dat",check_write=True)
                #filename = QtWidgets.QFileDialog.getSaveFileName(self, 
                #            "Export data", obj.text()+".dat", "*.dat")
                
                # cancel was clicked
                if not filename:
                    return
                
              
                with open(filename,'wb') as f:
                    for i, lf in enumerate(obj.getLFObj()):
                        f.write("# Muon site {}\n".format(i).encode('ascii'))
                        f.write("# Lorentz Dipolar Contact Total\n".format(i).encode('ascii'))
                        a = np.hstack([lf.L, lf.D,lf.C,lf.T])
                        if a.ndim == 1:
                            a = a.reshape((1,)+a.shape)
                        np.savetxt(f,a,delimiter=' ',fmt='%f')
                        
                log.info("Data saved")
            else:
                return
        except ValueError:
            log.error("Could not save")
            
    def _open_view_dialog(self, LFObj):
        d = DataPlotDialog()
        d.set_LFs(LFObj)
        d.exec_()
    
    def add_result(self, LFObj):
        self.res_counter += 1
        parent = LFObjListItem(LFObj,"Simulation {}".format(self.res_counter))
        self._model.appendRow(parent)
        # span container columns
        #self.setFirstColumnSpanned(i, view.rootIndex(), True)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # expand third container
        #index = model.indexFromItem(parent)
        #self.expand(index)        
        
        
class SampleTab(QtWidgets.QWidget):
    
    sampleNameChanged = QtCore.pyqtSignal([str])
    threadRunning = QtCore.pyqtSignal([bool])
    newLorentzSphere = QtCore.pyqtSignal([float])
    
    def __init__(self, parent = None, sample=None):
        super(SampleTab,self).__init__(parent=parent)
        
        self.gridLayout = QtWidgets.QVBoxLayout(self)
        
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        
        self.LblSampleDesc = QtWidgets.QLabel(self)
        
        
        self.horizontalLayout.addWidget(self.LblSampleDesc)
        self.LESampleDesc = QtWidgets.QLineEdit(self)
        
        
        self.horizontalLayout.addWidget(self.LESampleDesc)
        
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        
        self.BtnCreate = QtWidgets.QPushButton(self)
        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.BtnCreate.sizePolicy().hasHeightForWidth())
        self.BtnCreate.setSizePolicy(sizePolicy)
        
        self.horizontalLayout.addWidget(self.BtnCreate)
        
        self.BtnLoad = QtWidgets.QPushButton(self)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.BtnLoad.sizePolicy().hasHeightForWidth())
        
        self.BtnLoad.setSizePolicy(sizePolicy)

        self.horizontalLayout.addWidget(self.BtnLoad)
        
        
        self.BtnSave = QtWidgets.QPushButton(self)
        self.BtnSave.setEnabled(False)
        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.BtnSave.sizePolicy().hasHeightForWidth())
        self.BtnSave.setSizePolicy(sizePolicy)
        
        
        self.horizontalLayout.addWidget(self.BtnSave)
        
        
        self.gridLayout.addLayout(self.horizontalLayout)
        self.VLSampleInfoAndResults = QtWidgets.QVBoxLayout()

        self.SampleTab = QtWidgets.QTabWidget(self)
        self.SampleTab.setEnabled(True)

        self.VLSampleInfoAndResults.addWidget(self.SampleTab)
        self.LblResults = QtWidgets.QLabel(self)

        self.VLSampleInfoAndResults.addWidget(self.LblResults)
        
        self.results = ResultsTreeView(self)
        self.VLSampleInfoAndResults.addWidget(self.results)
        
        self.gridLayout.addLayout(self.VLSampleInfoAndResults)

        self.SampleTab.setCurrentIndex(-1)        
        
        self.LblSampleDesc.setText("Sample description:")
        self.BtnCreate.setText( "Create Sample")
        self.BtnLoad.setText( "Load")
        self.BtnSave.setText( "Save")
        self.LblResults.setText( "Results")
        
        #connection and final initialization
        
        
        
        self.BtnLoad.clicked.connect(self._load_sample)
        self.BtnCreate.clicked.connect(self._new_empty_sample)
        self.BtnSave.clicked.connect(self._save_sample)
        
        self.SampleTab.setEnabled(False)
        self.LESampleDesc.setEnabled(False)
        self.LESampleDesc.textChanged.connect(self._set_sample_name)
        
        
        self.crystal_interface = CrystalTabEntry(self)
        # signal/slot connection
        self.crystal_interface.crysLoadCell.connect(self.load_structure_file)
        self.crystal_interface.crysShowCell.connect(self.show_cell)
        self.SampleTab.addTab(self.crystal_interface ,"Crystal Structure")
        
        self.magnetism_interface = MagnetismTabEntry(self)
        # signal/slot connection
        self.magnetism_interface.magShowMagStructure.connect(self.show_magnetic_supercell)
        self.magnetism_interface.magExportMagStructure.connect(self.export_magneitc_order)
        self.magnetism_interface.magDescChanged.connect(self.update_mm_desc)
        self.magnetism_interface.magFCChanged.connect(self.update_mm_FC)
        self.magnetism_interface.magKChanged.connect(self.update_mm_K)
        self.magnetism_interface.magSelectOrder.connect(self.select_mag_order)
        self.magnetism_interface.magAddOrder.connect(self.add_mag_order)
        
        self.SampleTab.addTab(self.magnetism_interface,"Magnetic Order(s)")
        # signal/slot connection
        self.muon_interface = MuonsTabEntry(self)
        self.muon_interface.muonAdd.connect(self.add_muons)
        self.SampleTab.addTab(self.muon_interface,"Muon Site(s)")
        
        self.calc_interface = CalcTabEntry(self)
        # signal/slot connection
        self.calc_interface.calcRequestLorSphere.connect(self.getLargestLorentzSphere)
        self.calc_interface.calcQueueSum.connect(self.queue_calculation)
        # add to tabs
        self.SampleTab.addTab(self.calc_interface,"Calculate")
        
        #queues to comounicate with the thread that is actually executing 
        # the hard work.
        self.calc_queue = Queue() # send calculation requests
        self.res_queue = Queue()  # receives the results
        
        # This timer collects the results and joina threads when they are
        # dead.
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._parse_results)
        self.timer.start(100)
        
        self.calc_thread = None
        
        
        if sample == None:
            self.name = "No name"
            self.sample = None
        else:
            self.sample = sample
            self._update_interface()
    
    
    def _set_sample_name(self, txt):
        if self.sample:
            self.sample.name = str(txt)
            self.name = str(txt)
            self.sampleNameChanged.emit(str(txt))
    
    def _activate_buttons(self):
        self.SampleTab.setEnabled(True)
        self.LESampleDesc.setEnabled(True)
        self.BtnCreate.setEnabled(False)
        self.BtnLoad.setEnabled(False)
        self.BtnSave.setEnabled(True)
        
        
    def _parse_results(self):
        """
        This function collects the results and join the threads when
        they finish the calculation
        """
        if self.calc_thread != None:
            if (not self.calc_thread.isAlive()):
                self.calc_thread.join()
                self.calc_thread = None
                self.threadRunning.emit(False)
        else:
            # is there something we should do in this case?
            pass
          
          
            
        
        if self.res_queue.empty():
            return
        try:
            r = self.res_queue.get()
            self.results.add_result(r)
            log.info("New results ready.")
        except:
            # other error possible?
            pass
                
        
    def show_cell(self):
        
        if self.sample:
            try:
                show_cell(self.sample)
            except CellError:
                log.error("Lattice structure not defined!")
                
            
    def show_magnetic_supercell(self):
        if self.sample:
            try:
                k = self.sample.mm.k + 1e-5 # trick to remove inf
                if (k > 1/CONST_MAX_SC_CELL_FOR_XCRYSDEN).any():
                    sc = np.ceil(1/k)
                    sc[np.where(sc>CONST_MAX_SC_CELL_FOR_XCRYSDEN)]=1
                    
                elif np.linalg.norm(self.sample.mm.k) < 1e-4: #this is a FM
                    sc = [1,1,1]
                else:
                    # this should never happen
                    sc = [1,1,1]
                    log.info("Supercell would be too large! Showing just unit cell.")
                    
                
                show_supercell(self.sample, sc, False)
            except CellError:
                log.error("Lattice structure not defined!")
            except MagDefError:
                log.error("Magnetic order structure not defined!")
            except:
                log.error("Runtime error.")
                log.debug("Could show supercell", exc_info=True)
            
        
    def export_magneitc_order(self):
        if self.sample:
            try:
                fname = select_file_dialog(self, action='save',check_write=True,ffilter="*.fst")
                if fname:
                    export_fpstudio(self.sample,fname)
                
            except CellError:
                log.error("Lattice structure not defined!")
            except MagDefError:
                log.error("Magnetic order structure not defined!")
            except:
                log.error("Runtime error.")
                log.debug("Could show supercell", exc_info=True)
        
    def check(self, what=None):
        
        
        if what == 'run_sum':
            
            if not self.sample.check_status(cell=True, 
                                            magdefs=True, 
                                            muon=True, 
                                            sym=False):
                if not self.sample.check_status(cell=True):
                    log.error("Cannot run simulation, lattice structure not defined")
                if not self.sample.check_status(magdefs=True):
                    log.error("Cannot run simulation, no magnetic structure selected")
                if not self.sample.check_status(muon=True):
                    log.error("Cannot run simulation, muon position(s) not specified")
                    
                return False
                
        return True
    
    def getLargestLorentzSphere(self, sc):
        if self.sample:
            try:
                v = find_largest_sphere(self.sample, sc)
                self.newLorentzSphere.emit(v)
                
            except:
                log.error("Could not find largest sphere. " + \
                            "Did you define the lattice structure and " + \
                            "the muon positions?")
                #return 0.
        else:
            log.debug("mumble mumble...Sample not defined in getLargestLorentzSphere")
            #return 0.
        
        
    def load_structure_file(self):
        fname =select_file_dialog(self,action='open',check_read=True)
        
        if not fname:
            return
        
        
        _, fext = os.path.splitext(fname)
        if fext.lower() == '.xsf':
            load_xsf(self.sample,fname)
        elif fext.lower() == '.cif':
            load_cif(self.sample,fname)
        else:
            log.error("Invalid file extension for file {}. Only XSF and CIF file supported.".format(fname[0]))
            return
        
        
        self.crystal_interface.update_data(self.sample)
        self.muon_interface.update_data(self.sample)
        self.magnetism_interface.update_data(self.sample)
        
    def add_muons(self, muon_list):
        self.sample._reset(muon=True)
        allok = True
        counter = 0
        
        for mp in muon_list:
            try:
                self.sample.add_muon(mp)
                counter += 1
            except CellError:
                log.error("Lattice structure not defined! Could not add muons.")
                allok = False
                break
            except:
                log.debug("Could not add muons", exc_info=True)
                allok = False
                
        if allok:
            log.info("Muon positions updated. {} positions defined.".format(counter))
            
    
    def add_mag_order(self):
        try:
            self.sample.new_mm()
        except CellError:
            log.error("Lattice structure not defined! Could not create magnetic model.")
        except:
            log.error("Could not create magnetic model",
                          exc_info=True)
        else:
            log.info("New magnetic order added.")
            self.magnetism_interface.update_data(self.sample)
        
    def select_mag_order(self,index):
        try:
            self.sample.current_mm_idx = index
        except:
            log.error("Could not select magnetic model %d", index,
                          exc_info=True)
        else:
            log.info("Select magnetic model %d: %s ", index+1, self.sample.mm.desc)
            self.magnetism_interface.update_data(self.sample)
    
    def update_mm_desc(self,index, text):
        if index != self.sample.current_mm_idx:
            raise ValueError
        else:
            self.sample.mm.desc = text
            
    def update_mm_FC(self, index, data):
        if index != self.sample.current_mm_idx:
            raise ValueError
        else:
            try:
                self.sample.mm.fc_set(*data)
            except:
                log.error("Could not set Fourier Components",
                          exc_info=True)
            else:
                log.info("Updated Fourier Components")
            
        self.magnetism_interface.update_data(self.sample)
    
    def update_mm_K(self, index, K):
        if index != self.sample.current_mm_idx:
            raise ValueError
        else:
            try:
                self.sample.mm.k = K
            except:
                log.error("Could not set Propagation Vector",
                          exc_info=True)
            else:
                log.info("Updated Propagation Vector")
            
    
    def queue_calculation(self, data):
        if self.check('run_sum'):
            self.calc_queue.put((self.sample,data))
        else:
            return
            
        if self.calc_thread == None:
            self.calc_thread = CalcThread(self.calc_queue,self.res_queue)
            self.calc_thread.start()
            self.threadRunning.emit(True)

            
    
    def _new_empty_sample(self):
        self.sample = Sample()
        self.name = self.sample.name
        self.LESampleDesc.setText(self.sample.name)
        self._activate_buttons()
        

    
    def _load_sample(self):
        fname = select_file_dialog(self, check_read=True)
        if fname:
            try:
                self.sample = load_sample(fname)
            except:
                log.error("Could not load sample.", exc_info=True)
                return
        else:
            return
        
        self._update_interface()
        
    def _update_interface(self):
        #update interface
        self._activate_buttons()
        
        #update tabs
        self.crystal_interface.update_data(self.sample)
        self.magnetism_interface.update_data(self.sample)
        self.muon_interface.update_data(self.sample)
        
        #update name
        self.LESampleDesc.setText(self.sample.name)
        self.name = self.sample.name
        
    
    def _save_sample(self):
        fname = select_file_dialog(self, action='save',check_write=True)
        if not fname:
            return
        try:
            f = open(fname,'w')
            save_sample(self.sample, fileobj=f)
            f.close()
        except:
            log.error("Could not save sample.",exc_info=True)
            return
    

        
class QPlainTextEditLogger(log.Handler):
    def __init__(self, parent):
        super(QPlainTextEditLogger, self).__init__()
        self.ui = QtWidgets.QPlainTextEdit(parent)
        self.ui.setReadOnly(True)    
        self.ui.destroyed.connect(self._removeMe)

    def emit(self, record):
        msg = self.format(record)
        self.ui.appendPlainText(msg)
    
    def _removeMe(self):
        print("Removing log handler")
        log.getLogger().removeHandler(self)
        

#class MainWindow(QtWidgets.QMainWindow):
class MainWindow(QtWidgets.QMainWindow):
    """Load .ui file example, using setattr/getattr approach"""
    def __init__(self, parent=None, input_q=None, output_q=None):
        QtWidgets.QWidget.__init__(self, parent)
        
        self._running_threads = 0
        self._queues = (input_q, output_q)
        
        
        self.setObjectName("MueSR")
        self.resize(800, 500)
        self.centralwidget = QtWidgets.QWidget(self)
        
        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        
        self.centralwidget.setSizePolicy(sizePolicy)
        self.centralwidget.setObjectName("centralwidget")
        self.CntWdgVLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.CntWdgVLayout.setObjectName("CntWdgVLayout")
        
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName("tabWidget")
        self.CntWdgVLayout.addWidget(self.tabWidget)
        self.setCentralWidget(self.centralwidget)
        
        self.menubar = QtWidgets.QMenuBar(self)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 19))
        self.menubar.setObjectName("menubar")
        self.setMenuBar(self.menubar)
        
        self.statusbar = QtWidgets.QStatusBar(self)
        self.statusbar.setObjectName("statusbar")
        self.setStatusBar(self.statusbar)
        
        self.toolBar = QtWidgets.QToolBar(self)
        self.toolBar.setEnabled(True)
        self.toolBar.setObjectName("toolBar")
        self.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)


        
        self._timer = QtCore.QTimer(self)
        self._timer.timeout.connect(self._add_dots)
        self._timer.start(1000)
        self._dots = ''
        
        self.createMenusAndToolbars()
        
        
        
        self.tabWidget.setTabsClosable(True)
        self.tabWidget.tabCloseRequested.connect(self._del_sample)
        

        
        if running_inside_mantid and input_q != None:
            if input_q.empty():
                self._new_sample()
            else:
                while not input_q.empty():
                    self._new_sample(input_q.get())
        else:
            self._new_sample()
        
        logTextBox = QPlainTextEditLogger(self)
        # Format
        logTextBox.setFormatter(log.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        log.getLogger().addHandler(logTextBox)
        
        # set debug level
        log.getLogger().setLevel(log.DEBUG)
        log.captureWarnings(True)
        
        self.CntWdgVLayout.addWidget(logTextBox.ui)
        
        log.info("MuESR Gui started.")

    def createMenusAndToolbars(self):
        fileMenu = self.menuBar().addMenu('File')
        fileMenu.addAction('Save All', self.save_all)
        fileMenu.addAction('Quit', self.quit)
        fileMenu = self.menuBar().addMenu('Help')
        fileMenu.addAction('About', self.about)
        
        self.toolBar.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        
        new_ico = QtGui.QPixmap()
        new_ico.loadFromData(get_icon("new_sample"))
        
        
        self.toolBar.addAction(QtGui.QIcon(new_ico), "New sample", self._new_sample)

    
    def _del_sample(self, idx):
        
        self.tabWidget.removeTab(idx)
        
    
    def _new_sample(self, sample_obj = None):
        initial_tab = SampleTab(self, sample_obj)
        initial_tab.sampleNameChanged.connect(self.setTabTitle)
        initial_tab.threadRunning.connect(self.updateRunningThreads)
        initial_tab.results.mantidExportResult.connect(self._append_to_out_queue)
        
        if sample_obj:
            self.tabWidget.addTab(initial_tab, sample_obj.name)
        else:
            self.tabWidget.addTab(initial_tab, "New Sample")
        
    def _append_to_out_queue(self, obj):
        if self._queues[1] != None:
            self._queues[1].put(obj)
        

    def _create_toolbutton(self, parent=None, text = None, shortcut = None,
                          icon = None, tip = None, toggled = None,
                          triggered = None, autoraise = True,
                          text_beside_icon = False):
        ''' create an toolbutton '''
        button = QtWidgets.QToolButton(parent)
        if text is not None:
            button.setText(text)
        if icon is not None:
            icon = QtWidgets.getIcon(icon)
            button.setIcon(icon)
        if text is not None or tip is not None:
            button.setToolTip(text if tip is None else tip)
        if text_beside_icon:
            button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        button.setAutoRaise(autoraise)
        if triggered is not None:
            button.clicked.connect(triggered)
        if toggled is not None:
            button.toggled.connect(toggled)
            button.setCheckable(True)
        if shortcut is not None:
            button.setShortcut(shortcut)
        return button


    
    def setTabTitle(self, TabTitle):
        idx = self.tabWidget.currentIndex()
        self.tabWidget.setTabText(idx, TabTitle)
    
    def updateRunningThreads(self, running):
        if running:
            self._running_threads += 1
        else:
            self._running_threads -= 1
        
        self.statusbar.showMessage("New thread started. {} currently running.".format(self._running_threads) if running else "Thread stopped. {} currently running.".format(self._running_threads),1000)
    
    def _add_dots(self):
        if self._running_threads > 0:
            self._dots += '.'
            self.statusbar.showMessage("{} thread(s) currently running.".format(self._running_threads)+self._dots,1000)
        else:
            self._dots = ''
        
    def about(self):
        QtWidgets.QMessageBox.about(self,"MuESR","Created by Pietro Bonfa. Icons by ///")
    
    def save_all(self):
        log.info("Not implemented.")
        pass
      
    def quit(self):
        self.close()


#utility functions

def select_file_dialog(parent, name="Choose File", action='open', initdir="", check_read=False,check_write=False, ffilter='', add_ext_if_missing=True):
    filename = ""
    if action == 'open':
        filename = QtWidgets.QFileDialog.getOpenFileName(parent=parent, 
                    caption=name, directory=initdir,filter=ffilter)
    elif action == 'save':
        filename = QtWidgets.QFileDialog.getSaveFileName(parent=parent,
                    caption=name, directory=initdir,filter=ffilter)
    else:
        log.error("Internal error")
    
    if type(filename) in (list, tuple) :
        filename = filename[0]
    
    if not filename:
        return
        
    if action == 'save' and add_ext_if_missing and ffilter != '':
        fname, fext = os.path.splitext(filename)
        if fext != os.path.splitext(ffilter)[1]:
            filename = fname + os.path.splitext(ffilter)[1]
        
    
    if check_read:
        try:
            open(filename,'r').close()
        except:
            log.error("Cannot read file: {}. Check permissions".format(filename))
            return
    if check_write:
        try:
            open(filename,'w').close()
        except:
            log.error("Cannot write file: {}. Check permissions".format(filename))
            return 
    
    return filename

def main():

    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()

    window.show()
    app.exec_()

    app.exit()

if __name__ == '__main__':
    main()
