
import os
from muesr.gui.mantidtools import running_inside_mantid

if running_inside_mantid:
    os.environ['QT_PREFERRED_BINDING'] = os.pathsep.join(['PyQt4'])
    
from .Qt import QtCore, QtWidgets
    
    
have_matplotlib = True
try:
    if QtCore.PYQT_VERSION < 0x50000: # Qt5
        from matplotlib.backends.backend_qt4agg \
            import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qt4agg \
            import NavigationToolbar2QT as NavigationToolbar
    else:
        from matplotlib.backends.backend_qt5agg \
            import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qt5agg \
            import NavigationToolbar2QT as NavigationToolbar

    from matplotlib.figure import Figure
except ImportError:
    have_matplotlib = False
    print("matplotlib not found")

import numpy as np


COLORS=["#1f271b","#da344d","#7fb069","#ffd166","#19647e"]

if have_matplotlib:
    class MplCanvas(FigureCanvas):
    
        def __init__(self, parent=None):
          
            self._last_handles = []
            
            self.fig = Figure()
            self.ax = self.fig.add_subplot(111)
            super(MplCanvas, self).__init__(self.fig)
            if parent:
                self.setParent(parent)
                
            
            #FigureCanvas.setSizePolicy(self,
            #        QSizePolicy.Expanding,
            #        QSizePolicy.Expanding)
            #FigureCanvas.updateGeometry(self)
            
            
        def clear(self):
            self.ax.clear()
        
        
        def text(self, txt):
            self.ax.text(0.5, 0.5,txt, horizontalalignment='center',
                          verticalalignment='center',transform=self.ax.transAxes)
    
        
        def plot_bars(self,x,y,c):
            self._last_handles = self.ax.bar(
                        x+np.linspace(-0.2,0.2,len(y),endpoint=False), # x vals, centered
                        y, color=c, width=0.4/len(y))
        
        def plot_lines(self,x,ys,c):
            self._last_handles = []
            for i, y in enumerate(ys):
                self._last_handles += self.ax.plot(x, y, 'o-',color=c[i])
        
        def set_xylabels_muon_sites(self, nmu):
            self.ax.set_ylabel("Field (T)")
            self.ax.set_xticks(np.arange(nmu)+1)
            self.ax.set_xticklabels(["Site {}".format(x) for x in range(1,nmu+1)])
        
        def legend(self,lgnd_list, minmax=False):
            if self._last_handles:
                if minmax:
                    self.ax.legend(self._last_handles[::2], lgnd_list, 
                              loc='upper center', bbox_to_anchor=(0.5, 1.07),
                              ncol=2, fancybox=True, shadow=True)
    
                else:
                    self.ax.legend(self._last_handles, lgnd_list, 
                              loc='upper center', bbox_to_anchor=(0.5, 1.07),
                              ncol=4, fancybox=True, shadow=True)
        
        
class MplWidget(QtWidgets.QWidget):
    def __init__(self, parent = None):

        QtWidgets.QWidget.__init__(self)
        self.canvas = MplCanvas()
        
        self.vbl = QtWidgets.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        
        # radio buttons with type of plot
        
        self.BtnBars = QtWidgets.QRadioButton('Bar Plot of |B|')
        self.BtnLines = QtWidgets.QRadioButton('Lines Plot of |B|')
        
        hbl = QtWidgets.QHBoxLayout()
        
        hbl.addWidget(QtWidgets.QLabel("Select plot",self))
        hbl.addWidget(self.BtnBars)
        hbl.addWidget(self.BtnLines)

        self.group = QtWidgets.QButtonGroup()
        self.group.addButton(self.BtnBars,0)
        self.group.addButton(self.BtnLines,1)  
        
        
        
        self.group.buttonClicked.connect(self._plot)
        
        self.vbl.addLayout(hbl)
        
        self.setLayout(self.vbl)
        
        self._LFObjs=None
        
    def set_LFs(self, LFObjs):
        self._LFObjs=LFObjs
        #self._plot()
    
    def _calc_norm(self, idx, minmax=False):
        if not self._LFObjs:
            return None
            
        if minmax:
            L = np.linalg.norm(self._LFObjs[idx].L, axis = 1)
            D = np.linalg.norm(self._LFObjs[idx].D, axis = 1)
            C = np.linalg.norm(self._LFObjs[idx].C, axis = 1)
            T = np.linalg.norm(self._LFObjs[idx].T, axis = 1)
                                            
            return [L.min(),L.max(), D.min(),D.max(), \
                    C.min(),C.max(),T.min(),T.max()]
        else:
            return [np.linalg.norm(self._LFObjs[idx].L), \
                            np.linalg.norm(self._LFObjs[idx].D), \
                            np.linalg.norm(self._LFObjs[idx].C), \
                            np.linalg.norm(self._LFObjs[idx].T)]
    
    def _plot(self, btnObj):
        
        plotnum = self.group.checkedId()

        if self._LFObjs == None:
            return

        nmu = len(self._LFObjs)
        if nmu == 0:
            return

        #get dimensionality
        ndim = self._LFObjs[0].T.ndim
        print('ndim',ndim)
        if plotnum == 0:
            if nmu > 30:
                self.canvas.clear()
                self.canvas.text("Too many mu sites")
                self.canvas.draw()
                return
                
            self.canvas.clear()
            x_pos = np.arange(nmu) + 1.
            
            # one dimensional plot
            minmax = False
            if ndim == 1:
                minmax = False 
            elif ndim == 2:
                minmax = True # plot min and max if many angles are present
            else:
                raise ValueError # what the hell is going on!!
                
            #plot all components for each site
            for i in range(nmu):
                self.canvas.plot_bars(x_pos[i],
                                      self._calc_norm(i, minmax),
                                      [COLORS[int(i/2)] for i in range(8)] if minmax else COLORS[:4])
            if minmax:
                self.canvas.legend(["Lor. (min/max)","Dip. (min/max)","Cont. (min/max)","Tot. (min/max)"],True)
            else:
                self.canvas.legend(["Lor.","Dip.","Cont.","Tot."],False)
            self.canvas.set_xylabels_muon_sites(nmu)
            self.canvas.draw()

        # line plot            
        elif plotnum == 1:
            
            if nmu > 100:
                self.canvas.clear()
                self.canvas.text("Too many mu sites")
                self.canvas.draw()
                return
            
            # should plot minmax?
            minmax = False
            if ndim == 1:
                minmax = False 
            elif ndim == 2:
                minmax = True # plot min and max if many angles are present
            else:
                raise ValueError # what the hell is going on!!
                
            x_pos = np.arange(nmu) + 1.
            if minmax:
                L = np.zeros([nmu,2])
                D = np.zeros([nmu,2])
                C = np.zeros([nmu,2])
                T = np.zeros([nmu,2])
                
                for i in range(nmu):
                    L[i,0],L[i,1],D[i,0],D[i,1], \
                    C[i,0],C[i,1],T[i,0],T[i,1] = self._calc_norm(i, True)
                
            else:
                L = np.zeros_like(x_pos)
                D = np.zeros_like(x_pos)
                C = np.zeros_like(x_pos)
                T = np.zeros_like(x_pos)
                
                for i in range(nmu):
                    L[i],D[i],C[i],T[i] = self._calc_norm(i)
                    
                  
                
            self.canvas.clear()
            
            if minmax:
                self.canvas.plot_lines(x_pos,(L[:,0],L[:,1],D[:,0],D[:,1],C[:,0],C[:,1],T[:,0],T[:,1]),[COLORS[int(i/2)] for i in range(8)])
                self.canvas.legend(["Lor. (min/max)","Dip. (min/max)","Cont. (min/max)","Tot. (min/max)"], True)
            else:
                self.canvas.plot_lines(x_pos,(L,D,C,T),COLORS[:4])
                self.canvas.legend(["Lor.","Dip.","Cont.","Tot."], False)
            self.canvas.set_xylabels_muon_sites(nmu)
            self.canvas.draw()

                    
                    
            

class DataPlotDialog(QtWidgets.QDialog):
    def __init__(self, parent = None):
        super(DataPlotDialog, self).__init__(parent)

        
        self.btnExport = QtWidgets.QPushButton('Export as Figure')
        self.btnClose = QtWidgets.QPushButton('Close')
        self.btnClose.clicked.connect(self.close)
        
        if have_matplotlib:
            self.canvas = MplWidget(self)
        else:
            self.canvas = QtWidgets.QLabel("Matplotlib not found. No visualization possible.")
            self.btnExport.setEnabled(False)
        
        HL = QtWidgets.QHBoxLayout()
        HL.addWidget(self.btnClose)
        HL.addWidget(self.btnExport)
        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Preferred)
        self.btnClose.setSizePolicy(sizePolicy)
        self.btnExport.setSizePolicy(sizePolicy)
        
        VL = QtWidgets.QVBoxLayout()
        VL.addWidget(self.canvas)
        VL.addLayout(HL)

        self.setLayout(VL)
        
    def set_LFs(self, LFObj):
        self.canvas.set_LFs(LFObj)


# following stuff  is just for develop and debug.
class MW(QtWidgets.QMainWindow):
    def __init__(self, parent = None):
        super(MW, self).__init__(parent)
        self.button = QtWidgets.QPushButton('Plot')
        self.button.clicked.connect(self._openWin)
        
        self.main_widget = QtWidgets.QWidget(self)
        self.setCentralWidget(self.main_widget)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.button)
        self.main_widget.setLayout(layout)
        
        
    def _openWin(self):
        w = DataPlotDialog(self)
        w.exec_()
        
        

if __name__ == "__main__":
    import sys
    qApp = QtWidgets.QApplication(sys.argv)
    main_window = MW()
    main_window.show()
    sys.exit(qApp.exec_())
