import sys
import vtk
import os 
from PyQt5 import QtCore, QtGui
from PyQt5 import Qt
import numpy 
# noinspection PyUnresolvedReferences
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import (
    vtkFloatArray,
    vtkIdList,
    vtkPoints
)
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkPolyData
)
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkCamera,
    vtkPolyDataMapper,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


class VtkQtFrame(Qt.QMainWindow):
    def __init__(self):
        Qt.QMainWindow.__init__(self)    
        self.layout = Qt.QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)

        self.setSizePolicy(Qt.QSizePolicy.Preferred, Qt.QSizePolicy.Preferred)

        self.VTKRenderer = vtkRenderer()   
        self.VTKRenderer.SetBackground(1, 1, 1)
        self.VTKRenderer.SetViewport( 0, 0, 1, 1)

        self.VTKRenderWindow = vtkRenderWindow()
        self.VTKRenderWindow.AddRenderer(self.VTKRenderer)    
        self.VTKRenderWindowInteractor = QVTKRenderWindowInteractor(self,rw=self.VTKRenderWindow)
        self.layout.addWidget(self.VTKRenderWindowInteractor)

        self.VTKCamera = vtkCamera()
        self.VTKCamera.SetClippingRange(0.1,1000)
        self.VTKRenderer.SetActiveCamera(self.VTKCamera)

        self.VTKInteractorStyleSwitch = vtk.vtkInteractorStyleSwitch()
        self.VTKInteractorStyleSwitch.SetCurrentStyleToTrackballCamera()
        self.VTKRenderWindowInteractor.SetInteractorStyle(self.VTKInteractorStyleSwitch)
        self.VTKRenderWindowInteractor.Initialize()
        self.VTKRenderWindowInteractor.Start()
        self.VTKRenderWindowInteractor.ReInitialize()
        


class SphereActor(vtkActor):
    def __init__(self,rad,res,r,c):
        self.pos = numpy.array(r)
        self.source = vtk.vtkSphereSource()
        self.source.SetRadius(rad)
        self.source.SetPhiResolution(res)
        self.source.SetThetaResolution(res)
        self.source.SetCenter(r[0],r[1],r[2])
        self.Mapper = vtkPolyDataMapper()
        #self.Mapper.SetInput(self.source.GetOutput())
        

        self.Mapper.SetInputConnection(self.source.GetOutputPort())

        ## Create source
        #source = vtk.vtkSphereSource()
        #source.SetCenter(0, 0, 0)
        #source.SetRadius(5.0)

        ## Create a mapper
        #mapper = vtk.vtkPolyDataMapper()
        #mapper.SetInputConnection(source.GetOutputPort())

        ## Create an actor
        #actor = vtk.vtkActor()
        #actor.SetMapper(mapper) 
                
        self.SetMapper(self.Mapper)
        self.GetProperty().SetColor((c[0],c[1],c[2]))
    def move_to(self,r):
        self.pos = numpy.array(r)
        self.source.SetCenter(r[0],r[1],r[2])    
    def set_color(self,color):
        self.GetProperty().SetColor(color)
    def set_rad(self,rad):
        self.source.SetRadius(rad)    
    def get_pos(self):
        return self.pos
    
    
def main():
    app = Qt.QApplication(sys.argv)

    #create our new Qt MainWindow
    window = Qt.QMainWindow()

    #create our new custom VTK Qt widget
    render_widget = VtkQtFrame()

    for i in range(0,10):
        # random 3D position between 0,10
        r = numpy.random.rand(3)*10.0
        # random RGB color between 0,1
        c = numpy.random.rand(3)
        # create new sphere actor
        my_sphere = SphereActor(1.0,20,r,c)
        # add to renderer
        render_widget.VTKRenderer.AddActor(my_sphere)

    # reset the camera and set anti-aliasing to 2x
    render_widget.VTKRenderer.ResetCamera()
    #render_widget.VTKRenderWindow.SetAAFrames(2)

    # add and show
    window.setCentralWidget(render_widget)
    window.show()

    # start the event loop
    try:
        sys.exit(app.exec_())
    except SystemExit as e:
        if e.code != 0:
            raise()
        os._exit(0)
    return 

if __name__ == '__main__':
    main()