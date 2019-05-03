import vtk
from vtk import *
import string
import sys

f = open(sys.argv[1],'r')
outfile = sys.argv[2]

l = f.readline()
l.strip()
tokens = l.split(' ')

Nnodes = int(tokens[0])
Nelts = int(tokens[1])

Points = vtk.vtkPoints()

print("Reading %d nodes..." % Nnodes)
for i in range(Nnodes):
    line = f.readline()
    line.strip()
    d = line.split(' ')
    id = Points.InsertNextPoint(float(d[0]),float(d[1]),float(d[2]))

print("Done") 

Triangles = vtk.vtkCellArray()
Triangle = vtk.vtkTriangle()

print("Reading %d elements..." % Nelts)
for i in range(Nelts):
    line = f.readline()
    d = line.split()
    Triangle.GetPointIds().SetId(0,int(d[0])-1)
    Triangle.GetPointIds().SetId(1,int(d[1])-1)
    Triangle.GetPointIds().SetId(2,int(d[2])-1)
    Triangles.InsertNextCell(Triangle)

print("Done") 
f.close()


polydata = vtk.vtkPolyData()
polydata.SetPoints(Points)
polydata.SetPolys(Triangles)
polydata.Modified()
if vtk.VTK_MAJOR_VERSION <= 5:
    polydata.Update()
 
writer = vtk.vtkXMLPolyDataWriter();
writer.SetFileName(outfile);
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(polydata)
else:
    writer.SetInputData(polydata)
writer.SetDataModeToBinary()
writer.Write()
