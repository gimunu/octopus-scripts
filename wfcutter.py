#!/usr/bin/env python

# Copyright (C) 2015 Umberto De Giovannini 
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

from __future__ import division 

import argparse
import sys

try:
    import numpy as np
except ImportError:
    print "numpy is not installed"

try:
    from vtk import vtkStructuredPointsReader
except ImportError:
    print "vtk is not installed"
    
from vtk.util import numpy_support as VN


import vtk
import pylab as pl

# try:
#     from scipy import interpolate
# except ImportError:
#     print "scipy is not installed"


def display_vtk(filename):
    
    # Prepare to read the file
    reader = vtk.vtkStructuredPointsReader()
    reader.SetFileName( filename )
    # must call Update() before we fetch the dimensions
    reader.Update()

    # just for illustration:
    # get the extent of the data and print it
    W,H,D = reader.GetOutput().GetDimensions()
    # string formatting is similar to the sprintf style in C
    print "Reading '%s', width=%i, height=%i, depth=%i" %(filename, W, H, D)

    # create an outline of the dataset
    outline = vtk.vtkOutlineFilter()
    # outline.SetInput( reader.GetOutput() )
    outline.SetInputConnection( reader.GetOutputPort() )
    outlineMapper = vtk.vtkPolyDataMapper()
    # outlineMapper.SetInput( outline.GetOutput() )
    outlineMapper.SetInputConnection( outline.GetOutputPort() )
    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper( outlineMapper )

    # the actors property defines color, shading, line width,...
    outlineActor.GetProperty().SetColor(0.0,0.0,1.0)
    outlineActor.GetProperty().SetLineWidth(2.0)

    # Get the isosurface of our data
    isosurfExtractor = vtk.vtkContourFilter()
    isosurfExtractor.SetInputConnection(reader.GetOutputPort())
    isosurfExtractor.SetValue(0, 0.1)
    isosurfNormals = vtk.vtkPolyDataNormals()
    isosurfNormals.SetInputConnection(isosurfExtractor.GetOutputPort())
    isosurfNormals.SetFeatureAngle(60.0)
    isosurfMapper = vtk.vtkPolyDataMapper()
    isosurfMapper.SetInputConnection(isosurfNormals.GetOutputPort())
    isosurfMapper.ScalarVisibilityOff()
    isosurf = vtk.vtkActor()
    isosurf.SetMapper(isosurfMapper)

    # Add sphere

    # renderer and render window 
    ren = vtk.vtkRenderer()
    ren.SetBackground(.8, .8, .8)
    renWin = vtk.vtkRenderWindow()
    renWin.SetSize( 600, 600 )
    renWin.AddRenderer( ren )

    # render window interactor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow( renWin )

    # add the actors
    ren.AddActor( outlineActor )
    ren.AddActor( isosurf )
    renWin.Render()

    # create window to image filter to get the window to an image
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(renWin)

    # create png writer
    wr = vtk.vtkPNGWriter()
    # wr.SetInput(w2if.GetOutput())
    wr.SetInputConnection(w2if.GetOutputPort())

    # Python function for the keyboard interface
    # count is a screenshot counter
    count = 0
    def Keypress(obj, event):
        global count, iv
        key = obj.GetKeySym()
        if key == "s":
            renWin.Render()     
            w2if.Modified() # tell the w2if that it should update
            fnm = "screenshot%02d.png" %(count)
            wr.SetFileName(fnm)
            wr.Write()
            print "Saved '%s'" %(fnm)
            count = count+1
            # add your keyboard interface here
            # elif key == ...

    # add keyboard interface, initialize, and start the interactor
    iren.AddObserver("KeyPressEvent", Keypress)
    iren.Initialize()
    iren.Start()   
 
    

def read_vtk(filename, skipGrid = False):

    reader = vtkStructuredPointsReader()
    reader.SetFileName(filename)
    # reader.SetDataByteOrderToBigEndian()
    # reader.ReadAllVectorsOn()
    # reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()

    dim = data.GetDimensions()
    vec = list(dim)
    vec = [i-1 for i in dim]
    vec.append(3)
    print dim
    print vec
    print data.GetPointData()

    u = VN.vtk_to_numpy(data.GetPointData().GetArray('scalar_field'))
    # nodes_nummpy_array = VN.vtk_to_numpy(data.GetPoints().GetData())
    # print nodes_nummpy_array
    
    u = u.reshape(dim,order='F')

    x = np.zeros(data.GetNumberOfPoints())
    y = np.zeros(data.GetNumberOfPoints())
    z = np.zeros(data.GetNumberOfPoints())

    if not skipGrid:
        # Looping over all the point takes a while so we can do ot only once
        for i in range(data.GetNumberOfPoints()):
                x[i],y[i],z[i] = data.GetPoint(i)

    x = x.reshape(dim,order='F')
    y = y.reshape(dim,order='F')
    z = z.reshape(dim,order='F')
    
    
    return (u, [x,y,z])

######################## MAIN ################################

def main(args):

    VERSION = '0.0(alpha)'
    
    desc="""This utility allow to analyze wavefunctions performing different post-processing 
routines.        
"""

    epilog="""Examples:

To process input file(s):
\'%(prog)s file1 [file2 ...]\'
\n\n
    """

    parser = argparse.ArgumentParser(version='%s version %s' %(sys.argv[0],VERSION),
                                     description=desc,
                                     epilog=epilog,
                                     formatter_class=argparse.RawTextHelpFormatter)


#     parser.add_argument('-c', '--cut', action='store', metavar='p1x,p1y,p1z;p2x,p2y,p2z[;p3 ...]', default= None,
#     help=
# """Cut the bands structure along a segmented line connecting
# two, or more, points of the Brillouin zone in sequence.
# Each segment is defined by two points p1=(p1x,p1y,p1z) and
# p2=(p2x,p2y,p2z) and is parametrized as following:
# v = (p2-p1) * t + p1.
# The coordinates of each point are passed as a string in a
# column-separated triplet. For example:
#
# -c "p1x,p1y,p1z;p1x,p1y,p1z;[...]"
#
# Point coordinates can also be expressed as mathematical
# formula such as:
#
# -c "0, 0, 0 ; 1/2, cos(pi/4), sin(pi/4); 1, 1, 1"
#
# """)
#
#     parser.add_argument('-d', '--spacing', type=float, metavar='float', default= None,
#     help=
# """Define the output grid spacing. The default value is
# the minimum spacing of the input data.
# """
#     )
#
    parser.add_argument('-d','--display', action="store_true", default=False,
    help="Display the data isosurface and the integration volumes.")
    
    

    parser.add_argument('file', nargs='+')
    

    # parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    args = parser.parse_args()
    

    
    
    for file in args.file: 
            print file
            display_vtk(file)
            # (wf, grid) =read_vtk(file)
            # print wf
            # print grid
        


if __name__ == "__main__":
   main(sys.argv[1:])


