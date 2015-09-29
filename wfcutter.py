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


from math import *

import vtk
import pylab as pl


# list of supported geometries
geometries = ["sphere", 
              "cylinder", 
              "parallelepiped"]
              
parameters = [{'center':[0,3],'radius': 3 }, 
              {'center':[0,3],'axis'  :[3,6],'radius': 6 , 'height': 7},
              {'center':[0,3],'axis'  :[3,6],'sides': [6,9]}] 


def display_vtk(filename, geometry = None, isolevel = 0.1):
    
    # Prepare to read the file
    reader = vtk.vtkStructuredPointsReader()
    reader.SetFileName( filename )
    # must call Update() before we fetch the dimensions
    reader.Update()

    # just for illustration:
    # get the extent of the data and print it
    W,H,D = reader.GetOutput().GetDimensions()
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
    isosurfExtractor.SetValue(0, isolevel)
    isosurfNormals = vtk.vtkPolyDataNormals()
    isosurfNormals.SetInputConnection(isosurfExtractor.GetOutputPort())
    isosurfNormals.SetFeatureAngle(60.0)
    isosurfMapper = vtk.vtkPolyDataMapper()
    isosurfMapper.SetInputConnection(isosurfNormals.GetOutputPort())
    isosurfMapper.ScalarVisibilityOff()
    isosurfActor = vtk.vtkActor()
    isosurfActor.SetMapper(isosurfMapper)

    # Add geometries
    geoActors= []
    for geo in geometry:
        par = parameters[geometries.index(geo['shape'])]     

        if geo['shape'] == 'sphere':            
            source = vtk.vtkSphereSource()
            source.SetCenter(geo['parameters'][par['center'][0]:par['center'][1]])
            source.SetRadius(geo['parameters'][par['radius']])
            source.SetThetaResolution(50)
            source.SetPhiResolution(50)
            
            # mapper
            mapper = vtk.vtkPolyDataMapper()
            if vtk.VTK_MAJOR_VERSION <= 5:
                mapper.SetInput(source.GetOutput())
            else:
                mapper.SetInputConnection(source.GetOutputPort())
            # actor
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetOpacity(0.5)
            actor.GetProperty().SetColor(1.0,0.5,0.5)
            geoActors.append(actor)

        if geo['shape'] == 'cylinder':
            source = vtk.vtkCylinderSource()
            source.SetCenter(geo['parameters'][par['center'][0]:par['center'][1]])
            source.SetRadius(geo['parameters'][par['radius']])
            source.SetHeight(geo['parameters'][par['height']])
            source.SetResolution(100)
            
            # mapper
            mapper = vtk.vtkPolyDataMapper()
            if vtk.VTK_MAJOR_VERSION <= 5:
                mapper.SetInput(source.GetOutput())
            else:
                mapper.SetInputConnection(source.GetOutputPort())
            # actor
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetOpacity(0.5)
            actor.GetProperty().SetColor(1.0,0.5,0.5)
            geoActors.append(actor)    

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
    ren.AddActor( isosurfActor )
    for act in geoActors:
        ren.AddActor( act )
        
    renWin.Render()

    # create window to image filter to get the window to an image
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(renWin)

    # create png writer
    wr = vtk.vtkPNGWriter()
    wr.SetInputConnection(w2if.GetOutputPort())

    # Python function for the keyboard interface
    # count is a screenshot counter
    count = 0
    def Keypress(obj, event):
        global count
        key = obj.GetKeySym()
        if key == "s":
            renWin.Render()     
            w2if.Modified() # tell the w2if that it should update
            fnm = "screenshot%02d.png" %(count)
            wr.SetFileName(fnm)
            wr.Write()
            print "Saved '%s'" %(fnm)
            count = count+1
        if key == "q":
            #quit
            print "Exit"
            exit()                

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
    # print dim
    # print vec
    # print data.GetPointData()

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

def get_geometry(string):
    """ Parse geometry properties from sting"""
    
    # generate the short names
    geometries_abbr= [x[0:3] for x in geometries]

    #lowercase and split on ";"
    shapetokens = np.array(string.lower().split(';'))
    geometry = []
    for string in shapetokens:
        # get the name of the geometrical shape
        geotype = string[0:string.find("(")].strip()
        # is this a recognized geometry?
        if not geotype in geometries and not geotype in geometries_abbr:
            raise Exception("\"%s\" is not a recognized geometry."%(geotype))
        
        # get the string between rounded brackets
        values = string[string.find("(")+1:string.find(")")].split(',')
        values = [eval(val) for val in values] # evaluate expressions contained in each component

        # get full-name if short
        if geotype in geometries_abbr:
            geotype = geometries[geometries_abbr.index(geotype)]
            
        geo = {"shape": geotype, "parameters": values}
        geometry.append(geo)
        
    return geometry    
    
def integrateOverVolume(func, grid, geometry):
    print grid[0].shape[:]
    print grid[1].shape[:]
    print grid[2].shape[:]

    print grid[0][1][0][0],  grid[0][0][0][0]
    print grid[1][0][1][0],  grid[1][0][0][0]
    print grid[2][0][0][1],  grid[2][0][0][0]

    
    spacing=np.zeros(3)
    for idim in range(3):
        spacing[idim] = abs(grid[idim][1][0][0]-grid[idim][0][0][0])
        
    print spacing
    assert all(x == spacing[0] for x in spacing), "Cannot integrate! The grid appears to be not equally spaced."
    
    res = np.zeros(len(geometry))
    print res
    
    print    grid[0][grid[0].shape[0]/2]

    
    tot = sum(res[:])
    
    return (tot, res)
    
    
##############################################################    
######################## MAIN ################################
##############################################################    

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


    parser.add_argument('-g', '--geometry', action='store', metavar='sphere(cx,cy,cz,R) [;sphere ...]', default= None,
    help=
""" Define the geometry where the charge is integrated. 
The available geometries are sphere, cylinder and parallelepiped.

""")

    parser.add_argument('--isolevel', type=float, metavar='float', default= 0.1,
    help=
"""Define the value of at wich calculate the iso-surface (default 0.1).
"""
    )

    parser.add_argument('-d','--display', action="store_true", default=False,
    help="Display the data isosurface and the integration volumes.")
    
    

    parser.add_argument('file', nargs='+')
    

    # parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    args = parser.parse_args()
    
    geometry = None
    if args.geometry:
        geometry = get_geometry(args.geometry)
        
    for geo in geometry:    
        # print geo
        par = parameters[geometries.index(geo['shape'])]        
        # print par['center'][0]
        # print geo['parameters'][par['center'][0]:par['center'][1]]
        # print geo['parameters'][par['radius']]

    
    multiplefiles = False 
    for file in args.file: 
            print file
            if args.display:
                display_vtk(file, geometry, args.isolevel)
            if multiplefiles:    
                wf =read_vtk(file, skipGrid = True)
            else:    
                (wf, grid) =read_vtk(file)
                
            integral = integrateOverVolume(wf, grid, geometry)
                
            multiplefiles = True


if __name__ == "__main__":
   main(sys.argv[1:])


