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
    from scipy import interpolate
except ImportError:
    print "scipy is not installed"


from math import *
import os

import warnings
warnings.simplefilter('error', UserWarning)


def have_kpoint_symmetries(file):
    path = os.path.dirname(os.path.abspath(file))
    if "static" in path:
        path = path +"/../"


    def search_string_in_file(string, f):
        if f is not None:
            for line in f:
                if string.lower() in line.lower():
                    return line
    
        return None

    # Look for clues on symmetry being used to generate the k-point mesh
    try: 
        filename = path + "/inp"
        f = open(filename)
    except:
        f = None
    if search_string_in_file("KpointsUseSymmetries",f):
        return True

    try:
        filename = path + "/static/info"
        f = open(filename)
    except:
        f = None

    if search_string_in_file("symmetry-reduced k-points",f):
        return True

    return False     
        
def import_bands_file(fname, reduced = False):

    bands = np.loadtxt(fname)

    # detect the dimensionality from the number of columns
    dim = int((bands.shape[1] - 1)/2)
    kcol0 = 0
    # kcol1 = dim
    if reduced:
        kcol0 = dim 
    kcol1 = kcol0 + dim    
        
    ecol = bands.shape[1]-1
        
    
    kk    = np.array(bands[:,kcol0:kcol1])
    Etmp  = np.array(bands[:,ecol])

    nbands=0
    for l in range(bands.shape[0]):
        if ((bands[0,kcol0:kcol1]==bands[l,kcol0:kcol1]).all()): 
            nbands += 1
    # get unique k-point elements along the axes
    kx = np.unique(bands[:, kcol0 + 0])
    if dim > 1:
        ky = np.unique(bands[:, kcol0 + 1])
        if dim > 2: 
            kz = np.unique(bands[:, kcol0 + 2])
        else:
            kz = np.array([0.0],dtype= float)            
    else:
        ky = np.array([0.0],dtype= float)
        kz = np.array([0.0],dtype= float)
        
    if (kz.shape[0] == 1): dim = 2 # we are effectively 2D 
    if (ky.shape[0] == 1 and kz.shape[0] == 1): dim = 1 # we are effectively 1D 

    print "Detected a %d-dimensional k-space with %d bands "%(dim, nbands)        


    E  =  np.zeros((kx.shape[0],ky.shape[0],kz.shape[0],nbands))
    nnE = np.zeros((kx.shape[0],ky.shape[0],kz.shape[0],nbands),dtype= int)
    nnE[:,:,:] = nbands-1

    for l in range(bands.shape[0]):
        k=bands[l,kcol0:kcol1]
        i = np.where(kx == k[0])
        if dim > 1:
            j = np.where(ky == k[1])
            if dim > 2:
                m = np.where(kz == k[2])
            else:
                m = 0                 
        else:
            j = 0
            m = 0  
        E[i,j,m, nnE[i,j,m]] = bands[l,ecol]
        nnE[i,j,m] -= 1
        


    for i in range(kx.shape[0]):
        for j in range(ky.shape[0]):
            for m in range(kz.shape[0]):
                E[i,j,m,:].sort()

    
    return (E, [kx, ky, kz], dim)

def write_gpl(fname, E, kk, nk, header = None, append = False):
        if append:
            out = open( fname+".gpl","a")
        else:
            out = open( fname+".gpl","w")
                        
        if header is not None:
            out.write("%s"%header)

        for i in range(nk[0]):
            if nk[1] > 1:
                for j in range(nk[1]):
                    out.write("%e\t%e"%(kk[i,0], kk[j,1]))
                    for ie in range(E.shape[2]):
                        out.write("\t%e"%(E[i,j,ie]))
                    out.write("\n")
            else:
                out.write("%e"%(kk[i,0]))
                for ie in range(E.shape[2]):
                    out.write("\t%e"%(E[i,0,ie]))
                out.write("\n")
            
            
            if nk[1] > 1 : out.write("\n") 
            
        out.write("\n")             
        out.close()
    

def line(v, p, u):
    return (v[0]*u + p[0], v[1]*u + p[1], v[2]*u + p[2])
    
def find_line_bounds(v, p, grid, dx):
    bounds = np.zeros([2,2])
    ubounds = np.zeros(2)
    
    bounds[0,:] =[min(grid[0][:]), max(grid[0][:])]
    bounds[1,:] =[min(grid[1][:]), max(grid[1][:])]
    

    iu = 0
    while True:
        u = iu * dx
        pnt = np.array(line(v,p,u))
        if ((pnt[:2] <= bounds[:2,0]).any() or (pnt[:2] >= bounds[:2,1]).any() ):
            break
        ubounds[1] = u
        # print pnt[:2],bounds[:2,0],bounds[:2,1]
        iu +=1

    iu = 0
    while True:
        u = iu * dx
        pnt = np.array(line(v,p,u))
        if ((pnt[:2] <= bounds[:2,0]).any() or (pnt[:2] >= bounds[:2,1]).any() ):
            break
        ubounds[0] = u
        # print pnt[:2],bounds[:2,0],bounds[:2,1]
        iu -=1
    return ubounds    

def grid_bounds(grid, dim):
    bounds = np.zeros([dim,2])

    for idim in range(dim):
        bounds[idim,:] =[min(grid[idim][:]), max(grid[idim][:])]
    
    return bounds
    
def points_in_bounds(pnts, bounds):
    
    for i in range(pnts.shape[0]):
        pnt = pnts[i,:]
        if ((pnt[:] < bounds[:,0]).any() or (pnt[:] > bounds[:,1]).any() ):
            return False
        if ((pnt[:] < bounds[:,0]).any() or (pnt[:] > bounds[:,1]).any() ):
            return False


    return True

def slice_on_line(E, kmesh, dim, nk, p1, p2, len0, spacing = None):
    # We assume in the following that the kspace is 2D
        
    nbands = E.shape[3]

    #Same spacing as the original grid unless otherwise specified
    dk = np.zeros(dim-1)
    for idim in range(dim-1):
        dk[idim] = np.abs(kmesh[idim][1]-kmesh[idim][0])
        
    
    
    v = (p2-p1)
    length = (np.sqrt(np.dot(v,v)))
    du = min(dk[:])
    if spacing:
        du = spacing
    u = np.linspace(0, 1, num = length/du)
    # print u
    
    Eout = np.zeros([u.shape[0],1, nbands])

    # Generate the grid-points on the line 
    kk = np.zeros([u.shape[0],dim])
    for iu in range(u.shape[0]):
        kk[iu,:] = np.array(line(v, p1 ,u[iu])[0:dim])

    # Check that the requested points are within the bounds of the input grid
    if (not points_in_bounds(kk,grid_bounds(kmesh,dim))):
        print('Error:Some of the requested output grid point is outside the bounds of the input grid.\nCurrent bounds are:\n%s'%grid_bounds(kmesh, dim))
        sys.exit(1)

    #Loop over all the bands
    for ib in range(nbands):
        if dim < 3:
            f = interpolate.RectBivariateSpline(kmesh[0][:], kmesh[1][:], E[:,:,0,ib])
            Eout[:,0,ib] = f.ev(kk[:,0],kk[:,1])
        else:
            # print (kmesh[0][:], kmesh[1][:], kmesh[2][:])
            # print kk
            
             
            f = interpolate.RegularGridInterpolator((kmesh[0][:], kmesh[1][:], kmesh[2][:]), E[:,:,:,ib])
            Eout[:,0,ib] = f(kk)
    
    len1 = len0 + length
    u = u * length + len0
    # print len0, len1, du
    # print u
    
    return (Eout, u, len1)      
        

######################## MAIN ################################

def main(args):

    VERSION = '0.0(alpha)'
    
    desc="""This utility converts the bands structure file generated by Octopus 
into a human-readable and easy-to-plot format.    
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


    parser.add_argument('-c', '--cut', action='store', metavar='p1x,p1y,p1z;p2x,p2y,p2z[;p3 ...]', default= None,
    help=
"""Cut the bands structure along a segmented line connecting 
two, or more, points of the Brillouin zone in sequence. 
Each segment is defined by two points p1=(p1x,p1y,p1z) and 
p2=(p2x,p2y,p2z) and is parametrized as following: 
v = (p2-p1) * t + p1.
The coordinates of each point are passed as a string in a 
column-separated triplet. By default reduced coordinate
units are use so each point coordinate should vary 
between []-1/2,1,2]. 
For example:

-c "p1x,p1y,p1z;p1x,p1y,p1z;[...]" 

Point coordinates can also be expressed as mathematical 
formula such as:

-c "0, 0, 0 ; 1/2, cos(pi/4)/2, sin(pi/4)/4; 0.5, 0.5, 0.5"

""")                      

    parser.add_argument('-d', '--spacing', type=float, metavar='float', default= None, 
    help=
"""Define the output grid spacing. The default value is 
the minimum spacing of the input data.
"""    
    )
    
    parser.add_argument('--absolute', action="store_true", default=False, 
    help="Use absolute coordinates in reciprocal space.")
    
    

    parser.add_argument('file', nargs='+')
    

    args = parser.parse_args()
    

    
    
    for file in args.file: 

        if have_kpoint_symmetries(file):
            warnings.warn(
""" 
 It seems that the k-point grid has been generated using symmetries.
 Since symmetries are not supported by this program the result are likely to 
 be garbage.
"""                
            )
        
        # (E,kx,ky) = import_file(file, write_out = (dir is None))
        (E,kmesh,dim) = import_bands_file(file,not args.absolute)
        nk = np.zeros(3)
        nk = [kmesh[0].shape[0],kmesh[1].shape[0],kmesh[2].shape[0]]

        # parse the cut nodal points
        if args.cut is not None:
            tokens = np.array(args.cut.split(';'))
            pts = np.zeros((tokens.shape[0],3))
            for i in range(tokens.shape[0]):
                coords  = tokens[i].split(',')
                for j in range(len(coords)):
                   pts[i,j] = eval(coords[j])  
            
            # print pts.shape[:]
            # dir = np.array(args.cut.split(','),dtype=float)
        
                
        header = "# This file contains %s bands.\n"%(E.shape[3])
        if args.cut is not None:
            header = "%s# Bands are evaluated along the following path in the Brillouin zone:\n# %s\n"%(header,args.cut)
            file = file+".cut"
            append = False
            lenght = 0.0
            for i in range(pts.shape[0]-1):
                 p1 = pts[i,:]     
                 p2 = pts[i+1,:]                      
                 print "Segment: %s --> %s"%(p1[0:dim],p2[0:dim])
                 header = "%s#\n# Slice on a line segment connecting %s and %s \n#\n"%(header,p1[0:dim],p2[0:dim])
                 (E_,kx_,lenght) = slice_on_line(E, kmesh, dim, nk, p1, p2, lenght, spacing = args.spacing)
                 nk_ = np.array((kx_.shape[0],0,0))
                 kk = np.zeros([max(nk_[:]),3])
                 kk[0:nk_[0],0] = kx_[0:nk_[0]]
                 kk[0:nk_[1],1] = kmesh[1][0:nk_[1]]
                 
        
                 write_gpl(file, E_, kk, nk_, header, append = append)
                 header = ""
                 append = True
        elif dim <= 2:
            E_ = E[:,:,0,:]        
            kk = np.zeros([max(nk[:]),3])
            kk[0:nk[0],0] = kmesh[0][0:nk[0]]
            kk[0:nk[1],1] = kmesh[1][0:nk[1]]

            write_gpl(file, E_, kk, nk, header)
        
        

    


if __name__ == "__main__":
   main(sys.argv[1:])