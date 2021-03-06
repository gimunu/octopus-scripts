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
import time

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
# warnings.simplefilter('error', UserWarning)


# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def update_progress(progress):
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


def get_prj_root_from_file(file):
    path = os.path.dirname(os.path.abspath(file))
    if "static" in path:
        path = path +"/../"

    return path

def search_string_in_file(string, f, rewind = True, commentchar = None):
    if f is not None:
        for il, line in enumerate(f):
            if line[0] == commentchar:
                continue
            if string.lower() in line.lower():
                if rewind:
                    f.seek(0, 0)
                return il, line

        if rewind:
            f.seek(0, 0)
        return None

    return None

def have_kpoint_symmetries(file):

    path = get_prj_root_from_file(file)

    # Look for clues on symmetry being used to generate the k-point mesh
    try:
        filename = path + "/static/info"
        f = open(filename)
    except:
        f = None
    if search_string_in_file("symmetry-reduced k-points",f):
        il, line = search_string_in_file("symmetry-reduced k-points",f)
        nsymm =  int(line.split()[-1])
        il, line = search_string_in_file("Total number of k-point",f)
        nk =  int(line.split()[-1])
        if nk != nsymm :
            f.close()
            return True

    try: 
        filename = path + "/inp"
        f = open(filename)
    except:
        f = None
    if search_string_in_file("KpointsUseSymmetries",f):
        f.close()
        return True

    return False     


def get_kweight_from_info(path):
    fname = path +"./static/info"
    try:
        f = open(fname)
    except:
        return None 
        

    # Don't proceed any further if symmetries have been used
    if search_string_in_file("symmetry-reduced k-points",f):
        f.close()
        return None
        
    il, line= search_string_in_file("Total number of k-points", f)
    nk = int(line.split('=')[1].strip(' \n'))
    
    kweight = np.array([], dtype=np.float)
    ilstart = il + 5 
    ilend   = ilstart + nk + 1
    for i, line in enumerate(f):
        if i > ilstart and i < ilend :
            kweight = np.append(kweight, np.float(line.strip(' \n').split()[4]))
    f.close()
    return kweight

def get_kweight_from_inp(path):
    fname = path +"./inp"
    
    try:
        f = open(fname)
    except:
        return None 
    try:
        istart, line= search_string_in_file("%KPointsReduced", f, commentchar = '#',rewind = False)
    except:
        f.seek(0, 0)    
        try: 
            istart, line= search_string_in_file("%KPoints", f, commentchar = '#',rewind = False)
        except: 
            return None
        
    
    iend, line= search_string_in_file("%", f)

    iend = istart + iend + 1

    kweight = np.array([], dtype=np.float)
    for i, line in enumerate(f):  
        #skip comment  
        if line.strip().startswith("#"):
            continue
        if i > istart and i < iend :
            kweight = np.append(kweight, np.float(line.strip(' \n').split('|')[0]))
    f.close()
    
    kweight = kweight if kweight.shape[0] > 0 else None
    return kweight
     


def find_custom_kpoint_indices(fname, nk):
    kweight = None
    path = get_prj_root_from_file(fname)
    
    if 'info' in fname:
        kweight = get_kweight_from_info(path)
    else:     
        kweight = get_kweight_from_inp(path)
        
    # if kweight is not None:
    #     # if np.count_nonzero(kweight) == 0 :
    #     #     # it means that the userdefined kpoints only contain the zero-weight path
    #     #     n0 = kweight.shape[0]
    #     #     kweight = np.zeros(nk)
    #     #     kweight[0:nk-n0] = 1
    #     pass
    # else:
        # kweight = get_kweight_from_info(path)

    idx0 = np.nonzero(kweight == 0)[0]

    return idx0.astype(np.int)+1
    
def import_eigenvalues_file(fname, abs_coords, use_mesh_grid = False):
    
    use_mesh_grid = use_mesh_grid or abs_coords
    
    f = open(fname)
    lstart, line = search_string_in_file("#k =", f, rewind = False)

    nbands, line = search_string_in_file("#k =", f)

    # Figure out spin polarization
    for il, line in enumerate(f):
        if il == lstart + 1:

            if len(line.split()) >= 7:
                #spinors
                 sidx = [4,7]
            else:
                if line.split()[1] == "--":
                    #spin unpolarized
                    sidx = []
                else:
                    #spin polarized
                    sidx = [1,2]
    f.seek(0,0)

    
    lstart, line = search_string_in_file("#k =", f, rewind = False)
    k = line.strip(' \n)').split("=")[2].split()[1:4]
    dim = len(k)
    for i,kk in enumerate(k):
        k[i] = np.float(kk.strip(','))
    kmesh  = np.array(k, dtype=np.float)
    ik = np.int(line.strip(' \n)').split('=')[1].split()[0].strip(','))
    kidx = np.array([ik], dtype = np.int)
    
    E  = np.array([], dtype=np.float)
    nkpt = 1

    ist = nbands
    skip = True
    for il, line in enumerate(f):
        if line == "\n" or "Fermi" in line:
            break
        if il == 0:
            if len(sidx) > 0 :
                spin = np.array(line.split()[sidx[0]:sidx[1]], dtype=np.float)
            else: 
                spin = None 
                
        if ist > 0:
            ist -= 1

            E = np.append(E, np.float(line.split()[2]))
            if len(sidx) > 0 and il > 0 :
                    spin = np.vstack((spin, line.split()[sidx[0]:sidx[1]] ))
            if skip:
                skip = False
            else:
                kmesh = np.vstack((kmesh, k))
                kidx = np.append(kidx, ik)
        else:     
            nkpt += 1
            k = line.strip(' \n)').split("=")[2].split()[1:4]
            for i,kk in enumerate(k):
                k[i] = np.float(kk.strip(','))
            ik = np.int(line.strip(' \n)').split('=')[1].split()[0].strip(','))
            ist = nbands 

    if abs_coords:        
        f.seek(0, 0)
        lstart, line = search_string_in_file("Reciprocal-Lattice Vectors", f, rewind = False)
        rlattice = np.zeros([3,3])
    
        for il, line in enumerate(f):
            if il > 2:
                break
            rvec = line.strip('\n').split()
            for j,val in enumerate(rvec):
                rlattice[il,j]=np.float(val)
        # multiply each kpoint for the reciprocal lattice vectors
        kmesh = np.dot(kmesh, rlattice)


    f.close()
    
    
    
    # find_custom_kpoint_indices(get_prj_root_from_file(fname))
    idx0  = find_custom_kpoint_indices(fname, kmesh.shape[0])
    if len(idx0)>0:
        # generate a mask indicating the indices of kidx that are contained in idx0 
        mask = np.in1d(kidx, idx0)
        E0 = E[mask]
        kmesh0 = kmesh[mask,:]
        if len(sidx):
            spin0 = spin[mask,:]
        else:
            spin0 = None    
        mask = np.logical_not(mask) 
        E  = E[mask]
        kmesh = kmesh[mask,:]
        if len(sidx):
            spin = spin[mask,:]

            
    if dim == 1:
        (kx,ky,kz,dim) = refine_dims_and_ks(dim, kmesh[:,0])
    elif dim == 2:
        (kx,ky,kz,dim) = refine_dims_and_ks(dim, kmesh[:,0], ky=kmesh[:,1])
    elif dim == 3:
        (kx,ky,kz,dim) = refine_dims_and_ks(dim, kmesh[:,0], ky=kmesh[:,1],kz=kmesh[:,2])
    
    
    if use_mesh_grid:
        # select from kmesh with nbands stride
        kx, ky, kz = kmesh[1::nbands,0], kmesh[1::nbands,1], kmesh[1::nbands,2]
        E = get_Ei_from_bands (kx,ky,kz, kmesh, nbands, E, dim, spin)  
    else:
        E = get_Eijm_from_bands (kx,ky,kz, kmesh, nbands, E, dim, spin)
    
    
    if len(idx0)>0:
        if spin0 is not None :
            for idim in range(3):
                E0 = np.vstack((E0,spin0[:,idim]))
            E0 = np.transpose(E0)
        else:
            E0tmp = np.zeros([E0.shape[0],1])
            E0tmp[:,0] =E0[:]
            E0 = E0tmp 
            
        return (E, [kx, ky, kz], dim, E0, kmesh0, nbands)

    return (E, [kx, ky, kz], dim)
  


def unique_or_vec(vec, use_mesh_grid, nkpt):
    if (use_mesh_grid):
        return vec[0:nkpt]    
    else:
        return np.unique(vec)
  
def refine_dims_and_ks(dim, kx, ky=None, kz=None, use_mesh_grid = False, nkpt=None):
    
    # get unique k-point elements along the axes
    kx = unique_or_vec(kx, use_mesh_grid, nkpt) #np.unique(kx)
    if dim > 1:
        ky = unique_or_vec(ky, use_mesh_grid, nkpt) #np.unique(ky)
        if dim > 2: 
            kz = unique_or_vec(kz, use_mesh_grid, nkpt) #np.unique(kz)
        else:
            kz = np.array([0.0],dtype= np.float)            
    else:
        ky = np.array([0.0],dtype= np.float)
        kz = np.array([0.0],dtype= np.float)
                    
    if (kz.shape[0] == 1 or all(kzi == 0.0 for kzi in kz)): dim = 2 # we are effectively 2D
    if (ky.shape[0] == 1 and kz.shape[0] == 1): dim = 1 # we are effectively 1D 
                
    return (kx,ky,kz, dim)

def get_Eijm_from_bands (kx,ky,kz,kmesh, nbands, bands, dim, other = None ):
    
    if other is None:
        nother = 0
    else:
        assert(other.shape[0] == bands.shape[0])
        nother = other.shape[1]
    
    E  =  np.zeros((kx.shape[0],ky.shape[0],kz.shape[0],nbands, nother+1))
    nnE = np.zeros((kx.shape[0],ky.shape[0],kz.shape[0],nbands),dtype= int)
    nnE[:,:,:] = nbands-1

    for l in range(bands.shape[0]):
        k=kmesh[l,:]
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
        E[i,j,m, nnE[i,j,m], 0] = bands[l]
        for n in range(nother):
            # add other (i.e. spin) columns 
            E[i,j,m, nnE[i,j,m], n+1] = other[l,n]
        nnE[i,j,m] -= 1
    


    for i in range(kx.shape[0]):
        for j in range(ky.shape[0]):
            for m in range(kz.shape[0]):
                # we have to sort all the other columns contaitning the spin 
                # according only to energy ordering
                E[i,j,m,:,:] = E[i,j,m, np.argsort(E[i,j,m,:,0])]

    return E
   
def get_Ei_from_bands (kx,ky,kz,kmesh, nbands, bands, dim, other = None ):

    if other is None:
        nother = 0
    else:
        assert(other.shape[0] == bands.shape[0])
        nother = other.shape[1]
    
    E  =  np.zeros((kx.shape[0],nbands, nother+1))
    nnE = np.zeros((kx.shape[0],nbands),dtype= int)
    nnE[:] = nbands-1

    if dim == 1: kkmat = np.array([kx])
    if dim == 2: kkmat = np.array([kx,ky])
    if dim == 3: kkmat = np.array([kx,ky,kz])

    kkmat=np.transpose(kkmat)
    for l in range(bands.shape[0]):
        kvec=kmesh[l,0:dim]
        i = np.argwhere( np.all((kkmat - kvec)==0, axis=1))[0][0]

        E[i, nnE[i], 0] = bands[l]
        for n in range(nother):
            # add other (i.e. spin) columns 
            E[i, nnE[i], n+1] = other[l,n]
        nnE[i] -= 1
    


    for i in range(kx.shape[0]):
        # we have to sort all the other columns containing the spin 
        # according only to energy ordering
        E[i,:,:] = E[i, np.argsort(E[i,:,0])]

        update_progress(i/(kx.shape[0]-1))
    
    print ""
    return E
    

   
def get_nkp_from_bands_file(fname):
    f = open(fname)
    if f is not None:
        for il, line in enumerate(f):
            if line[0] == "#":
                continue
            if line== '\n':
                f.close()
                return il-1
        f.close()        
        return None

    return None
          
def import_bands_file(fname, abs_coords, use_mesh_grid = False):

    use_mesh_grid = use_mesh_grid or abs_coords

    bands = np.loadtxt(fname)

    # detect the dimensionality from the number of columns
    dim = int((bands.shape[1] - 1)/2)
    kcol0 = 0
    # kcol1 = dim
    if not abs_coords:
        kcol0 = dim 
        
    kcol1 = kcol0 + dim    
        
    ecol = bands.shape[1]-1

    nbands=0
    for l in range(bands.shape[0]):
        if ((bands[0,kcol0:kcol1]==bands[l,kcol0:kcol1]).all()):
            nbands += 1
    
    nkpt = get_nkp_from_bands_file(fname)
        
                      
    if dim == 1:
        (kx,ky,kz,dim) = refine_dims_and_ks(dim, bands[:, kcol0 + 0], use_mesh_grid=use_mesh_grid, nkpt=nkpt)
    elif dim == 2:    
        (kx,ky,kz,dim) = refine_dims_and_ks(dim, bands[:, kcol0 + 0], ky=bands[:, kcol0 + 1], use_mesh_grid=use_mesh_grid, nkpt=nkpt)
    elif dim == 3:   
        (kx,ky,kz,dim) = refine_dims_and_ks(dim, bands[:, kcol0 + 0], ky=bands[:, kcol0 + 1],kz=bands[:, kcol0 + 2], use_mesh_grid=use_mesh_grid, nkpt=nkpt)  

    
    if use_mesh_grid:
        E = get_Ei_from_bands (kx,ky,kz, bands[:,kcol0:kcol1], nbands, bands[:,ecol], dim)
    else:
        E = get_Eijm_from_bands (kx,ky,kz, bands[:,kcol0:kcol1], nbands, bands[:,ecol], dim)
    
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
                    for io in range(E.shape[3]):
                        for ie in range(E.shape[2]):
                            out.write("\t%e"%(E[i,j,ie,io]))
                    out.write("\n")
            else:
                out.write("%e"%(kk[i,0]))
                for io in range(E.shape[3]):
                    for ie in range(E.shape[2]):
                        out.write("\t%e"%(E[i,0,ie,io]))
                out.write("\n")
            
            
            if nk[1] > 1 : out.write("\n") 
            
        out.write("\n")             
        out.close()

def write_gpl_meshgrid(fname, E, kk, header = None, append = False, kk_reduced = None):
        if append:
            out = open( fname+".gpl","a")
        else:
            out = open( fname+".gpl","w")
                        
        if header is not None:
            out.write("%s"%header)

        # print kk

        if kk_reduced is not None:
            idx= np.lexsort((kk_reduced[:,1],kk_reduced[:,0]))
        else:
            idx= np.lexsort((kk[:,1],kk[:,0]))

        # print kk[idx,0:2]

        vprev = np.zeros(2) 
        v     = np.zeros(2) 
        vprev_r = np.zeros(2) 
        v_r     = np.zeros(2) 
        for ik in range(E.shape[0]):
            v = kk[idx[ik],0:2]
            
            if kk_reduced is not None:
                v_r = kk_reduced[idx[ik],0:2]
                if v_r[0] != vprev_r[0]: out.write("\n")
                vprev_r = v_r
            else:
                if v[0] != vprev[0]: out.write("\n")
                vprev = v
                
                
            out.write("%e\t%e"%(v[0], v[1]))
            for ie in range(E.shape[1]):
                for io in range(E.shape[2]):
                    out.write("\t%e"%(E[idx[ik],ie,io]))
                    
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


def generate_kpath(pts, nk, spacing = None):
    # nk = 100 if nk is None else nk

    length = np.zeros(pts.shape[0]-1)
    nki = np.zeros(pts.shape[0]-1)
    

    for i in range (pts.shape[0]-1):
        p1 = pts[i,:]
        p2 = pts[i+1,:]
                
        v = (p2-p1)
        length[i] = (np.sqrt(np.dot(v,v)))

    du = length.sum()/nk  if spacing is None else spacing


    kpath = []
    upath = []
    for i in range (pts.shape[0]-1):
        p1 = pts[i,:]
        p2 = pts[i+1,:]
        
        v = (p2-p1)
        length = (np.sqrt(np.dot(v,v)))

        u = np.linspace(0, 1, num = length/du)
    
        # Generate the grid-points on the line 
        kk = np.zeros([u.shape[0],3])
        for iu in range(u.shape[0]):
            kk[iu,:] = np.array(line(v, p1 ,u[iu])[:])
        kpath.append(kk)
        upath.append(u)    
            
    
    return kpath, upath
    

def slice_on_line(E, kmesh, dim, nk, p1, p2, len0, spacing = None, kkin = None, uin = None):
    # We assume in the following that the kspace is 2D
        
    nbands = E.shape[3]

    #Same spacing as the original grid unless otherwise specified
    # dk = np.zeros(dim-1)
    # for idim in range(dim-1):
    #     dk[idim] = np.abs(kmesh[idim][1]-kmesh[idim][0])
    #
    #
    #
    # v = (p2-p1)
    # length = (np.sqrt(np.dot(v,v)))
    # du = min(dk[:])
    # if spacing:
    #     du = spacing
    # u = np.linspace(0, 1, num = length/du)
    # # print u
    #
    # Eout = np.zeros([u.shape[0],1, nbands, E.shape[4]])
    #
    # # Generate the grid-points on the line
    # kk = np.zeros([u.shape[0],dim])
    # for iu in range(u.shape[0]):
    #     kk[iu,:] = np.array(line(v, p1 ,u[iu])[0:dim])

    ########
    

    u = uin 
    Eout = np.zeros([u.shape[0],1, nbands, E.shape[4]])
    kk = np.zeros([u.shape[0],dim])
    kk[:,0:dim] = kkin[:,0:dim]
    
    v = kk[-1,:]-kk[0,:]
    length = (np.sqrt(np.dot(v,v)))
    

    # Check that the requested points are within the bounds of the input grid
    if (not points_in_bounds(kk,grid_bounds(kmesh,dim))):
        print('Error:Some of the requested output grid point is outside the bounds of the input grid.\nCurrent bounds are:\n%s'%grid_bounds(kmesh, dim))
        sys.exit(1)

    #Loop over all the bands
    for ib in range(nbands):
        for io in range(E.shape[4]):
            if dim < 3:
                f = interpolate.RectBivariateSpline(kmesh[0][:], kmesh[1][:], E[:,:,0,ib,io])
                Eout[:,0,ib,io] = f.ev(kk[:,0],kk[:,1])
            else:
                # print (kmesh[0][:], kmesh[1][:], kmesh[2][:])
                # print kk
            
             
                f = interpolate.RegularGridInterpolator((kmesh[0][:], kmesh[1][:], kmesh[2][:]), E[:,:,:,ib,io])
                Eout[:,0,ib,io] = f(kk)
    
    len1 = len0 + length
    u = u * length + len0

    return (Eout, u, len1)      
        

def import_file(file, abs_coords, use_mesh_grid = False):
    
    if "bands-gp.dat" in file.lower():
        imported = import_bands_file(file, abs_coords=abs_coords, use_mesh_grid=use_mesh_grid)
        
    elif "eigenvalues" in file.lower() :     
        imported = import_eigenvalues_file(file, abs_coords=abs_coords, use_mesh_grid=use_mesh_grid)            
    elif "info" in file.lower():     
        imported = import_eigenvalues_file(file, abs_coords=abs_coords, use_mesh_grid=use_mesh_grid)            
    else:
        try: 
            imported = import_bands_file(file,not args.absolute)
        except:
            print("""
Error: Unrecognized input file %s. The supported files 
are \'bands-gp.dat\', \'info\' and \'eigenvalues\'."""%(file))
            sys.exit(1)
    
    return imported
    


######################## MAIN ################################

def main(args):

    VERSION = '0.1'
    
    desc="""This utility converts the bands structure file generated by Octopus 
into a human-readable and easy-to-plot format.  
The supported files include 'bands-gp.dat', 'eigenvalues', and 'info'.
"""

    epilog="""Examples:

To process input file(s):
\'%(prog)s file1 [file2 ...]\'
\n\n
    """

    # parser = argparse.ArgumentParser(version='%s version %s' %(sys.argv[0],VERSION),
    #                                  description=desc,
    #                                  epilog=epilog,
    #                                  formatter_class=argparse.RawTextHelpFormatter)
    parser = argparse.ArgumentParser(description=desc,
                                     epilog=epilog,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # parser.add_argument(action='version', version='%s version %s' %(sys.argv[0],VERSION))

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
between [-1/2,1,2]. 
For example:

-c "p1x,p1y,p1z;p1x,p1y,p1z;[...]" 

Point coordinates can also be expressed as mathematical 
formula such as:

-c "0, 0, 0 ; 1/2, cos(pi/4)/2, sin(pi/4)/4; 0.5, 0.5, 0.5"

optionally the string can be stored on a file and given as 
option:

-c path_file 

""")                      

    parser.add_argument('--spacing', type=float, metavar='float', default= None, 
    help=
"""Define the grid spacing on the cut line. The default value is 
the minimum spacing of the input data.
"""    
    )

    parser.add_argument('--npoints', type=int, metavar='int', default= 100, 
    help="Number of points on the cut line. By default npoints = 100."    
    )

    
    parser.add_argument('--absolute', action="store_true", default=False, 
    help="Use absolute coordinates in reciprocal space.")
    
    parser.add_argument('--outpoints', action="store_true", default=False, 
    help=
"""Prints the cut line as a zero-weight k-point path suitable 
for octopus input file. 
Note: positional argument 'file' is ignored with this option."""    
    )
    

    # parser.add_argument('file', nargs='+')
    # parser.add_argument('file', nargs='?')
    parser.add_argument('file', nargs='*')
    

    args = parser.parse_args()
    
    
    # parse the cut nodal points
    if args.cut is not None:
        try:
            # read from file
            f = open(args.cut)
            path_string = f.readline()
            f.close()
        except:
            # read from command line  
            path_string = args.cut  
        tokens = np.array(path_string.split(';'))
        pts = np.zeros((tokens.shape[0],3))
        for i in range(tokens.shape[0]):
            coords  = tokens[i].split(',')
            for j in range(len(coords)):
               pts[i,j] = eval(coords[j])  
               
        kpath, upath = generate_kpath(pts, nk = args.npoints, spacing = args.spacing)       
        
        # print pts.shape[:]
        # dir = np.array(args.cut.split(','),dtype=float)    
    
    # Write the points along the path
    if args.outpoints:
        if args.cut is None:
            raise Exception("You have to specify a cut path.")
        
        # print kpath
        for kk in enumerate(kpath):
            # print kk[1]
            # for i in range(kk[1].shape[0]):
            #     print ("%2.6f | %2.6f | %2.6f | %2.6f "%(0.0, kk[i,0], kk[i,1], kk[i,2]))
            for i in range(kk[1].shape[0]):
            # for k in enumerate(kk[1]):
                k =  kk[1][i]
                print ("%2.6f | %2.6f | %2.6f | %2.6f "%(0.0, k[0], k[1], k[2]))
        
        exit(0)    
    
    
    # check if there are input files 
    if args.file is None or len(args.file) == 0:
        parser.print_usage()
        print "%s: error: too few arguments."%(os.path.basename(sys.argv[0]))
        print "You must specify a file if not using '--outpoints'."        
        exit(1)
    
    # loop over all the input files
    for file in args.file: 

        if have_kpoint_symmetries(file):
            warnings.warn(
""""
 It seems that the k-point grid has been generated using symmetries.
 Since symmetries are not supported by this program the results are likely to
 be garbage.
"""
            )
        if args.absolute:
            imported = import_file(file, abs_coords = args.absolute)
            kmesh_abs = imported[1]            
            imported = import_file(file, abs_coords = False, use_mesh_grid = True)
        else:        
            imported = import_file(file, abs_coords=args.absolute)
        
#         if "bands-gp.dat" in file.lower():
#             if args.absolute:
#                 imported = import_bands_file(file, abs_coords = args.absolute)
#                 kmesh_abs = imported[1]
#                 imported = import_bands_file(file, abs_coords = False, use_mesh_grid = True)
#             else:
#                 imported = import_bands_file(file, abs_coords = args.absolute)
#
#         elif "eigenvalues" in file.lower() :
#             imported = import_eigenvalues_file(file)
#         elif "info" in file.lower():
#             imported = import_eigenvalues_file(file)
#         else:
#             try:
#                 imported = import_bands_file(file,not args.absolute)
#             except:
#                 print("""
# Error: Unrecognized input file %s. The supported files
# are \'bands-gp.dat\', \'info\' and \'eigenvalues\'."""%(file))
#                 sys.exit(1)

        E = imported[0]  
        kmesh = imported[1]            
        dim = imported[2] 
        meshdim = len(E.shape )-1
        
        
        if not args.absolute:
            print "Detected a %d-dimensional k-space with %d (%d x %d x %d) k-points and %d bands "%(dim,np.prod(E.shape[0:2]),E.shape[0],E.shape[1],E.shape[2],  E.shape[3])
            header = "# This file contains %s bands"%(E.shape[3])

            if E.shape[4]>1:
                print "with %d spin components"%(E.shape[4]-1)
                header = "%s with %d spin components (reduced coordinates).\n"%(header,E.shape[4]-1)
            else:
                header = "%s.\n"%(header)
        else:
            print "Detected a %d-dimensional k-space with %d k-points and %d bands "%(dim,E.shape[0],E.shape[1])
            header = "# This file contains %s bands"%(E.shape[1])
            if E.shape[2]>1:
                print "with %d spin components"%(E.shape[2]-1)
                header = "%s with %d spin components (absolute coordinates).\n"%(header,E.shape[2]-1)
            else:
                header = "%s.\n"%(header)
            
                 
            
        nk = np.zeros(3)
        nk = [kmesh[0].shape[0],kmesh[1].shape[0],kmesh[2].shape[0]]

        
                
            
        if args.cut is not None:
            header = "%s# Bands are evaluated along the following path in the Brillouin zone:\n# %s\n"%(header,path_string.strip("\n"))
            file = file+".cut"
            append = False
            lenght = 0.0
            for i in range(pts.shape[0]-1):
                 p1 = pts[i,:]     
                 p2 = pts[i+1,:]                      
                 print "Segment: %s --> %s"%(p1[0:dim],p2[0:dim])
                 header = "%s#\n# Slice on a line segment connecting %s and %s \n#\n"%(header,p1[0:dim],p2[0:dim])
                 (E_,kx_,lenght) = slice_on_line(E, kmesh, dim, nk, p1, p2, lenght, spacing = args.spacing, kkin = kpath[i], uin = upath[i])
                 nk_ = np.array((kx_.shape[0],0,0))
                 kk = np.zeros([max(nk_[:]),3])
                 kk[0:nk_[0],0] = kx_[0:nk_[0]]
                 kk[0:nk_[1],1] = kmesh[1][0:nk_[1]]
                 
        
                 write_gpl(file, E_, kk, nk_, header, append = append)
                 header = ""
                 append = True
        elif dim <= 2:
            if not args.absolute :
                E_ = E[:,:,0,:,:]        
                kk = np.zeros([max(nk[:]),3])
                kk[0:nk[0],0] = kmesh[0][0:nk[0]]
                kk[0:nk[1],1] = kmesh[1][0:nk[1]]
                
                write_gpl(file, E_, kk, nk, header)
            else:
                E_ = E[:,:,:]        
                kk_red = np.zeros([kmesh[0].shape[0],3])
                kk_red[0:nk[0],0] = kmesh[0][0:nk[0]]
                kk_red[0:nk[1],1] = kmesh[1][0:nk[1]]

                kk_abs = np.zeros([kmesh[0].shape[0],3])
                kk_abs[0:nk[0],0] = kmesh_abs[0][0:nk[0]]
                kk_abs[0:nk[1],1] = kmesh_abs[1][0:nk[1]]

                
                write_gpl_meshgrid(file, E_, kk_abs, header, kk_reduced = kk_red)

        
        
        if len(imported) > 3:
            # We have a zero-weight path
            file = file+".cut0"

            E0     = imported[3]
            kmesh0 = imported[4]
            nbands = imported[5]
                        
            nk_ = np.array([int(kmesh0.shape[0]/nbands), 0, 0])
            kk = np.zeros([nk_[0],3])
            E_ = np.zeros([nk_[0],1, nbands, E0.shape[1]])
            
            def ll(i,j):
                return i*nbands + j
            
            for i in range(nk_[0]):
                for j in range(nbands):
                    l = i*nbands + j
                    E_[i,0,j,:] = E0[ll(i,j),:]    
                if i > 0:                    
                    dk = np.sqrt(np.sum((kmesh0[ll(i,j),:]-kmesh0[ll(i-1,j),:])**2))
                    kk[i] = kk[i-1] + dk

            
            
            print "Detected a zero-weight path (%d points) in k-space with %d bands"%(kk.shape[0], nbands )
            header = "# Zero-weight path in k-space (%d points) with %d bands\n"%(kk.shape[0], nbands)
            
            write_gpl(file,E_, kk, nk_, header)

    


if __name__ == "__main__":
   main(sys.argv[1:])