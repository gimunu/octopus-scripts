import numpy as np
from scipy.spatial import Voronoi


# dim = 1


def get_weight(points, rho):
    """docstring for get_weight"""


    dim = rho.shape[1]-1
    if dim == 1:
        (ww,vol) = get_weight_1D(points, rho)
    else: 
        # need to implement properly
        ww = np.zeros(points.shape[0])
        ww = ww + rho[:,dim].sum()*(rho[1,dim-1]-rho[0,dim-1])**dim/points.shape[0]
    #
    # if dim==2:
    #     vor = Voronoi(points)
    #     print vor


    sumW = ww.sum()
    sumN = rho[:,dim].sum()*(rho[1,dim-1]-rho[0,dim-1])**dim
    # Checksum 
    if abs(sumW-sumN)> 1E-4:
        print "Weights density checksum failed: %e instead of %e"%(sumW, sumN)

    
    
    return (ww,vol)
    


def get_weight_1D(points, rho):
    """Integrate on 1D Voronoi cells"""
    
    Lmax = rho[:,0].max()
    Lmin = rho[:,0].min()
    
    idx=np.argsort(points)
    idx_inv=np.argsort(idx)
    pts = points[idx]

    Dleft=(pts - np.roll(pts,1))/2.
    # Dleft[0] =pts[0] - Lmin
    Dleft[0] = 0

    Dright=(np.roll(pts,-1)-pts)/2.
    # Dright[-1] =Lmax - pts[-1]
    Dright[-1] =0
        
    weight = np.zeros(pts.shape[0])
    
    # for ip in range(pts.shape[0]):
    #     # the cells at the extreme are usually large and
    #     # need more points for the integral
    #     if ip == 0 or ip == pts.shape[0]-1:
    #         nsamples = 1E5
    #     else:
    #         nsamples = 1e3
    #     x1 = pts[ip]-Dleft[ip]
    #     x2 = pts[ip]+Dright[ip]
    #     xx=np.linspace(x1,x2,nsamples, endpoint= False)
    #     dx = xx[1]-xx[0]
    #     weight[ip] = np.sum(np.interp(xx, rho[:,0], rho[:,1]))*dx

    vol = abs(Dright[:]+Dleft[:])
    for ip in range(pts.shape[0]):
        weight[ip] = np.interp(pts[ip], rho[:,0], rho[:,1])

    
    print np.sum(weight[1:-1]*vol[1:-1])
    print np.sum(weight[:]*vol[:])
    #come back to the original point ordering 
    return (weight[idx_inv], vol[idx_inv])

###################################
###################################
###################################
if __name__ == "__main__":
    data=np.loadtxt('td.general/trajectories')
    points = data[0, 2:]

    rho=np.loadtxt('static/density.y=0,z=0')

    ww = get_weight(points,rho)
    print ww    

    # print data.shape[:], points.shape[:]

    