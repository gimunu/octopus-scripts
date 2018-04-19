import numpy as np
import traj_weigth as weight



data=np.loadtxt('td.general/trajectories')
traj = data[:, 2:]

rho=np.loadtxt('static/density.y=0,z=0')
points = data[0, 2:]
(weights,vol0) = weight.get_weight(points,rho)

ntr = (data.shape[1]-2)
dt = abs(data[1,1]- data[0,1])

it = 500

points = data[it, 2:]
(trash,volt) = weight.get_weight(points,rho)

dx = abs(rho[1,0]-rho[0,0])
hist = np.zeros(rho.shape[0])
npt = rho.shape[0]
for itr in range(ntr):
    ip = int((traj[it,itr]- rho[0,0])/dx)
    if (ip < npt):
        hist[ip] +=weights[itr]*vol0[itr]/dx
    
    
print hist.sum()*dx    
idx = np.where((rho[:,0] > -5)  & ( rho[:,0]<5))
print rho[idx,1].sum()*dx   

f = open('rho.trj','w')
for ii in range(hist.shape[0]):
   f.write("%1.3e\t %1.3e\n"%(rho[ii,0],hist[ii]))       
   
f.close()      

f = open('rho.trj.raw','w')
idx=np.argsort(traj[it,:])
for itr in range(ntr):
   f.write("%1.3e\t %1.3e\n"%(traj[it,idx[itr]],weights[idx[itr]]*vol0[itr]/volt[itr]))        
f.close()      