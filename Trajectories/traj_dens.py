import numpy as np
import traj_lib as tlib

dim  = 1

data=np.loadtxt('td.general/trajectories')
ntr = int((data.shape[1]-2)/dim)
ntr = ntr/2
traj = data[:, 2:ntr+2]
Jacobian = data[:, ntr+2:]
itern = data[:,0]

rho=np.loadtxt('static/density.y=0,z=0')
points = traj[0, :]
weights = tlib.get_weight(points,rho)

dt = abs(data[1,1]- data[0,1])

itt = 15000
it = np.where(itern == itt)[0][0]


dx = abs(rho[1,0]-rho[0,0])
hist = np.zeros(rho.shape[0])
npt = rho.shape[0]
for itr in range(ntr):
    ip = int((traj[it,itr]- rho[0,0])/dx)
    if (ip < npt):
        hist[ip] +=weights[itr]*vol0[itr]/dx

#polinomial least square fitting 
polfit = np.poly1d(np.polyfit(traj[it,:], weights[:],40, w= vol0[:]))     
rhopf = polfit(rho[:,0])
    
print hist.sum()*dx    
idx = np.where((rho[:,0] > -5)  & ( rho[:,0]<5))
print rho[idx,1].sum()*dx   

f = open('rho.trj','w')
for ii in range(hist.shape[0]):
   f.write("%1.3e\t %1.3e\t %1.3e\n"%(rho[ii,0],hist[ii], rhopf[ii]))       
   
f.close()      

f = open('rho.trj.raw','w')
idx=np.argsort(traj[it,:])
for itr in range(ntr):
   rhotrj = weights[idx[itr]]/Jacobian[it,idx[itr]] 
   f.write("%1.3e\t %1.3e\n"%(traj[it,idx[itr]],rhotrj))        
f.close()      

