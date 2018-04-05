import numpy as np
from traj_weigth import *



data=np.loadtxt('td.general/trajectories')
traj = data[:, 2:]

rho=np.loadtxt('static/density.y=0,z=0')
points = data[0, 2:]
weights = get_weight(points,rho)

ntr = (data.shape[1]-2)
dt = abs(data[1,1]- data[0,1])

it = 6000

dx = abs(rho[1,0]-rho[0,0])
hist = np.zeros(rho.shape[0])
print hist.shape, dx
for itr in range(ntr):
    ip = int(np.rint((traj[it,itr]- rho[0,0])/dx))
    hist[ip] +=weights[itr]/dx
    
    
print hist.sum()    
f = open('rho.trj','w')
for ii in range(hist.shape[0]):
   f.write("%1.3e\t %1.3e\n"%(rho[ii,0],hist[ii]))       
   
f.close()      