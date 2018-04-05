from __future__ import division
import numpy as np
import pandas as pd
import traj_weigth as weight


############
Radius = 50

ne = 200
kmax = 2
Emax = kmax**2/2.

selectV = [-0.12, -0.05]#[0.05, 0.12]

have_vel = False
#############


def get_pes(traj, weights, laser, dt, v_traj=None, vtarget = None):
    """Photoelectron velocity distribution"""
    use_vel = False
    if v_traj is not None:
        use_vel = True
    print use_vel
    
    ntr = traj.shape[1]

    velocities = []
    ww = []    
    targ_trj  =[]
    targ_vtrj =[]
    
    nocross = 0

    if False:
        Rthr= abs(rho[1,0]-rho[0,0])*3 #grid spacing 
        for itr in range(ntr):
            trj  = traj[:,itr]
            trjR = trj[:]**2
            # print trjR
            # idx  = np.where( abs(trjR - Radius**2)<=Rthr**2 )
            idx  = np.where( trjR > Radius**2)
            # Skip the trajectories that don't cross the border
            if not len(idx[0]):
                nocross += 1
                continue
            idx=idx[0]
        
            #identify the crossing time indices idcr 
            edges=np.roll(idx,-1)-idx
            idcr = np.where( edges > 1)
        
            if not len(idcr[0]):
                idx=np.array([idx[0]])
            else:
                # print itr,idx
                # print len(idx)
                # print edges
                idcr=np.array(idcr)
                # print idcr,idcr+1
                tmp = np.array([idx[idcr],idx[idcr+1]]).flatten()
                # print tmp
                idx= np.append(np.array(idx[0]),tmp).flatten()
                # print idx
          
            if len(idx)>1:    
                print ">>>>>>>",idx,len(idx)

            vv=np.zeros(dim) 
            for it in range(len(idx)):
                it = idx[it]

                # vv +=  (trj[it]- trj[it-1])/dt   + laser[it]/137
                if (use_vel):
                    vv = v_traj[it,itr] + laser[it]/137
                else:
                    vv =  (3/2*trj[it] -2*trj[it-1] +1/2*trj[it-2])/dt + laser[it]/137#/weights[itr]
                
                if (laser[it]/137 > 0):
                    print itr, it, trj[it]
                    print "vv=",vv, "laser = ", laser[it]/137

                if vtarget is not None:
                    if (vv <= vtarget[1] and vv >= vtarget[0]):
                        targ_trj.append(trj)
                    
                
            velocities.append(vv)
            ww.append(weights[itr])
        

    else:
    #analyze the current outside the region r>R at max time
        nabegde=0 
        maxtrj = 0 
        it = traj.shape[0]-1
        # it = 7000
        for itr in range(ntr):
            trj  = traj[:,itr]
            trjR = trj[:]**2
            maxtrj = np.sqrt(trjR[it]) if (np.sqrt(trjR[it])>maxtrj) else maxtrj 
            if trjR[it-1]>=Radius**2:
                if (use_vel):
                    vv = v_traj[it,itr]
                else:
                    vv =  (trj[it]- trj[it-1])/dt #- laser[it]/137
            
                # print vv - (trj[it]- trj[it-1])/dt
                if vtarget is not None:
                    if (vv <= vtarget[1] and vv >= vtarget[0]):
                        targ_trj.append(trj)                
   
                velocities.append(vv)
                ww.append(weights[itr])
            else:
                nocross += 1
    

            if trjR[it-1]>=170**2:
                nabegde+=1
        print "maxtrj =", maxtrj        
        print "number of trajectories in the boundary region =",nabegde,"(%2.2f)"%(nabegde*100.0/ntr)


    print "number of trajectories not crossing the border =",nocross,"(%2.2f)"%(nocross*100.0/ntr)


    velocities = np.array(velocities)
    ww = np.array(ww)
    targ_trj = np.array(targ_trj)
    
    return (velocities, ww, targ_trj)

###################################
###################################
###################################
if __name__ == "__main__":

    rho=np.loadtxt('static/density.y=0,z=0')
    dim = rho.shape[1]-1

    laserd=np.loadtxt('td.general/laser')
    laser=np.zeros((laserd.shape[0],dim))
    # get the full field
    for idim in range(dim):
        laser[:,idim] =laserd[:,idim+2::dim].sum(axis=1)

    # data = np.array(pd.read_csv('td.general/trajectories', comment="#",delim_whitespace=True))
 
    data=np.loadtxt('td.general/trajectories')
    ntr = int((data.shape[1]-2)/dim)
    if (have_vel):
        ntr = int(ntr/2) 
        traj = data[:, 2:ntr+2]
        v_traj = data[:, ntr+2:]
    else:
        traj = data[:, 2:]
        v_traj = None

    points = data[0, 2:]
    weights = weight.get_weight(points,rho)

    dt = abs(data[1,1]- data[0,1])


    velocities, ww, outtrj = get_pes(traj, weights, laser, dt, v_traj = v_traj, vtarget = selectV)


    # Energy distribution
    
    kinen = velocities**2/2.0 #electron mass is 1 in atomic units
    Egrid=np.linspace(0., Emax,ne)
    De = Egrid[1]-Egrid[0]

    hist = np.zeros(Egrid.shape[0])
    pes = np.zeros(Egrid.shape[0])

    #binning 
    for ii in range(kinen.shape[0]):
        idx = int((kinen[ii]-Egrid.min())/De)
        if (idx <= ne):
            hist[idx] += 1 
            pes[idx] += ww[ii]

    print "sum =", pes.sum()

    f = open('espect','w')
    for ii in range(hist.shape[0]):
       f.write("%1.3e\t %1.3e\t %1.3e\n"%(Egrid[ii],pes[ii],hist[ii]))       
   
    f.close()   

    # Velocity distribution

    Vgrid=np.linspace(-kmax, kmax,ne)
    Dv = Vgrid[1]-Vgrid[0]

    hist = np.zeros(Egrid.shape[0])
    pes = np.zeros(Egrid.shape[0])

    #binning 
    for ii in range(velocities.shape[0]):
        idx = int((velocities[ii] -Vgrid.min())/Dv)
        if (idx <= ne):    
            hist[idx] += 1 
            pes[idx] += ww[ii]

    print "sum =", pes.sum()

    f = open('vspect','w')
    for ii in range(hist.shape[0]):
       f.write("%1.3e\t %1.3e\t %1.3e\n"%(Vgrid[ii],pes[ii],hist[ii]))       
   
    f.close()   
    
    # Selected trajectories
    f = open('trj.sel','w')
    f.write("##############\n#\n")
    f.write("# Trajectories number = %d\n"%(outtrj.shape[0]))
    f.write("# Vrange = %s\n"%(selectV))
    f.write("#\n##############\n")
    
    for ii in range(data.shape[0]):
       f.write("%d\t %1.3e"%(data[ii, 0],data[ii, 1]))
       for itr in range(outtrj.shape[0]):
           f.write("\t %1.3e"%(outtrj[itr,ii]))      
       f.write("\n")
   
    f.close()   

    
    