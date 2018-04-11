from __future__ import division
import numpy as np
import pandas as pd
import traj_weigth as weight

import sys
import yaml

import time
import pprint


#---------------------------------------
def init_params(parameters_in = None, append_to_default = None):
    """docstring for init_params"""

    ########## DEFAULT PARAMETERS #########

    parameters  = {}
    parameters['radius']   = 50
    parameters['have_vel'] = False
    parameters['strategy'] = 'full'
    
    parameters['density_file']    = 'static/density.y=0,z=0'
    parameters['trajectory_file'] = 'td.general/trajectories'
    parameters['laser_file']      = 'td.general/laser'


    # output
    parameters['pmax']   = 2
    parameters['emax']   = parameters['pmax']**2/2.
    parameters['nume']   = 200

    parameters['selectV'] = None #[-0.12, -0.05]
    

    if append_to_default is not None:
        parameters.update(append_to_default)
    
    ##### input parameters ######
    if parameters_in is not None:
        for key in parameters_in:
            if key.lower() in parameters:
                parameters[key.lower()] = parameters_in[key]

    pprint.pprint('#####################################################')
    pprint.pprint(parameters)
    pprint.pprint('#####################################################')


    return parameters

#---------------------------------------
def parse_inp(inp_file):
    
    with open(inp_file, 'r') as stream:
        try:
            params= yaml.load(stream)
            for key in params:
                try:
                    params[key] = eval(params[key])

                except:
                    pass
                            
        except yaml.YAMLError as exc:
            print(exc)
            params = {}

    return params        


#---------------------------------------
def get_pes(parameters, traj, weights, laser, dt, v_traj=None, vtarget = None):
    """Photoelectron velocity distribution"""
    use_vel = False
    if v_traj is not None:
        use_vel = True
    
    Radius = parameters['radius']
    dim = parameters['dim']
    
    ntr = traj.shape[1]

    velocities = []
    ww = []    
    targ_trj  =[]
    targ_vtrj =[]
    
    nocross = 0

    if parameters['strategy'] == 'flux':
        Rthr= abs(rho[1,0]-rho[0,0])*3 #grid spacing 
        for itr in range(ntr):
            trj  = traj[:,itr,:]
            trjR = np.sum(trj[:,:]**2, axis=1)
        
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

                if (use_vel):
                    vv = v_traj[it,itr] + laser[it]/137
                else:
                    vv =  (3/2*trj[it] -2*trj[it-1] +1/2*trj[it-2])/dt + laser[it]/137
                
                if (laser[it]/137 > 0):
                    print itr, it, trj[it]
                    print "vv=",vv, "laser = ", laser[it]/137

                if vtarget is not None:
                    if (vv <= vtarget[1] and vv >= vtarget[0]):
                        targ_trj.append(trj)
                    
                
            velocities.append(vv)
            ww.append(weights[itr])
        

    elif parameters['strategy'] == 'full':
    #analyze the current outside the region r>R at max time
        nabegde=0 
        maxtrj = 0 
        it = traj.shape[0]-1
        # it = 7000
        for itr in range(ntr):
            trj  = traj[:,itr,:]
            trjR = np.sum(trj[:,:]**2, axis=1)
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
    

        print "trajectories maximum value =", maxtrj        

    else:
        raise NameError('Undefined strategy \''+ parameters['strategy']+'\'')

    print "number of trajectories not crossing the border =",nocross,"(%2.1f%%)"%(nocross*100.0/ntr)


    velocities = np.array(velocities)
    ww = np.array(ww)
    targ_trj = np.array(targ_trj)
    
    return (velocities, ww, targ_trj)

###################################
###################################
###################################
if __name__ == "__main__":


    if len(sys.argv) > 1:
        inp_file = sys.argv[1]

        params_in = parse_inp(inp_file)
    
    else:
        params_in={}
    
    # oct_parameters = parse_inp('exec/parser.log')
    # pprint.pprint(oct_parameters)
    # exit()
    
    parameters = init_params(params_in)


    selectV = parameters['selectV']
    have_vel = parameters['have_vel']

    rho=np.loadtxt(parameters['density_file'])
    dim = rho.shape[1]-1
    parameters['dim'] = dim
    print "Dimensions: %d"%(dim)

    laserd=np.loadtxt(parameters['laser_file'])
    laser=np.zeros((laserd.shape[0],dim))
    # get the full field
    for idim in range(dim):
        laser[:,idim] =laserd[:,idim+2::dim].sum(axis=1)

    # data = np.array(pd.read_csv('td.general/trajectories', comment="#",delim_whitespace=True))
 
    data=np.loadtxt(parameters['trajectory_file'])
    ntr = int((data.shape[1]-2)/dim)
    if (have_vel):
        ntr = int(ntr/2) 
        traj = data[:, 2:ntr+2]
        v_traj = data[:, ntr+2:]
        v_traj = v_traj.reshape((-1,ntr,dim))
    else:
        traj = data[:, 2:]
        v_traj = None
        
    traj = traj.reshape((-1,ntr,dim))
    
    print "Number of trajectories = %d"%(ntr)

    points = data[0, 2:]
    weights = weight.get_weight(points,rho)

    dt = abs(data[1,1]- data[0,1])


    velocities, ww, outtrj = get_pes(parameters, traj, weights, laser, dt, v_traj = v_traj, vtarget = selectV)


    #OutPut 

    # Energy distribution
    ne   = parameters['nume']
    Emax = parameters['emax']
    
    
    kinen = np.sum(velocities[:,:]**2,axis=1)/2.0 #electron mass is 1 in atomic units
    Egrid=np.linspace(0., Emax, ne)
    De = Egrid[1]-Egrid[0]

    hist = np.zeros(Egrid.shape[0])
    pes = np.zeros(Egrid.shape[0])

    #binning 
    for ii in range(kinen.shape[0]):
        idx = int((kinen[ii]-Egrid.min())/De)
        if (idx <= ne):
            hist[idx] += 1 
            pes[idx] += ww[ii]

    print "Pes spectrum integral =", pes.sum()*De

    f = open('espect.'+parameters['strategy'],'w')
    for ii in range(hist.shape[0]):
       f.write("%1.6e\t %1.6e\t %1.6e\n"%(Egrid[ii],pes[ii]*np.sqrt(2*Egrid[ii]),hist[ii]))       
   
    f.close()   

    # Velocity distribution
    kmax = parameters['pmax']
    Lgrid=np.linspace(-kmax, kmax,ne)
    if dim ==1:
        Vgrid = np.meshgrid(Lgrid)    
    elif dim ==2:
        Vgrid = np.meshgrid(Lgrid,Lgrid, indexing='ij')
    elif dim ==3:
        Vgrid = np.meshgrid(Lgrid,Lgrid,Lgrid,  indexing='ij')


    DL = Lgrid[1]-Lgrid[0]


    hist = np.zeros(Vgrid[0].shape[:])
    pes = np.zeros(Vgrid[0].shape[:])


    #binning 
    idxa=np.zeros(dim,dtype=int)
    for ii in range(velocities.shape[0]):
        for idim in range(dim):
            idxa[idim] = int((velocities[ii, idim] -Lgrid.min())/DL)    
        if (all(idxa[:] <= ne)):  
            idx = tuple(idxa)  
            # diff = velocities[ii,:] -[Vgrid[0][idx], Vgrid[1][idx]]
            # if (any(diff > DL)):
            #     print idx
            #     print dist, DL*np.sqrt(2),DL
            #     print diff
            #     print velocities[ii,:], Vgrid[0][idx], Vgrid[1][idx]
            hist[idx] += 1 
            pes[idx] += ww[ii]


    print "sum =", pes.sum()

    f = open('vspect.'+parameters['strategy'],'w')
    if   (dim==1):
        for ii in range(Lgrid.shape[0]):
                f.write("%1.6e\t"%(Vgrid[0][ii]))
                f.write("%1.6e\t %1.6e\n"%(pes[ii],hist[ii]))
            
    elif (dim==2):        
        for jj in range(Lgrid.shape[0]):
            for ii in range(Lgrid.shape[0]):
                for jdim in range(dim):
                    f.write("%1.6e\t"%(Vgrid[jdim][ii,jj]))
                    # f.write("%1.3e\t %1.3e\t %1.3e\n"%(Vgrid[ii],pes[ii],hist[ii]))
                f.write("%1.6e\t %1.6e\n"%(pes[ii,jj],hist[ii,jj]))
            f.write("\n")
   
    f.close()   

    # debug the velocity distribution
    f = open('vspect.raw','w')
    for ii in range(velocities.shape[0]):
        for idim in range(dim):
            f.write("%1.6e\t"%(velocities[ii, idim]))
        f.write("\n")    
    f.close() 
    
    # Selected trajectories
    if selectV is not None:
        f = open('trj.sel','w')
        f.write("##############\n#\n")
        f.write("# Trajectories number = %d\n"%(outtrj.shape[0]))
        f.write("# Vrange = %s\n"%(selectV))
        f.write("#\n##############\n")
    
        for ii in range(data.shape[0]):
           f.write("%d\t %1.6e"%(data[ii, 0],data[ii, 1]))
           for itr in range(outtrj.shape[0]):
               f.write("\t %1.6e"%(outtrj[itr,ii]))      
           f.write("\n")
   
        f.close()   

    
    