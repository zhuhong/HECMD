#!/usr/bin/env python

'''
author: zhuh
date: 2013.10.30
version: 0.0.1
'''

import sys
import os
import re
import math
import getopt
import struct
import numpy as np
import MDAnalysis
import MDAnalysis.core.qcprot as qcp
import time as Time

K= 1.3806504e-23
'''It's the Boltzmann constant'''
HBA= 1.05457148e-34
'''It's the Planck constant'''
H = 6.62606957e-34
RC= 8.3145
'''It's the idea gas constant'''
E= 2.71828183
'''natural exp poment e'''
MP=1.66055402e-27
'''1 a.m.u. = Proton mass [kg] '''
NA=6.0221415e23        
''' Avogadro's Number  [1/gmol] '''

#CYCLES = 20


class Index_class:
    '''
    Core class
    '''
    group_name = " "
    group_list = []
    def __init__(self , group_name = " " , group_list = []):
        ''' Initialzation. '''
        self.group_name = group_name
        self.group_list = group_list
    def __str__(self):
        return "Group %15s has %6d elements" %(self.group_name , len(self.group_list))
    def __repr__(self):
        return "Group %15s has %6d elements" %(self.group_name , len(self.group_list))
        


def Read_index_to_Inclass(filename):
    '''Read in a index.ndx file ,  and return a list of Index_class'''
    try:
        fp = open(filename , 'r')
    except:
        sys.exc_info()[1]      

    list=[]
    temp_Index_class=Index_class()
    for line in fp:
        if ("[" in line):
       #     print line
            if  len(temp_Index_class.group_list)!=0:
              #  temp=Index_class()
                temp=temp_Index_class
                list.append(temp)
                temp_Index_class=Index_class()
            temp_Index_class.group_name=(line.split(" "))[1]
            temp_Index_class.group_list=[]
        elif len(line)>1: 
            line_units= re.split("\s+",line.strip())
            for s in line_units:
                temp_Index_class.group_list.append(int(s))
        else:
            pass
    list.append(temp_Index_class)
    fp.close()
    return list


def Print_Index(index_list):
    '''Print the index list like Gromacs does.'''
    for i in range(len(index_list)):
        print "Group %4d (%15s) has %6d elements" %(i , index_list[i].group_name , len(index_list[i].group_list))


def QHE(traj_data2,group_mass,cycle_frame):
    traj_ave  = np.zeros((3*natoms),dtype=np.float32)
    covar_mat = np.zeros((3*natoms,3*natoms),dtype=np.float32)

    for i in range(cycle_frame):
        traj_ave += traj_data2[i]

    traj_ave = traj_ave / cycle_frame

    for i in range(cycle_frame):
        _delta_vect=traj_data2[i]-traj_ave
        _mass_delta_vect=_delta_vect*group_mass
        covar_mat += np.matrix(_mass_delta_vect).T * np.matrix(_mass_delta_vect)
        sys.stderr.write("\t Reading frame %8d\r" %(i+1))
        sys.stderr.flush()
        
    # covar_mat = covar_mat / ( nframes -1 )   # be careful!
    covar_mat = covar_mat / cycle_frame
     
    # TIME1=Time.time()
    # print "\t Time used %6.2f s" %(TIME1 - BEGIN_TIME)


    print "step 4: Diagnoalizing the covariance matrix, cycling %d" %cycle
    eigenvalues,eigenvectors = np.linalg.eig(covar_mat)
    #    print eigenvalues
    # np.savetxt("eigenvalues.dat",eigenvalues)
    # sys.exit()
    # TIME2=Time.time()    
    # print "\t Time used %6.2f s" %(TIME2 - TIME1)


    eigen_thr = 1e-5
    truncate_num = 0
    for i in range(3*natoms):
        if eigenvalues[i] < eigen_thr:
            truncate_num = i
            break
        elif eigenvalues[-1] > eigen_thr:
            truncate_num = 3*natoms
    print "\t Using %d eigenvalues to calculate entropy." %truncate_num
        
    
    print "step 5: Calculating the quasiharmonic entropy, cycling %d" %cycle
    global HBA
    global K
    global MP
    global RC
    alpha_const  =HBA*(1e10)/math.sqrt(K*temperature*MP)
    # eigval_class = eigenvalues[:truncate_num]
    # eigvec_class = eigenvectors[:truncate_num]
    alpha     =[alpha_const/math.sqrt(eigenvalues[i]) for i in range(truncate_num)]
    s_quantum =[alpha[i]/(math.exp(alpha[i])-1) - math.log(1-math.exp(-alpha[i])) for i in range(truncate_num)]
    # np.savetxt("quantum.dat",s_quantum)
    
    total_entropy=sum(s_quantum)*RC
    print "\t Entopy: %8.2f J/mol/K, %8.2f Kcal/mol." %(total_entropy,total_entropy/4186.8*temperature)
    fp.write("\t Entopy: %8.2f J/mol/K, %8.2f Kcal/mol.\n" %(total_entropy,total_entropy/4186.8*temperature))
    fp.flush()


######################
    print "step 6: Calculating the anhormonic entropy"
   
    _bp_temp = [eigenvectors[:,i]*group_mass for i in range(truncate_num)]   
    # print np.shape(_bp_temp)
    _delta_temp= np.zeros((3*natoms,nframes),dtype=np.float32)
    for i in range(nframes):
        for j in range(3*natoms):
            _delta_temp[j,i]=traj_data2[i][j]-traj_ave[j]
    bp=np.matrix(_bp_temp)*np.matrix(_delta_temp)
    for i in range(truncate_num):
        np.savetxt("dist_%d.dat" %(i+1),bp[i])
    # print np.shape(bp)
#    print bp[0]
#   print bp[10]
    anh=np.zeros([truncate_num],dtype=np.float32)
    _beta_const = 0.5*(1-math.log(H*HBA/K/temperature))
    BIN_NUM = 100
    rxy=np.zeros([truncate_num],dtype=np.float32)
    # print _beta_const
    for m in range(truncate_num):
        _hist,_bins=np.histogram(bp[m],bins=np.linspace(np.amin(bp[m]),np.amax(bp[m]),BIN_NUM+1))   
        _prob = [(_item+0.01)/sum([_hist[i]*(_bins[i+1]-_bins[i]) for i in range(len(_hist))]) for _item in _hist]
        _int_prob =sum([math.log(_prob[i]*1e10/math.sqrt(MP))*_prob[i]*(_bins[i+1]-_bins[i]) for i in range(len(_hist))]) #care the unit
        # print _prob

        # for f(a)
        
        _fa = [ math.sqrt(abs(-1*eigenvalues[m]*math.log(_prob_item*math.sqrt(2*math.pi*eigenvalues[m])))) for _prob_item in _prob ]
        _x = np.zeros([len(_prob)],dtype=np.float32)
        _y = np.zeros([len(_prob)],dtype=np.float32)
        for kk in range(len(_prob)):
            _x[kk] = 0.5*(_bins[kk]+_bins[kk+1])
            if _x[kk] < 0 :
                _y[kk] = -1* _fa[kk]
            else:
                _y[kk] = _fa[kk]
        _x_ave = sum(_x) / BIN_NUM
        _y_ave = sum(_y) / BIN_NUM
        rxy[m] = sum([(_x[kk]-_x_ave)*(_y[kk]-_y_ave) for kk in range(BIN_NUM)])\
        /math.sqrt(sum([(_x[kk]-_x_ave)**2 for kk in range(BIN_NUM)])\
            *(sum([(_y[kk]-_y_ave)**2 for kk in range(BIN_NUM)])))
        
        # end of f(a)
        # sys.exit()


        anh[m] = _beta_const - _int_prob
    # print anh
    np.savetxt("R_value.dat",rxy)
    print "Writed R vaule to R_value.dat"
    # np.savetxt('anh.dat',anh)     
    _delta = 0.
    for m in range(10):
        _delta += s_quantum[m] - anh[m]
    Delta_S = _delta * RC
    print "\t Delta_anh: %8.2f J/mol/K, %8.2f Kcal/mol." %(Delta_S, Delta_S/4186.8*temperature)
    print "After correction: %8.2f J/mol/K, %8.2f Kcal/mol." %(total_entropy- Delta_S, (total_entropy-Delta_S)/4186.8*temperature)
    # bp = numpy.zeros((truncate_num,nframes),dtype=numpy.float32)
    # for i in range(truncate_num):
    #         bp[i,k]

    # print bq[0]




def Main(top_file,traj_file,index_file,temperature=300,skip=1,begin=0,end=-1,CYCLES=1):
    '''
    Add some words here.
    '''
    #step 1: reading the trajectory.    
    BEGIN_TIME=Time.time()
    print "step 1: Reading trajectory"
    U       =MDAnalysis.Universe(top_file,traj_file)
    index   =Read_index_to_Inclass(index_file)
    Print_Index(index)
    while True:
        try:
            solute_index=int(raw_input("Choose the group for entropy calculation:"))
            break
        except:
            print "You should input a number."
            continue
    chose_group=index[solute_index].group_list
#    print chose_group

    if end== -1:
        nframes=U.trajectory.numframes-begin
    elif end > begin:
        nframes=end-begin
    elif begin > U.trajectory.numframes:
        print "The begin is out of the range of trajectory frames."
        sys.exit()
    else:
        print "The begin and the end of the frame seem not correct."
        sys.exit()

    natoms  =len(chose_group)
    print "\t Reading %d frames from trajectory file: %s" %(nframes,traj_file)

    
    #step 2: put data to traj_data matrix. get eigenvalues and eigenvectors.
    
    traj_data =np.zeros((nframes,natoms,3),dtype=np.float32)
    traj_data2=np.zeros((nframes,3*natoms),dtype=np.float32)
    traj_ave  =np.zeros((3*natoms),dtype=np.float32)
    covar_mat = np.zeros((3*natoms,3*natoms),dtype=np.float32)

    # sqrt_mass = [ math.sqrt(U.atoms[i].mass) for i in range(U.trajectory.numatoms)] # list for all atoms
    group_mass =np.repeat([math.sqrt(U.atoms[i].mass) for i in chose_group],3) 
    # list for choosing group. each atom repeat three times.
    # print group_mass
    # sys.exit()

    

    POINT=0
    ##  POINT used to replace the U.trajectory.frame in this part.
    for ts in U.trajectory:
#        temp_vect = np.zeros(3*natoms,dtype=np.float32)
        POINT += 1
        if POINT > end : 
            break
        elif POINT < begin:
            continue
#        elif POINT > end and end > 0:
#            break
#        else:
#            print "Note: reach here. begin=%6d,end=%6d,POINT=%6d" %(begin,end,POINT)
#            sys.exit()

        sys.stderr.write("\t Reading frame %8d\r" %POINT)
        sys.stderr.flush()

        for i,ai in list(enumerate(chose_group)):
            # traj_data[POINT-begin-1,i,0] = ts._x[ai-1] * sqrt_mass[ai-1]
            # traj_data[POINT-begin-1,i,1] = ts._y[ai-1] * sqrt_mass[ai-1]
            # traj_data[POINT-begin-1,i,2] = ts._z[ai-1] * sqrt_mass[ai-1]
            traj_data[POINT-begin-1,i,0] = ts._x[ai-1]
            traj_data[POINT-begin-1,i,1] = ts._y[ai-1]
            traj_data[POINT-begin-1,i,2] = ts._z[ai-1]

 ########################
 #Add fitting code here. from matrix traj_data to a new matrix.
    # print "\nstep 2: Fitting the trajectory."

    ref_coor=traj_data[0,:,:]

    ref_com = np.array([sum(ref_coor[:,0])/natoms,sum(ref_coor[:,1])/natoms,sum(ref_coor[:,2])/natoms])
    #ref_com means center of coordinate the reference.
    ref_coordinates = ref_coor - ref_com

    for k in range(natoms):
        for l in range(3):
            traj_data2[0,3*k+l]=ref_coordinates[k,l]    

    # traj_coordinates = traj_atoms.coordinates().copy()

    # nframes = len(frames)
    rmsd = np.zeros((nframes,))

    rot = np.zeros(9,dtype=np.float64)      # allocate space for calculation
    R = np.matrix(rot.reshape(3,3))

    for k in range(1,nframes):
        # shift coordinates for rotation fitting
        # selection is updated with the time frame

        sys.stderr.write("\t Fitting frame %8d\r" %(k+1))
        sys.stderr.flush()

        traj_coor=traj_data[k,:,:]
        x_com = np.array([sum(traj_coor[:,0])/natoms,sum(traj_coor[:,1])/natoms,sum(traj_coor[:,2])/natoms])
        traj_coordinates = traj_coor - x_com

        rmsd[k] = qcp.CalcRMSDRotationalMatrix(ref_coordinates.T.astype(np.float64),
                                               traj_coordinates.T.astype(np.float64),
                                               natoms, rot, None)
#        print rmsd
#a       sys.exit(0)
        R[:,:] = rot.reshape(3,3)

        # Transform each atom in the trajectory (use inplace ops to avoid copying arrays)
        # (Marginally (~3%) faster than "ts._pos[:] = (ts._pos - x_com) * R + ref_com".)
        # ts._pos   -= x_com
        # ts._pos[:] = ts._pos * R # R acts to the left & is broadcasted N times.
        # ts._pos   += ref_com
        new_coor=np.array(np.matrix(traj_coordinates)*np.matrix(R))
        for i in range(natoms):
            for j in range(3):
                traj_data2[k,3*i+j]=new_coor[i,j]

    del traj_data
    # np.savetxt("rmsd.dat",rmsd)
    # print traj_data2[:10,0]
    # sys.exit(0)
 ########################
    print "\nstep 3: Creating covariance matrix"
    print "\t Generazing the (%d X %d) covariance matrix" %(3*natoms,3*natoms)

#        traj_data[ts.frame-1]=temp_vect


##########################################
    fp = open("S-t.dat",'w')

    for cycle in range(CYCLES):
        QHE(traj_data2,group_mass,nframes/CYCLES*(cycle+1))

    fp.close()
        # sys.exit()

##########################################
    



def Usage():
    print "Entro_Analysis.py -f <traj_file> -p <top_file> -n <index_file> -t [temperature] -b [begin frame] -e [end frame] -k [skip] -c [cycles]"

def CheckArgs():
    if len(sys.argv) > 1:
        traj_file   =""
        top_file    =""
        index_file  =""
        temperature =300
        skip        =1
        begin_time  = 0
        end_time    = 0
        cycles      = 1
        opts,args=getopt.getopt(sys.argv[1:],"f:p:t:k:n:b:e:")
        for a,b in opts :
            if a =="-f":
                traj_file = b
            elif a == "-p":
                top_file = b
            elif a == "-t":
                temperature=int(b)
            elif a == "-k":
                skip = int(b)
            elif a == "-n":
                index_file=b
            elif a == "-b":
                begin_time = int(b)
            elif a == "-e":
                end_time = int(b)
            elif a == "-c":
                cycles = int(b)


        if os.path.isfile(traj_file) and os.path.isfile(top_file) and os.path.isfile(index_file):
            Main(top_file,traj_file,index_file,temperature,skip,begin_time,end_time,cycles)
    else:
        Usage()
        sys.exit()

if __name__=="__main__":
#    QH_entro("traj-closest.pdb","permute.xtc","index.ndx")
    CheckArgs()
