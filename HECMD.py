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
import numpy as np
import MDAnalysis
import MDAnalysis.core.qcprot as qcp
import time as Time

K= 1.3806504e-23
'''It's the Boltzmann constant'''
HBA = 1.05457148e-34
H   = 6.62606957e-34
'''It's the Planck constant'''

RC= 8.3145
'''It's the idea gas constant'''
MP=1.66055402e-27
'''1 a.m.u. = Proton mass [kg] '''



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

def histogram_2D(list_1,list_2,BINS):
    length= len(list_1)
    # print list_1
    # print length
    # print np.shape(list_1)
    # sys.exit()
    min_1 = np.amin(list_1)
    min_2 = np.amin(list_2)
    max_1 = np.amax(list_1)
    max_2 = np.amax(list_2)
    _interval_1 = (max_1 - min_1) / BINS
    _interval_2 = (max_2 - min_2) / BINS
    # print _interval_1,_interval_2
    _hist = np.zeros((BINS,BINS),dtype = np.int32   )
    _prob = np.zeros((BINS,BINS),dtype = np.float32 )
    for i in range(length):
        # print int((list_1[i]-min_1)/_interval_1),int((list_2[i]-min_2)/_interval_2)
        try:
            _hist[int((list_1[i]-min_1)/_interval_1),int((list_2[i]-min_2)/_interval_2)] += 1
        except:
            pass
    for i in range(BINS):
        for j in range(BINS):
            _prob[i,j] = (_hist[i,j]) / (length*_interval_1*_interval_2)
            # if _prob[i,j] == 0:
            #     _prob[i,j] = 1e-10

    int_prob = 0.0
    # print _prob
    # sys.exit()
    for i in range(BINS):
        for j in range(BINS):
            if _prob[i,j] > 0.0:
                int_prob += math.log(_prob[i,j]*1e20/(MP))*_prob[i,j] * _interval_1 * _interval_2

    return int_prob


def R_VALUE(prob,eig_value,bin_list,BINS):

    _fa = [ math.sqrt(abs(-1*eig_value*math.log(_prob_item*math.sqrt(2*math.pi*eig_value)))) for _prob_item in prob ]
    _x = np.zeros([len(prob)],dtype=np.float32)
    _y = np.zeros([len(prob)],dtype=np.float32)
    for kk in range(len(prob)):
        _x[kk] = 0.5*(bin_list[kk]+bin_list[kk+1])
        if _x[kk] < 0 :
            _y[kk] = -1* _fa[kk]
        else:
            _y[kk] = _fa[kk]
    _x_ave = sum(_x) / BINS
    _y_ave = sum(_y) / BINS
    _rxy = sum([(_x[kk]-_x_ave)*(_y[kk]-_y_ave) for kk in range(BINS)])\
    /math.sqrt(sum([(_x[kk]-_x_ave)**2 for kk in range(BINS)])\
        *(sum([(_y[kk]-_y_ave)**2 for kk in range(BINS)])))   

    return _rxy


def QHE(traj_data2,group_mass,temperature,cycle_frame,eigen_thr = 1e-7,corr_items=10):
    global HBA, H, K, MP, RC

    _mat_size = len(group_mass)
    traj_ave  = np.zeros((_mat_size),dtype=np.float32)
    covar_mat = np.zeros((_mat_size,_mat_size),dtype=np.float32)

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
     

    print "step 4: Diagnoalizing the covariance matrix"
    eigenvalues,eigenvectors = np.linalg.eig(covar_mat)
    #    print eigenvalues
    # np.savetxt("eigenvalues.dat",eigenvalues)

    truncate_num = 0
    for i in range(_mat_size):
        if eigenvalues[i] < eigen_thr:
            truncate_num = i
            break
        elif eigenvalues[-1] > eigen_thr:
            truncate_num = _mat_size
    print "\t Using %d eigenvalues for calculating entropy." %truncate_num
        
    
    print "step 5: Calculating the quasiharmonic entropy" 

    alpha_const  =HBA*(1e10)/math.sqrt(K*temperature*MP)
    alpha     =[alpha_const/math.sqrt(eigenvalues[i]) for i in range(truncate_num)]
    s_quantum =[alpha[i]/(math.exp(alpha[i])-1) - math.log(1-math.exp(-alpha[i])) for i in range(truncate_num)]
    # np.savetxt("quantum.dat",s_quantum)
    
    total_entropy=sum(s_quantum)*RC
    print "\t Entopy: %8.2f J/mol/K, %8.2f Kcal/mol." %(total_entropy,total_entropy/4186.8*temperature)
    # fp.write("\t Entopy: %8.2f J/mol/K, %8.2f Kcal/mol.\n" %(total_entropy,total_entropy/4186.8*temperature))
    # fp.flush()


######################
    print "step 6: Calculating the anhormonic entropy"
   
    _bp_temp = [eigenvectors[:,i]*group_mass for i in range(corr_items)]   
    _delta_temp= np.zeros((_mat_size,cycle_frame),dtype=np.float32)
    for i in range(cycle_frame):
        for j in range(_mat_size):
            _delta_temp[j,i]=traj_data2[i][j]-traj_ave[j]
    bp=np.matrix(_bp_temp)*np.matrix(_delta_temp)
    bp = np.array(bp)
    # for i in range(truncate_num):
        # np.savetxt("dist_%d.dat" %(i+1),bp[i])
    anh=np.zeros([corr_items],dtype=np.float32)

    _beta_const = 0.5*(1-math.log(H*HBA/K/temperature))
    BIN_NUM = 100
    # rxy=np.zeros([truncate_num],dtype=np.float32)
    
    # print _beta_const
    for m in range(corr_items):
        _hist,_bins=np.histogram(bp[m],bins=np.linspace(np.amin(bp[m]),np.amax(bp[m]),BIN_NUM+1))   
        _sum_hist= sum([(_hist[i]+1.0/BIN_NUM)*(_bins[i+1]-_bins[i]) for i in range(BIN_NUM)])
        _prob = [(_item+1.0/BIN_NUM)/_sum_hist for _item in _hist]
        # _int_prob = 0.0
        # for i,_prob_i in enumerate(_prob):
        #     if _prob_i > 0.0:
        #         _int_prob += math.log(_prob_i*1e10/math.sqrt(MP))*_prob_i *(_bins[i+1]-_bins[i])
        _int_prob =sum([math.log(_prob[i]*1e10/math.sqrt(MP))*_prob[i]*(_bins[i+1]-_bins[i]) for i in range(len(_hist))]) #care the unit

        # rxy[m] = R_VALUE(_prob,eigenvalues[m],_bins,BIN_NUM)      

        anh[m] = _beta_const - _int_prob
    # print anh
    # np.savetxt("R_value.dat",rxy)
    # print "Writed R vaule to R_value.dat"
    # np.savetxt('anh.dat',anh)     
    _delta_ah = 0.
    for m in range(corr_items):
        _delta_ah += anh[m] - s_quantum[m]
    Delta_ah = _delta_ah * RC
    print "\t Delta_anh: %8.2f J/mol/K, %8.2f Kcal/mol." %(Delta_ah, Delta_ah/4186.8*temperature)
    print "\t After correction: %8.2f J/mol/K, %8.2f Kcal/mol." %(total_entropy + Delta_ah, (total_entropy + Delta_ah)/4186.8*temperature)


    print "step 7: Calculating the pairwise correlation"

    delta_pch = np.zeros((corr_items,corr_items),dtype=np.float32) # pairwise correlation for each unit
    for m in range(corr_items-1):
        for n in range(m+1,corr_items):
            delta_pch[m,n] = ( 2* _beta_const - histogram_2D(bp[m],bp[n],20) ) - anh[m] - anh[n]
            # print ( 2* _beta_const - histogram_2D(bp[m],bp[n],20) ),anh[m],anh[n], delta_pch[m,n]

    Delta_pch = np.sum(delta_pch) * RC
    print "\t Delta_pch: %8.2f J/mol/K, %8.2f Kcal/mol." %(Delta_pch, Delta_pch/4186.8*temperature)
    print "\t After correction: %8.2f J/mol/K, %8.2f Kcal/mol." %(total_entropy + Delta_ah+Delta_pch, \
        (total_entropy + Delta_ah + Delta_pch)/4186.8*temperature)




def Main(top_file,traj_file,index_file,temperature=300,begin=0,end=-1,CYCLES=1):
    '''
    Add some words here.
    '''
    #step 1: reading the trajectory.    
    # BEGIN_TIME=Time.time()
    print "step 1: Reading trajectory"
    U       =MDAnalysis.Universe(top_file,traj_file)
    index   =Read_index_to_Inclass(index_file)
    Print_Index(index)

    while True:
        try:
            solute_index=int(raw_input("Choosing the group for entropy calculation:"))
            break
        except:
            print "You should input a number."
            continue
    chose_group=index[solute_index].group_list
#    print chose_group
    print "\t Loading the trajectory file %s, please wait..." %(traj_file)
    if end== -1:
        nframes=(U.trajectory.numframes-begin)
    elif end > begin:
        nframes=(end-begin)
    elif begin > U.trajectory.numframes:
        print "The begin is out of the range of trajectory frames."
        sys.exit()
    else:
        pass
        # print "The begin and the end of the frame seem not correct."
        # sys.exit()

    natoms  =len(chose_group)
    print "\t Reading %d frames from trajectory file: %s" %(nframes,traj_file)
    
    traj_data =np.zeros((nframes,natoms,3),dtype=np.float32)
    traj_data2=np.zeros((nframes,3*natoms),dtype=np.float32)
    

    POINT=0
    ##  POINT used to replace the U.trajectory.frame in this part.
    for ts in U.trajectory:
        

        if (ts.frame > end-1) and end > 0: 
            break
        elif ts.frame < begin:
            continue
        # elif (ts.frame-begin) % skip != 0:
        #     continue

        sys.stderr.write("\t Reading frame %8d\r" %ts.frame)
        sys.stderr.flush()

        for i,ai in list(enumerate(chose_group)):
            traj_data[POINT,i,0] = ts._x[ai-1]
            traj_data[POINT,i,1] = ts._y[ai-1]
            traj_data[POINT,i,2] = ts._z[ai-1]
        POINT += 1

 ########################
 #Add fitting code here. from matrix traj_data to a new matrix.
    print "\nstep 2: Fitting the trajectory."

    ref_coor=traj_data[0,:,:]
    ref_com = np.array([sum(ref_coor[:,0])/natoms,sum(ref_coor[:,1])/natoms,sum(ref_coor[:,2])/natoms])
    ref_coordinates = ref_coor - ref_com
    # print ref_coordinates

    for k in range(natoms):
        for l in range(3):
            traj_data2[0,3*k+l]=ref_coordinates[k,l]    

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
        R[:,:] = rot.reshape(3,3)
        new_coor=np.array(np.matrix(traj_coordinates)*np.matrix(R))
        for i in range(natoms):
            for j in range(3):
                traj_data2[k,3*i+j]=new_coor[i,j]

    del traj_data
    # np.savetxt("rmsd.dat",rmsd)

 ########################
    print "\nstep 3: Creating covariance matrix"
    print "\t Generazing the (%d X %d) covariance matrix" %(3*natoms,3*natoms)

    group_mass =np.repeat([math.sqrt(U.atoms[i-1].mass) for i in chose_group],3) 
    for cycle in range(CYCLES):
        QHE(traj_data2,group_mass,temperature,nframes/CYCLES*(cycle+1))


def File_input():
    print "%10s   %8s   %20s    %s " %("Option","Type","Filename","Description")
def Print_line():
    print "-"*65
def Show(Option,Type,value,description):
    print "%10s   %8s   %20s    %s" %(Option,Type,value,description)

# def Usage():
#     print "Entro_Analysis.py -f <traj_file> -p <top_file> -n <index_file> -t [temperature] -b [begin frame] -e [end frame] -c [cycles]"

def Usage(coor_file="coor_file",traj_file="traj_file",ndx_file="index.ndx",log_file="log_file",\
    eigen_file = "",\
    temperature=300.0,begin=0,end=-1,cycles=1):
    '''
    print the usage information.
    '''
    print ""
    File_input()
    Print_line()
    Show("-p","Input",coor_file,"Structure file: gro pdb etc.")
    Show("-f","Input",traj_file,"Trajectory: xtc trr.")
    Show("-n","Input",ndx_file, "Index file.")
    Show("-g","Output",log_file,"Log file.")
    Show("-e","Output",eigen_file,"eigenvalues file.")

    print 
    Show("--begin","int",str(begin),"begin time.")
    Show("--end","int",str(end),"end time.")
    Show("--temp","float",str(temperature),"Kelvin temperature")
    Show("--cycle","int",str(cycles),"Cycles for entropy calculation")

    print ""

    if os.path.isfile(traj_file) and os.path.isfile(coor_file) and os.path.isfile(ndx_file):
        Main(coor_file,traj_file,ndx_file,temperature,begin,end,cycles)



def CheckArgs():
    if len(sys.argv) > 1:
        traj_file   = ""
        top_file    = ""
        index_file  = ""
        log_file    = ""
        eigen_file  = ""

        temperature =300
        # skip        =1
        begin_time  = 0
        end_time    = -1
        cycles      = 1
        opts,args=getopt.getopt(sys.argv[1:],"f:p:g:t:n:b:e:",["begin=","end=","temp=","cycle="])

        for a,b in opts :
            if a =="-f":
                traj_file = b
            elif a == "-p":
                top_file = b
            elif a == "--temp":
                temperature=float(b)
            # elif a == "-k":
            #     skip = int(b)
            elif a == "-n":
                index_file=b
            elif a == "--begin":
                begin_time = int(b)
            elif a == "--end":
                end_time = int(b)
            elif a == "--cycle":
                cycles = int(b)

        Usage(top_file,traj_file,index_file,log_file,eigen_file,temperature,begin_time,end_time,cycles)

    else:
        Usage()
        sys.exit()

if __name__=="__main__":
    CheckArgs()
