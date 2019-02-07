from .utils import Atom, Residue, ActiveSite

$ pip freeze > requirements.txt #Miriam Trick
#functions for similarity metrics

'''
Bio.PDB Package Citation

Bio.PDB does XXXX

'''
#packages
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
import sys
import numpy as np
import math
import os
import glob
import itertools
import re
import pickle
import pandas as pd


def similarity_metric(file_location): ### PUT THIS ALL TOGETHER
    PDB_files=glob.glob(file_location+"*.pdb") #grabbing all of the files within the PDB folder
    p = PDBParser(PERMISSIVE=1)
sim_dict = {}
for i in PDB_files:
    name=(re.findall('\d+.pdb$', i)[0])
    print(name)
    structure = p.get_structure(i, i) #grabbing biopython structure (generator object)
    model=structure[0] #transitioning to model ????
    residues=model.get_residues() #getting each residue. For my similarity metric, I am grabbing the coordinates of the alphaCarbon atom
    x_res=[]
    y_res=[]
    z_res=[]
    residues2,residues=itertools.tee(residues) #this is storing 2 generator objects. I know this is not great spacewise, but in most PDBs there are very few residues
    for a in residues:
        x_res.append(a["CA"].get_vector()[0])
        y_res.append(a["CA"].get_vector()[1])
        z_res.append(a["CA"].get_vector()[2])
    centroid_x, centroid_y,centroid_z = np.mean(x_res), np.mean(y_res), np.mean(z_res)
    dist=[]
    for b in residues2:
        #print(i)
        dist.append(math.sqrt(((centroid_x-(b["CA"].get_vector()[0]))**2)+((centroid_y-(b["CA"].get_vector()[1]))**2)+((centroid_z-(b["CA"].get_vector()[2]))**2)))
    PDB_3D_dist=np.mean(dist) #put this into a dictionary/hash table
    print(PDB_3D_dist)
    sim_dict[PDB_3D_dist]=name
print(sim_dict)
with open('PDB_HW2_sim_metric.pickle', 'wb') as handle:
    pickle.dump(sim_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

##################PARTITIONING ALGO################################

#seperate functions
def select_initial_values(dist_list, k): #k is the number of clusters #figure out how to re-iterate over this for k_i
    k1=float(list(dist_list)[0]) #assuming this value is stored in pickle keys
    k2=float(list(dist_list)[1])
    k3=float(list(dist_list)[2])
    print('k1,2,3:')
    print(k1, k2, k3)
    return k1,k2,k3

def assign_new_clusters(dist_list,k1,k2,k3):
    print('k1')
    print(k1)
    c_1=[k1]
    c_2=[k2]
    c_3=[k3]
    for i in range(len(dist_list)-1): #deal with length of pickle
        i=i+1
        dist_k1=math.sqrt((float(list(dist_list)[i])-k1)**2) ##sqrt distance
        dist_k2=math.sqrt((float(list(dist_list)[i])-k2)**2)
        dist_k3=math.sqrt((float(list(dist_list)[i])-k3)**2)
        min_dist=min(dist_k1,dist_k2,dist_k3)
        #print(min_dist)
        if min_dist==dist_k1:
            c_1.append(float(list(dist_list)[i]))#assign to c_1
        elif min_dist==dist_k2:
            c_2.append(float(list(dist_list)[i]))#assign to c_1
        else:
            c_3.append(float(list(dist_list)[i]))
    k1=np.median(c_1)
    k2=np.median(c_2)
    k3=np.median(c_3)
    print('K1')
    print(k1)
    return k1,k2,k3,c_1,c_2,c_3

##MAIN FUNCTION##
def k_means(pickle_input, k):
    k1,k2,k3=select_initial_values(pickle_input,k) #figure out input/output
    for i in range(100): #iterate until ?
    #while k1_new not equal to k1...  #confidence metric #need to output name with distance metric
        k1,k2,k3,c_1,c_2,c_3=assign_new_clusters(pickle_input,k1,k2,k3)
    print('c_1')
    print(c_1)
    return c_1, c_2, c_3

    

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!

    return 


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
