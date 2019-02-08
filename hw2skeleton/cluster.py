from .utils import Atom, Residue, ActiveSite

#pip freeze > requirements.txt #Miriam Trick
#functions for similarity metrics

'''
Bio.PDB Package: Hamelryck, T., Manderick, B. (2003) PDB parser and structure class implemented in Python. Bioinformatics 19: 2308–2310

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
    print(file_location+"*.pdb")
    PDB_files=glob.glob(file_location+"*.pdb") #grabbing all of the files within the PDB folder
    p = PDBParser(PERMISSIVE=1)
    sim_dict = {}
    for i in PDB_files: #we are iterating over each PDB file in our directory
        name=(re.findall('\d+.pdb$', i)[0])
        structure = p.get_structure(i, i) #grabbing biopython structure (generator object)
        model=structure[0] #this is extracting the model. Some structures have multiple models, but all of our structures only have one model
        residues=model.get_residues() #getting each residue. For my similarity metric, I am grabbing the coordinates of the alphaCarbon atom
        x_res=[]
        y_res=[]
        z_res=[]
        residues2,residues=itertools.tee(residues) #this is storing 2 generator objects. I know this is not great spacewise, but in most PDBs there are very few residues
        for a in residues: # we are extracting the x,y,z location of the alpha carbon for each residue.
            x_res.append(a["CA"].get_vector()[0]) 
            y_res.append(a["CA"].get_vector()[1])
            z_res.append(a["CA"].get_vector()[2])
        centroid_x, centroid_y,centroid_z = np.mean(x_res), np.mean(y_res), np.mean(z_res) #finding the centroid of the enzyme binding space
        dist=[]
        for b in residues2: #for each residue in each PDB, calculating the eucledian distance between the CA on each residue and the centriod
            dist.append(math.sqrt(((centroid_x-(b["CA"].get_vector()[0]))**2)+((centroid_y-(b["CA"].get_vector()[1]))**2)+((centroid_z-(b["CA"].get_vector()[2]))**2)))
        PDB_3D_dist=np.mean(dist) #this is the distance or space measurement for this PDBs enzyme active site.
        #print(PDB_3D_dist)
        sim_dict[PDB_3D_dist]=name #adding the distance metric to a dictionary corresponding to the name of the PDB 
    #print(sim_dict)
    with open('PDB_HW2_sim_metric.pickle', 'wb') as handle: #exporting this dictionary as a pickle file 
        pickle.dump(sim_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

##################PARTITIONING ALGO################################

import os
import pandas as pd
import pickle
import math
import numpy as np

'''
I have set up this function to take in a dictionary file which has keys==the volume similarity metric
(or some numerical representation of the similarity between PDBs) and values==PDB name.

Initially find 3 starting points to start as the cluster means. We are choosing the first three items in the list.
We then iterate through all items in the list, determine the square root distance away from each cluster mean.
For each item, we assign it to the cluster that it is closest to the mean.
Once everything has been assigned, we find a new mean for that cluster.


'''


def select_initial_values(dist_list): #k is the number of clusters #figure out how to re-iterate over this for k_i
    '''This function is finding & stating the initial mean values we will have for kmeans.
    We are using the first three items in the dictionary.
    Input: dictionary with PDB=sim_metric
    Output: 3 Initial values
    '''
    k1=float(list(dist_list)[0]) #grabbing the similarity metric from the keys in the dictionary
    k2=float(list(dist_list)[1])
    k3=float(list(dist_list)[2])
    print('k1,2,3:')
    print(k1, k2, k3)
    return k1,k2,k3


def assign_new_clusters(dist_list,k1,k2,k3):
    '''This funciton will be used once we iterate over all the items in the entire list, assign them to a cluster and then get the mean of the newly created cluster.
    Input: dictionary with PDB=sim_metric, 3 mean values for each cluster
    Output: 3 means for new cluster, list of 3 clusters

    '''
    c_1=[] #creating each cluster
    c_2=[]
    c_3=[]
    for i in range(len(dist_list)-1): #for length of the similarity metric,
        i=i+1 # i is going to grab the next item in the list
        dist_k1=math.sqrt((float(list(dist_list)[i])-k1)**2) ##sqrt distance of each item in the list from the median of each cluster
        dist_k2=math.sqrt((float(list(dist_list)[i])-k2)**2)
        dist_k3=math.sqrt((float(list(dist_list)[i])-k3)**2)
        min_dist=min(dist_k1,dist_k2,dist_k3) #for i, determine what cluster mean it is closest to
        #assign each value to the cluster it belongs to (based on the distance)
        if min_dist==dist_k1:
            c_1.append(float(list(dist_list)[i]))#assign to c_1
        elif min_dist==dist_k2:
            c_2.append(float(list(dist_list)[i]))#assign to c_1
        else:
            c_3.append(float(list(dist_list)[i]))#assign to c_3
    k1=np.median(c_1) #find the median for each cluster
    k2=np.median(c_2)
    k3=np.median(c_3)
    return k1,k2,k3,c_1,c_2,c_3 #return the clusters and the new medians & DO IT AGAIN!

def k_median(pickle_input):
    k1,k2,k3=select_initial_values(pickle_input) #figure out input/output
    for i in range(10): #iterate until ?
        k1,k2,k3,c_1,c_2,c_3=assign_new_clusters(pickle_input,k1,k2,k3)
    return c_1, c_2, c_3

'''
All of the clustering above was done on the similarity metrics
'''
#connect values back to pdb ids
def replace_values(cluster, pickle_input):
    new_list = []
    for i in cluster:
        if isinstance(i, list):
              new_list.append(replace_vals(i, pickle_input))
        elif i in pickle_input:
            new_list.append(pickle_input.get(i))
        else:
            new_list.append(i)
    return new_list
