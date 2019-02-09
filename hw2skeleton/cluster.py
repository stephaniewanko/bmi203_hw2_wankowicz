#from .utils import Atom, Residue, ActiveSite I did not use any of these functions. 

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



#def testing_test_function(num): #SO META
#    print(num)


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

#This function will perform the k median function to return 3 clusters. 
    k1,k2,k3=select_initial_values(pickle_input) #figure out input/output
    for i in range(1000): 
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


#######hierarchical_clustering#####

'''
We are going to to do agglomerative clustering where each item (PDB) starts in its own cluster.
We will then iteratively find the most similar clusters, combine them until everything is combined.
We are going to use complete linkage, which measures closesness of two clusters by the most dissimilar members.
'''

'''
1)iterate over pickle an assign each to its own cluster
2)for i in clusters: calculate the distance from each other cluster via . Find the smallest difference. # do you need to calculate these all? Create distance matrix
bind those two together, and assign it to cluster 1a
3)input the next set of clusters and do it all over again.

At the end of the clustering algorithm, there is a function that will return the distance values and the pdb names once you give the number of clusters (aka the cut off you want to look at)
'''

def mediod_link(i,j):
    '''
    The medido link looks at the median eucledian distance between two clusters
    i and j should be lists of values that correspond with a cluster
    Input: 2 list of similarity metrics
    Output: The eucledian distance between the medians of both lists
    '''
    dist=math.sqrt((np.median(j)-np.median(i))**2)
    return(dist)

def dict_to_list(input_pickle):
    '''
    Since we take in a dictionary, we want to convert that dictionary to a list of the similarity metrics to do the clustering.
    At the end, we will use the dictionary to hook things back up.
    '''
    dist_metrics=list(input_pickle.keys()) #putting the keys into a list
    return dist_metrics #returning that list

def most_sim_clusters(master_list):
    '''
    This is going to look at the mediod distance between each existing cluster and determine which two existing clusters are the closest.
    We are deciding which ones are most similar based on the mediod linkage.
    INPUT: List of clusters
    OUTPUT: most closely related clusters
    '''
    min=None #We are initiating this by stating that the minimum is empty.
    for i in master_list: #for every cluster in the master list...
        for j in master_list: #compare to every other cluster in the master list
            if any(j) not in i: #determine if j is in i, if so, skip
                min_tmp=(mediod_link(i,j)) #determine the mediod link and assign it to the temporary minimum
                if i==j: #since we are looping over the list, we want to make sure i and j are not the same value
                    continue
                else:
                    if min is None: #if this is the first mediant distance you are looking at, assign the min value to the mediod linkage calculation, and i&j
                        min=min_tmp #assign min value to the temp if the min is empty
                        i_min=i
                        j_min=j
                    elif min_tmp<min: #for each mediod linkage value, if it is less than the stored min value, reassign min value, i & j
                        min=min_tmp
                        i_min=i
                        j_min=j
                    else: #if the mediod value not smaller than min, keep looping through
                        continue
    return i_min, j_min #return the lists that have the minimum mediod linkage


def merging_clusters(i,j,input_list):
    '''
    This will merege two clusters
    INPUT: The two clusters you want to merge (i&j) and the master cluster list
    OUTPUT: An updated master cluster list with the newly merged cluster
    '''
    new_cluster=i+j #creating new cluster from the two most similar clusters
    input_list.remove(i) #remove item from the master cluster list (you are going to put it in the new cluster)
    input_list.remove(j) #remove item from the master cluster list (you are going to put it in the new cluster)
    input_list.append(new_cluster) #put the newly merged cluster in with the master cluster list
    return input_list #return the new master cluster list

def hier_cluster(input_dict,num_clusters):
    '''
    This funciton will call the input dictionary and the number of clusters you want to cluster to.
    If you want to cluster all the way up the dendrogram, then select the numb of clusters as 1.
    It will output a dictionary. To get the values and corresponding PDB names from a specific number of clusters, call the return cluster function.

    INPUT: The similarity dictionary, the number of clusters you want to cluster to
    OUTPUT: Dictionary with the number of clusters:list of clusters
    '''

    master_list=dict_to_list(input_dict) #this is function that will convert the input dictionary to a list
    output_dict={} #we are going to store the number of clusters and cluster values into a dictionary
    clusters = [[master_list[i]] for i in range(len(master_list))] #this is going to split every item in the master list into their own cluster since we are doing agglomerative clsutering.
    cluster_num=len(clusters) #we need to determine the number of clusters that we start with
    while cluster_num>num_clusters: #we are iterating over this until we have combined all of the values into one cluster
        output_dict[cluster_num]=clusters #we are going to store the number of clusters and cluster values into a dictionary
        i,j=most_sim_clusters(clusters) #determining which clusters are closest
        clusters=merging_clusters(i,j,clusters) #merging the two closest clusters together
        cluster_num=len(clusters) # we want to keep track of the number of cluster we have. This will determine when the while loop stops.
    return output_dict



def replace_values(cluster, pickle_input):
    '''
    This is going to take a list and replace the similarity metric with the PDB name from the original dictionary
    INPUT: The cluster you want to get the PDB names from, the dictionary with the PDB names and similarity metric
    OUTPUT: The cluster with the PDB names
    '''
    new_list = []
    for i in cluster: #for each item in the cluster
        if isinstance(i, list): #if there are nested list, unnest the list
            new_list.append(replace_values(i, pickle_input)) #replace the similarity metric with the PDB name
        elif i in pickle_input:
            new_list.append(pickle_input.get(i)) #add values onto the new list
        else:
            new_list.append(i)
        return new_list


def return_cluster(distance_dict, num_clusters, pickle_input):
    '''
    This funciton will give you back the data (both PDB names and similairty values) seperated into the number of clusters you want to look at.
    INPUT: The dictionary with # of clusters:clusters, the number of clusters you want to look at, and the dictionary with the PDB names and similarity metric.
    OUTPUT: List of the cluster with the similarity metrics and PDB name.
    '''
    for key, value in distance_dict.items():
        if key==num_clusters:
            cluster=value
        pdb_names=replace_values(cluster, pickle_input)
    return pdb_names, cluster

###comparison functions###
def jacard_index(cluster1, cluster2):
    '''
    To compare the two clustering algorithms, we are using the Jacard Index, which determines the length of the union of the two clusters over the intersection.
    We are going to the three cluster level.
    The higher the Jacard index, the better the clusters compare (more intersection per item in the two indexes).

    Input: one representative cluster from each algorithm
    Output: Jacard Index
    '''
    union = list(set(cluster1).union(cluster2))
    intersection = list(set(cluster1) & set(cluster2))
    return len(intersection)/len(union)


def silhouette(cluster1,cluster2):
    '''
    '''
    mean_cluster1=mean(cluster1) #we want to find the closest cluster to compare this to
    mean_cluster2=mean(cluster2)
    mean_cluster3=mean(cluster3)
    #find the closest cluster
    clu1_2=math.mean(sqrt((mean_cluster1-mean_cluster2)**2))
    clu1_3=math.mean(sqrt((mean_cluster1-mean_cluster3)**2))
    if min(clu1_2,clu1_3)==clu1_3:
        cluster2=cluster3
    for i in cluster1:
        a_list=[]
        b_list=[]
        for j in cluster1:
            a_list.append(math.mean(sqrt((j-i)**2))) #finding the eucledian distance between each values in the same cluster
        for j in cluster2:
            b_list.append(math.mean(sqrt((j-i)**2))) #finding the eucledian distance between each values in the other cluster

