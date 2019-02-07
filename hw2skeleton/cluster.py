from .utils import Atom, Residue, ActiveSite

#functions for similarity metrics

'''
Bio.PDB Package Citation

Bio.PDB does XXXX

'''
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



def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!

    return []


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
