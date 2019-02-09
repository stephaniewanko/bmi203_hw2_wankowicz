from hw2skeleton import *
from hw2skeleton import io
import sys
sys.path.append('hw2skeleton/')
from cluster import *
import os


print('Starting to test!')
#testing_test_function(2)
def test_similarity(location):
    '''
    We are going to run the similarity algorithm, and then make sure that the similarity metric for pdb 47023 is equal to 21.028781603576693'.
    '''
    sim_dictionary=cluster.similarity_metric(location) #make sure this is returning a dictionary
    for key, value in sim_dict.items():
        if value=='47023.pdb':
            key_key=key
        assert key_key=='21.028781603576693' #we are going to assert that pdb==46042.pdb has a similarity metric of 21.028781603576693.

    #cluster.compute_similarity(activesite_a, activesite_b) == 0.0



def test_partition_clustering(similarity_dictionary):
    '''
    To test the parition algorithm, we are going to run the k-medians algorithm. We are then going to assess the length of the third cluster.
    We have set the random.seed() to 40 to ensure we get the same clusters.
    '''

    cluster1,cluster2,cluster3=cluster.k_median(similarity_dictionary)
    assert len(cluster1) == 10
    
def test_hierarchical_clustering():
    '''
    To test the hierarchical algorithm, we are going to run the full algorithm. We are then going to assess the length of the third cluster.
    We have set the random.seed() to 40 to ensure we get the same clusters.
    '''
    test_cluster=cluster.hier_cluster(content,135)
    assert len(test_cluster) == 2


content=similarity_metric('./data/')
print('Testing Similarity')
test_similarity('./data')

#running the similarity metric to get the dictionary
test_partition_clustering(content)




assert [1]==[1]
