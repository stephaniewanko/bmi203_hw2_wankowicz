###comparison functions###
def silhouette(cluster1,cluster2):
    '''
    The silhouette score compares how similar each item is to all the other items in the cluster it belongs to compared to another cluster. 
    INPUT: 3 clusters as lists
    OUTPUT: Silhouette Score
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
    s=(mean(a_list)-mean(b_list))/(min(a_list)-min(b_list)
    #return s 

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
