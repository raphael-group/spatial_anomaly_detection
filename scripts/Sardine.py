import matplotlib.pyplot as plt
import numpy as np
import itertools
import seaborn as sns
import math
from sklearn.cluster import KMeans
import meld 
from sklearn.decomposition import PCA
from scipy.spatial import KDTree
import pandas as pd

def create_spatialKDTrees(coordinate_set_1, coordinate_set_2):
    # construct the KDTrees for later use
    KDTrees_for_replicates = []   
    for coordinates_1 in coordinate_set_1:
        coordinate_KDTree = KDTree(coordinates_1)                          
        KDTrees_for_replicates.append(coordinate_KDTree)
    for coordinates_2 in coordinate_set_2:
        coordinate_KDTree = KDTree(coordinates_2)                          
        KDTrees_for_replicates.append(coordinate_KDTree)
        
    return KDTrees_for_replicates

def calculate_anomaly_score_localManifold(coordinates, 
                              spatial_KDTrees, 
                              gene_feature_replicate_matrices, 
                              num_replicates=1,
                              neighborhood_size=20,
                              knn=8,
                              beta=15): 

    """
    Run Sardine - local Manifold version

    Parameters:
    coordinates: a numpy array of size n x 2
        The numpy array that contains the points to comute the anomaly score at.
    spatial_KDTrees: list of KDTree objects
        This list contains 1 KDTree object per replicate, structured as:
        [condition_0_rep_1_KDTree, condition_0_rep_2_KDTree, ..., condition_1_rep_k_KDTree]
    gene_feature_replicate_matrices: list of gene expression matrices
        A list of gene expression matrices, 1 for each coordinate, matching the structure of the spatial_KDTrees list
    num_replicates: int
        The number of replicates used, assumed to be same for both conditions
    neighborhood_size: int
        The r parameter that controls how many spots from each condition Sardine should include the spatial windows R_i
    beta: int
        hyperparamter for the MELD plug in density estimator, which controls the Laplacian smoothing term

    Returns:
    anomaly_scores: list of floats
        The output of one float anomaly score between 0 and 1, corresponding to the coordinate matrix input.
    """
    
    count = 0
    anomaly_scores_local_manifold = []
    joint_feature_mat = np.vstack(gene_feature_replicate_matrices)
    for row in range(coordinates.shape[0]):
        # 
        x,y = coordinates[row] # current spot location
        neighborhood_expression = []
        # in practice often the replicates are collapsed into a single tuple and the 
        # num_replicates is set to 1, but this handles the case where we want to 
        # sample equally from each replicate
        for i in range(0,2*num_replicates):
            spatial_kdtree = spatial_KDTrees[i]
            distances, indices = spatial_kdtree.query((x,y), k=neighborhood_size)

            replicate_gene_matrix = gene_feature_replicate_matrices[i]
            neighborhood_gene_expression = replicate_gene_matrix[indices] # get the gene feature matrix

            neighborhood_expression.append(neighborhood_gene_expression)

        # run MELD 
        neighborhood_mat = np.asarray(neighborhood_expression).reshape(2*num_replicates*neighborhood_size, -1)
        condition_vector = np.concatenate((np.full(neighborhood_size * num_replicates, 1),
                                           np.full(neighborhood_size * num_replicates, 2)))

        meld_op = meld.MELD(beta=beta, knn=knn)
        sample_densities = meld_op.fit_transform(neighborhood_mat, condition_vector)
        sample_likelihoods = meld.utils.normalize_densities(sample_densities)

        # get the likelihoods for the perturbed condition
        second_condition_likelihoods_for_second_samples = sample_likelihoods[2].to_numpy()[neighborhood_size*num_replicates:]

        # store the resulting 
        local_manifold_mean = np.mean(second_condition_likelihoods_for_second_samples)
        anomaly_scores_local_manifold.append(local_manifold_mean)

    return anomaly_scores_local_manifold



def calculate_anomaly_score(coordinates, 
                               spatial_KDTrees, 
                               gene_feature_replicate_matrices, 
                               num_replicates=1, 
                               neighborhood_size=20,
                              knn=8,
                              beta=15):    

    """
    Run Sardine

    Parameters:
    coordinates: a numpy array of size n x 2
        The numpy array that contains the points to comute the anomaly score at.
    spatial_KDTrees: list of KDTree objects
        This list contains 1 KDTree object per replicate, structured as:
        [condition_0_rep_1_KDTree, condition_0_rep_2_KDTree, ..., condition_1_rep_k_KDTree]
    gene_feature_replicate_matrices: list of gene expression matrices
        A list of gene expression matrices, 1 for each coordinate, matching the structure of the spatial_KDTrees list
    num_replicates: int
        The number of replicates used, assumed to be same for both conditions
    neighborhood_size: int
        The r parameter that controls how many spots from each condition Sardine should include the spatial windows R_i
    beta: int
        hyperparamter for the MELD plug in density estimator, which controls the Laplacian smoothing term

    Returns:
    anomaly_scores: list of floats
        The output of one float anomaly score between 0 and 1, corresponding to the coordinate matrix input.
    """

    anomaly_scores = []
    joint_feature_mat = np.vstack(gene_feature_replicate_matrices)
    
    meld_op = meld.MELD(beta=beta, knn=knn)
    meld_op.fit(joint_feature_mat)
    
    # 
    for row in range(coordinates.shape[0]):
        # 
        x,y = coordinates[row] # current spot
        sample1_indices = []
        sample2_indices = []
        for i in range(0,2*num_replicates):
            spatial_kdtree = spatial_KDTrees[i]
            distances, indices = spatial_kdtree.query((x,y), k=neighborhood_size)

            # save the indices (by replicates)
            start = 0

            # in practice often the replicates are collapsed into a single tuple and the 
            # num_replicates is set to 1, but this handles the case where we want to 
            # sample equally from each replicate
            for j in range(0,i):
                start+=gene_feature_replicate_matrices[j].shape[0]
            if i < num_replicates:
                sample1_indices.append(start+indices)
            else:
                sample2_indices.append(start+indices)

        sample1_indices = np.asarray(sample1_indices).flatten()
        sample2_indices = np.asarray(sample2_indices).flatten()

        # create the condition vector on only these samples
        condition_vector = np.zeros(joint_feature_mat.shape[0])
        condition_vector[sample1_indices] = 1
        condition_vector[sample2_indices] = 2

        sample_densities = meld_op.transform(condition_vector)
        sample_likelihoods = meld.utils.normalize_densities(sample_densities[[1,2]])
        print(f"Computing spot {row}'s' average anomaly score.")
        # get the likelihoods for the perturbed condition
        condition_two_likelihoods_for_cond_two = sample_likelihoods[2].to_numpy()[sample2_indices]

        # store the resulting averages
        condition_two_mean_for_neighborhood = np.mean(condition_two_likelihoods_for_cond_two)
        anomaly_scores.append(condition_two_mean_for_neighborhood)
        
    return anomaly_scores

