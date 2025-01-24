# Scripts directory

This directory contains code to run Sardine. 

A simple implementation of Sardine, using the plug in density estimator MELD [1] (https://github.com/KrishnaswamyLab/MELD), and defining the spatial anomaly family as balls around each spot.

The key function is `calculate_anomaly_score()`, which takes in the following parameters: 

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

We also include a function `create_spatialKDTrees` to create the spatial KDTrees list. An example script is also included: `example_notebook.ipynb`.

Sardine is under active development, with improvements to come.

[1] Burkhardt, D.B., Stanley, J.S., Tong, A. et al. Quantifying the effect of experimental perturbations at single-cell resolution. Nat Biotechnol 39, 619–629 (2021). https://doi.org/10.1038/s41587-020-00803-5
