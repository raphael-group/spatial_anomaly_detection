# Sardine - A method for detecting anomalous regions in spatial transcriptomics data 

Sardine is a method that quantifies how anomalous a particular region is given spatial transcriptomics (ST) data between two or more different conditions. 

The input to Sardine are the expression matrices and the spatial coordinates matices for the ST data from each condition. For ease of notation, we refer to one condition as the "healthy" slice and the other as a "perturbed" slice, but in practice the conditions are general and not restricted to comparing disesaed tissue. The spatial coordinates are assumed to be spatially registered onto the same common coordinate system before being passed into Sardine. Common coordinate system mapping can be done using the method of the user's choice. We've tested STalign (https://github.com/JEFworks-Lab/STalign) (Clifton, K., Anant, M., Aihara, G. et al. STalign: Alignment of spatial transcriptomics data using diffeomorphic metric mapping. Nat Commun 14, 8123 (2023). https://doi.org/10.1038/s41467-023-43915), but in practice any method will work. 

As output, Sardine returns an anomaly score for each spot in the dataset, corresponding to how anomalous the spatial region centered around that spot by comparing the spots in the spatial region on the "perturbed" slice to the spots in the spatial region in the "healthy" slice. The way Sardine does this is by estimating spatially restricted probability densities of expression for each condition, and then using these labels to compute the average likelihood of the condition given the observed expression in the perturbed spots. For a more detailed explanation, please see our manuscript. 

This repository contains 3 directories:
- `/data` which contains the synthetic data used (as well as links to the real data used in the manuscript)
- `/scripts` which contain code necessary to run Sardine as well as example notebooks on how to do so.
- `/results` which contain a record of the outputted results for easy comparison

each of these subdirectories also contains a README explaining the file contents.
