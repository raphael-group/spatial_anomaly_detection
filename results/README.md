# Results
This directory contains all of the results illustrated in the paper

This directory is broken up further into several more subdirectories
`/GSEA_outputs` contains the results of running the gene set enrichment analysis using webGestalt (https://2024.webgestalt.org/).

Alternatively, you can look at the GSEA outputs using the following links:
Stab wound analysis
	- STANDS https://2024.webgestalt.org/results/1737569257/#
	- Vespucci: https://2024.webgestalt.org/results/1737569712/#
	- Sardine: https://2024.webgestalt.org/results/1737569950/#

spinal cord GSEA:
	- Sardine: https://2024.webgestalt.org/results/1737570684/#
	- Vespucci: https://2024.webgestalt.org/results/1737570335/#

 MELD's results are absent as GSEA did not return any significant (FDR < 0.05) results. 

`/simulation_anomaly_scores` contains the anomaly scores outputted by Sardine and other methods on the simulation data.

`/visium_anomaly_scores` contains the anomaly scores outputted by Sardine and other methods on real data.

In each case, "anomaly score" reviews to the output of each method that quantifies how anomalous a particular spot is. All methods output a score between 0 and 1, where 1 implies an anomaly. The data is stored as an array of length equal to the number of spots in the perturbed condition, where the indices represent the same order as the original data.  
