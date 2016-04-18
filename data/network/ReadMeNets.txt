Readme for six biological networks:

Three PPI networks:

STRING: human PPIs with high confidence score (>600) from STRING, download from http://string-db.org/. Mapping protein ids to gene ids based on biomart. File: STRINGnetmap.txt

HPRD: STRING + HPRD PPI network used in MAGI (The discovery of integrated gene networks for autism and related disorders, http://genome.cshlp.org/content/25/1/142.short).
The PPI network is composed of the union of StringDB v9.05 (Szklarczyk et al. 2011) human (organism ID 9606) interactions that are experimentally verified (experimental scores > 400) and have high confidence scores (> 700), together with the complete HPRD database (http://www.hprd.org/). File: StringNew_HPRD_mnet.txt

iRefIndex: http://irefindex.org/wiki/index.php?title=iRefIndex. File: iRefIndexm.txt

Two co-expression networks were built based on the brain development expression data from Brainspan (http://www.brainspan.org/). The neocortical regions expression data is focused on mid-fetal development in 10-24 post-conception weeks with 139 samples. 

CORR: The first one is built based on a threshold to select edges from the pairwise Pearson correlation matrix. The correlation threshold r = 0.7 is used here to select edges. File: brainspan_net_cor.txt

CoEXP: a robust rank-based network construction method (A general co-expression network-based approach to gene expression analysis: comparison and applications. BMC systems biology, 2010) is used, named as CoEXP. Top five correlated genes are used as the connected neighbors for each gene. File: brainspan_net_top5.txt

CoPrePPI: combined the high-confidence predictions of domain-motif mediated interactions from PrePPI (https://honiglab.c2b2.columbia.edu/PrePPI/ref/DomainMotif/DomainMotif_new_predictions.txt) and CoEXP. Edges in CoPrePPI are the union of two integral networks. File: ComCo_PrePPI.txt

All network node and edge betweenness centralities are computed based on R package igraph with functions: betweenness; edge.betweenness

