# PLOGS(Pseudotemporal Locating and Ordering of Gene by Single-cell)

PLOGS is a tool to identify marker genes for cell clusters, calculate the pseudotime of genes and profile genes map accoding to cell map. It can be directly used in the [scanpy](https://scanpy.readthedocs.io/en/latest/) workflow. 

<p align="center"><img src="figures/PLOGS.png" alt="PLOGS" width="50%"></p>

PLOGSconsists of three core functions: 
* (i) get_DEG_single to find genes that are differentially expressed in only one cell type based on normalized raw counts of scRNA-seq data; * (ii) get_DEG_multiple to find genes that are differentially expressed in one or more  cell types based on normalized raw counts of scRNA-seq data; 
* (iii) get_genes_pseudotime_location to calculate the pseudotemporal expression of individual genes based on pseudotime ordering of cells and to locate genes into suitable coordinates based on the cells’ locations. 

