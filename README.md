### Hyperion

Useful scripts for the analysis of hyperion data. 
The segmentation and creation of masks has been precomputed and these scripts work on their output, which is a table of intensities per cell. 
The scripts work in the following order: 

1. preprocess the data and produce plots to examine density plots of the markers and potential batch effects,  
2. cluster the cells using FlowSOM and make all necessary plots for the annotation of the clusters, 
3. cluster the cells using Rphenograph.   
