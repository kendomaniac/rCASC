# CASC
CASC: Classification Analysis of Single Cell Sequencing Data

Genome-wide single-cell measurements such as transcriptome sequencing enable the characterization of cellular composition as well as functional variation in homogenic/heterogenic cell populations. An important step in the single-cell transcriptome analysis is to group cells that belong to the same sub-type based on gene expression patterns [1-3]. Critical issues in cell clustering are (i) cluster stability and (ii) feature selection, i.e. the identification genes playing the major role in cluster formation. To address the above issues, we have developed CASC, a tool implemented in a docker container, that uses as core application to detect cell clusters the “kernel based similarity learning” [4] and allows: (i) identification of the optimal number of clusters for cell partitioning using “silhouette method”. (ii) The evaluation of clusters stability, measuring the permanence of a cell in a cluster upon random removal of subsets of cells. (iii) Feature selection via “nearest shrunken centroid method” [5], applied to the gene Index Of Dispersion [6]. CASC was tested on previously published data sets [1-3]. CASC feature selection procedure efficiently allows the identification of the subpopulation of genes playing the main role in discriminating   between cell subpopulations. 

References
<br>[1] Usoskin et al. Nat. Neurosci. 2014, 18:145–153
<br>[2] Pollen et al. Nat. Biotechnol. 2014, 32:1–37
<br>[3] Kolodziejczyk et al. Cell Stem Cell 2015, 17:471–485
<br>[4] Wang et al. http://biorxiv.org/content/early/2016/05/09/052225
<br>[5] Tibshirani et al. PNAS 2002, 99:6567-6572
<br>[6] Diaz et al. Bioinformatics 2016, 32: 2219-2220

