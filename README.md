# CASC: Classification Analysis of Single Cell Sequencing Data
Genome-wide single-cell measurements such as transcriptome sequencing enable the characterization of cellular composition as well as functional variation in homogenic/heterogenic cell populations. An important step in the single-cell transcriptome analysis is to group cells that belong to the same sub-type based on gene expression patterns [1-3]. Critical issues in cell clustering are (i) cluster stability and (ii) feature selection, i.e. the identification genes playing the major role in cluster formation. To address the above issues, we have developed CASC, an analysis workflow implemented in docker containers. CASC uses as core application to detect cell clusters the “kernel based similarity learning” [4] and allows: (i) identification of the optimal number of clusters for cell partitioning. (ii) The evaluation of clusters stability, measuring the permanence of a cell in a cluster upon random removal of subsets of cells. 
CASC was tested on previously published data sets [1-3, ....]. 

References
1] Usoskin et al. Nat. Neurosci. 2014, 18:145–153
[2] Pollen et al. Nat. Biotechnol. 2014, 32:1–37
[3] Kolodziejczyk et al. Cell Stem Cell 2015, 17:471–485
[4] Wang et al. Nat Methods. 2017 14:414-416

## Installation

```
install.packages("devtools")
library(devtools)
install_github("kendomaniac/casc", ref="master")
```

## Requirements

You need to have docker installed on your linux machine, for more info see this document: https://docs.docker.com/engine/installation/. 

The functions in CASC package require that user is sudo or part of a docker group. See the following document for more info: https://docs.docker.com/engine/installation/linux/ubuntulinux/#/manage-docker-as-a-non-root-user

IMPORTANT The first time casc is installed the downloadContainers needs to be executed to download to the local repository the containers that are needed for the use of docker4seq

```
library(casc)
downloadContainers()
```

More info: [CASC vignette](http://rpubs.com/rcaloger/285423)

