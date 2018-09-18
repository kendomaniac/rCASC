#' @title Gene prioritization with SIMLR
#' @description This function executes a ubuntu docker that produces a specific number of permutation to evaluate clustering and identify the genes that play the major role in clustering.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param nPerm, number of permutations to be executed
#' @param permAtTime, number of permutations computed in parallel
#' @param percent, percentage of random cells removed in each permutation
#' @param nCluster, the number of clusters, where to run prioritization
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param logTen, 1 if the count matrix is already in log10, 0 otherwise
#' @param seed, important value to reproduce the same results with same input
#' @param sp, minimun number of percentage of cells that has to be in common in a cluster, between two permutations, default 0.8
#' @param clusterPermErr, probability error in depicting the number of clusters in each permutation, default = 0.05
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#' @return VioPlot of silhouette cells value for each number of cluster used,clusterP file with clustering results for each permutation, killedCell file with removed cells in each permutation, clustering.output a sommarize file with general information for each cells.
#' @examples
#' \dontrun{
#' system("wget http://130.192.119.59/public/section4.1_examples.zip")
#' unzip("section4.1_examples.zip")
#' setwd("section4.1_examples")
#' simlrFeatures(group="docker",scratch.folder="/data/scratch/",file=paste(getwd(), "bmsnkn_5x100cells.txt", sep="/"), nPerm=160, permAtTime=8, percent=10, nCluster=5, separator="\t", logTen=0, seed=111, sp=0.8, clusterPermErr=0.05)
#'}
#' @export

simlrFeatures <- function(group=c("sudo","docker"), scratch.folder, file, nPerm, permAtTime, percent, nCluster, separator, logTen=0, seed=111, rK=1, sp=0.8, clusterPermErr=0.05){
  
  permutationClustering(group=group, scratch.folder=scratch.folder, file=file, nPerm=nPerm, permAtTime=permAtTime, percent=percent, range1=nCluster, range2=nCluster, separator=separator, logTen=logTen, clustering="SIMLR", perplexity=10 , seed=seed, rK=1)
  
  
}
