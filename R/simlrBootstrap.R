#' @title Permutations and Clustering
#' @description This function executes a ubuntu docker that produces a specific number of permutation to evaluate clustering.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param nPerm, number of permutations to perform the pValue to evaluate clustering
#' @param permAtTime, number of permutations that can be computes in parallel
#' @param percent, percentage of random cells that has to be removed in each permutation
#' @param range1, first number of cluster for k means algorithm
#' @param range2, last number of cluster for k means algorithm
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param logTen, 1 if the count matrix is already in log10, 0 otherwise
#' @param seed, important value to reproduce the same results with same input
#' @param rK, 1 for rankGene algorithm 0 otherwise WARNING, very slow with this feature. This parameter is specific for SIMLR. Setting this parameter to 1 with other clustering methods will not give any different result to set the parameter to 0.
#' @param sp, minimun number of percentage of cells that has to be in common between two permutation to be the same cluster, default 0.8.
#' @param clusterPermErr, error that can be done by each permutation in cluster number depicting.Default = 0.05
#' @param maxDeltaConfidence, max value for Delta confidence for genes feature selection, only required if rK=1
#' @param minLogMean, min value for Log mean for genes feature selection, only required if rK=1


#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return VioPlot of silhouette cells value for each number of cluster used,clusterP file with clustering results for each permutation, killedCell file with removed cells in each permutation, clustering.output a sommarize file with general information for each cells.
#' @examples
#' \dontrun{
#' system("wget http://130.192.119.59/public/section4.1_examples.zip")
#' unzip("section4.1_examples.zip")
#' setwd("section4.1_examples")
#' simlrBootstrap(group="docker",scratch.folder="/data/scratch/",file=paste(getwd(), "bmsnkn_5x100cells.txt", sep="/"), nPerm=160, permAtTime=8, percent=10, range1=4, range2=6, separator="\t",logTen=0, seed=111, sp=0.8, clusterPermErr=0.05, maxDeltaConfidence=NULL, minLogMean=NULL)
#'}
#' @export

simlrBootstrap <- function(group=c("sudo","docker"), scratch.folder, file, nPerm, permAtTime, percent, range1, range2, separator, logTen=0, seed=111, rK=0, sp=0.8, clusterPermErr=0.05, maxDeltaConfidence=NULL, minLogMean=NULL){

  if(is.null(maxDeltaConfidence)){
    maxDeltaConfidence="NULL"
  }
  if(is.null(minLogMean)){
    minLogMean="NULL"
  }

  permutationClustering(group=group, scratch.folder=scratch.folder, file=file, nPerm=nPerm, permAtTime=permAtTime, percent=percent, range1=range1, range2=range2, separator=separator, logTen=logTen, clustering="SIMLR", perplexity=10 , seed=seed, rK=rK)
  permAnalysis(group=group, scratch.folder=scratch.folder,file=file, range1=range1, range2=range2, separator=separator, sp=sp, clusterPermErr=clusterPermErr, maxDeltaConfidence=maxDeltaConfidence, minLogMean=minLogMean)

}

