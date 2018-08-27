#' @title Permutations and Clustering
#' @description This function executes a ubuntu docker that produces a specific number of permutation to evaluate the range of optimal number of clusters using griph algorithm. For more info see CSC vignette.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param nPerm, number of permutations to perform the pValue to evaluate clustering
#' @param permAtTime, number of permutations that can be computes in parallel
#' @param percent, percentage of random cells that has to be removed in each permutation
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param logTen, 1 if the count matrix is already in log10, 0 otherwise
#' @param seed, important value to reproduce the same results with same input

#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return VioPlot of silhouette cells value for each number of cluster used,clusterP file with clustering results for each permutation, killedCell file with removed cells in each permutation, clustering.output a sommarize file with general information for each cells.
#' @examples
#' \dontrun{
#'  clusterNgriph("docker",scratch.folder="/data/scratch/",file="/data/test_griph/b.m.s.nk.n_5x100cells.txt", nPerm=16, permAtTime=8, percent=10, separator="\t",logTen=0, seed=111)
#'}
#' @export
clusterNgriph <- function(group=c("sudo","docker"), scratch.folder, file, nPerm, permAtTime, percent, separator, logTen=0, seed=111){

  permutationClustering(group=group, scratch.folder=scratch.folder, file=file, nPerm=nPerm, permAtTime=permAtTime, percent=percent,range1=1,
                        range2=1, separator=separator, logTen=logTen, clustering="griph", perplexity=0, seed=seed, rK=0)

}
