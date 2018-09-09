#' @title Bootstraps vide
#' @description This function generate a video showing the relocation of cells during the bootstraps
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param nCluster, which nCluster results to use for this analysis
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param framePP, Number of frames for each bootstrap
#' @param permutationNumber, Number of random bootstraps, should be less or equal to the total number of permutations used to generate tha clustering data

#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return a video
#' @examples
#'\dontrun{
#'  bootstrapsVideo(group,scratch.folder,file,nCluster,separator)#
#'}
#' @export
bootstrapsVideo <- function(group=c("sudo","docker"), scratch.folder, file, nCluster, separator, framePP=200, permutationNumber){

  cat("\nStart cluster identification\n")
  clusterIdentification(group=group, scratch.folder=scratch.folder, file=file, nCluster=nCluster, separator=separator)
  cat("\nStart movie generation\n")
  permutationMovie(group=group, scratch.folder=scratch.folder, file=file, nCluster=nCluster, separator=separator, framePP=framePP, permutationNumber=permutationNumber)

}
