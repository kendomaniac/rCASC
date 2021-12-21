#' @title wrapperMixModelsUmap
#' @description This function executes a ubuntu docker that performs wrapperMixModelsUmap
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the count matrix
#' @param geneList, a character string indicating the path of the geneList matrix (no header, no row names)
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param k, number of soft labels to divide the density
#' @param seed, number for reproducibility
#' @param epochs, epochs number for umap algorithm
#' @param finalName, name for output files
#' @param maxit, number of max iteration for mixmodels
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return plot
#' @examples
#'\dontrun{
#' dir.create("scratch")
#' wrapperMixModelsUmap(group=c("docker"), scratch.folder=paste(getwd(),"/scratch",sep=""), file=paste(getwd(),"/setA.csv",sep=""),separator=",",seed=111,epochs=1000,k=3,finalName="prova",geneList=paste(getwd(),"/geneList.csv",sep=""))
#' }
#' @export
wrapperMixModelsUmap <- function(group=c("sudo","docker"), scratch.folder, file,geneList,separator,k,seed,epochs,finalName,maxit=10000){
  dir.create("scratch")
  
  
  umap(group=group, scratch.folder=scratch.folder, file=file,separator=separator,seed=seed,epochs=epochs)
  name=tools::file_path_sans_ext(basename(file))
  geneVisualization(group=group, scratch.folder=scratch.folder, file=file,clustering.output=paste(getwd(),"/",name,"_fake_clustering.output.csv",sep=""),geneList=geneList,separator=separator,finalName=finalName)
  mixmodels(group=group, scratch.folder=scratch.folder, file=file,geneList=geneList,separator=separator,k=k,maxit=maxit)
}
