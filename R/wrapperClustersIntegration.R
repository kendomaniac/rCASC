#' @title integrationCircos
#' @description This function execute toprnk analysis which search for correspondence between clusters of two different experiments using clusters-pseudobulks, zscored on rows, and the cluster specific genes from comet analysis. Thus, the function clustersBulk and cometsc have to be run in the two datasets before their comparison.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file.matrix1, a character string indicating the path of the first matrix 
#' @param file.matrix2, a character string indicating the path of the second matrix 
#' @param separator1, separator used in count file, e.g. '\\t', ','
#' @param separator2, separator used in count file, e.g. '\\t', ','
#' @param cl1, path of clustering.output for file.matrix1
#' @param cl2, path of clustering.output for file.matrix2
#' @param permutation, number of permutation to be run
#' @param seed, integer file necessary for reproducibility
#' @param file.pblk1, a character string indicating the path of the pseudobulkRow file, with file name and extension included. 
#' @param file.pblk2, a character string indicating the path of the pseudobulkRow file, with file name and extension included. 
#' @param comet.folder1, path of Comet outputdata folder from experiment 1
#' @param comet.folder2, path of Comet outputdata folder from experiment 2
#' @param top.ranked, MAX number of top comet genes to be used for each cluster, default 320
#' @param file.total1, a character string indicating the path to the total.csv.for the 1st dataset to be integrated. Total.csv is generated with autoencoder4pseudoBulk. File, with file name and extension included. 
#' @param file.total2, a character string indicating the path to the total.csv.for the 2nd dataset to be integrated. Total.csv is generated with autoencoder4pseudoBulk. File, with file name and extension included. 
#' @param gsea, default msigdb.all, which includes all classes. List of the available GSEA classes: c1.all, c2.cgp, c2.cp.biocarta, c2.cp.kegg, c2.cp.pid, c2.cp.reactome, c2.cp.wikipathways, c3.all, c3.mir, c3.tft.gtrd, c3.tft, c4.cgn, c4.cm, c5.go.bp, c5.go.cc, c5.go.mf, c5.hpo, c6.all, c7.all, c8.all, h.all, msigdb.all. Please note that msigdb.all includes all gsea classes.
#' @param X, X parameter for the XLmHG, default 5, for more info please see XLmHG help: https://xl-mhg.readthedocs.io/en/latest/.
#' @param L, L parameter for the XLmHG, default 0.15, for more info please see XLmHG help: https://xl-mhg.readthedocs.io/en/latest/.
#' @param pvalue, XLmHG pvalue threshold,default 0.05
#' @param outputFolder, where results are placed
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return A picture called integrated_score.png and a file called integrated_score.csv and all the final_scores.csv used to produce the integrated results. 
#' @examples
#' \dontrun{
#'  library(rCASC)
#'  wrapperClustersIntegration(group="docker", 
#'         scratch.folder="/scratch", 
#'         file.matrix1="/data/clusters_association_paper/setA1_set1/setA1/VandE/VandE.csv",
#'         file.matrix2="/data/clusters_association_paper/setA1_set1/set1/VandE/VandE.csv", separator1=",",separator2=",",
#'         cl1="/data/clusters_association_paper/setA1_set1/setA1/VandE/Results/VandE/5/VandE_clustering.output.csv",
#'         cl2="/data/clusters_association_paper/setA1_set1/set1/VandE/Results/VandE/4/VandE_clustering.output.csv",
#'         file.pblk1="/data/clusters_association_paper/setA1_set1/setA1/VandE/VandE_bulkRow.csv",
#'         file.pblk2="/data/clusters_association_paper/setA1_set1/set1/VandE/VandE_bulkRow.csv", 
#'         comet.folder1="/data/clusters_association_paper/setA1_set1/setA1/VandE/Results/VandE/5/outputdata",
#'         comet.folder2="/data/clusters_association_paper/setA1_set1/set1/VandE/Results/VandE/4/outputdata", 
#'         file.total1="/data/clusters_association_paper/setA1_set1/setA1/VandE/Results/setA1/permutation/total.csv",
#'         file.total2="/data/clusters_association_paper/setA1_set1/set1/VandE/Results/set1/permutation/total.csv",
#'         permutation=100, seed=111, top.ranked=320, gsea="msigdb.all", X=5, L=0.15, pvalue=0.05,
#'         outputFolder="/data/clusters_association_paper/setA1_set1"
#'         )
#'}
#' @export
wrapperClustersIntegration <- function(group=c("sudo","docker"), scratch.folder, 
                                       file.matrix1, file.matrix2, 
                                       file.pblk1, file.pblk2, 
                                       file.total1, file.total2, 
                                       comet.folder1, comet.folder2,
                                       cl1, cl2, 
                                       separator1, separator2, 
                                       permutation=100, seed=111,
                                       top.ranked=320, 
                                       gsea="msigdb.all", X=5, L=0.15, pvalue=0.05,
                                       outputFolder){
  
  cat("\nRunning gseaXLmHG\n")
  gseaXLmHG <- function(group=group, scratch.folder=scratch.folder, xCometFolder=comet.folder1, yCometFolder=comet.folder1, 
                        gsea=gsea, X=X, L=L, pvalue=pvalue, separatorX=separator1, separatorY=separator1, outputFolder=outputFolder)
  
  cat("\nRunning seuratIntegrationPermutation\n")
  seuratIntegrationPermutation(group=group, scratch.folder=scratch.folder, 
                               file1=file.matrix1, file2=file.matrix2, separator1=separator1, 
                               separator2=separator1, cl1=cl1, cl2=cl2, permutation=permutation, 
                               seed=seed, outputFolder=paste(outputFolder, "ISC", sep="/"))
  
  cat("\nRunning toprnk\n")
  toprnk(group=group, scratch.folder=scratch.folder, fileX=file.pblk1, fileY=file.pblk2, 
         separatorX=separator1, separatorY=separator2, xCometFolder=comet.folder1, yCometFolder=comet.folder2, 
         top.ranked=top.ranked, outputFolder=outputFolder)
  
  cat("\nRunning integrationPblkae\n")
  integrationPblkae(group=group, scratch.folder=scratch.folder, fileX=file.total1, fileY=file.total2, outputFolder=outputFolder)
  
  cat("\nRunning integrationCircos\n")
  gsea.file=paste(outputFolder, "GSEA/final_score.csv", sep="/")
  isc.file=paste(outputFolder, "ISC/final_score.csv", sep="/")
  XYpb.file=paste(outputFolder, "/XYpb/XYpb_final_score.csv", sep="/")
  pblkae.file=paste(outputFolder, "pblkae/final_score.csv", sep="/")
  
  integrationCircos(group=group, scratch.folder=scratch.folder, gsea.file=gsea.file, isc.file=isc.file, XYpb.file=XYpb.file, pblkae.file=pblkae.file, outputFolder=outputFolder)
  
}
 
