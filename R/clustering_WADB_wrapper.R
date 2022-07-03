#' @title wrapperAutoencoder
#' @description This function executes the whole autoencoder pipeline
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param permutation, number of permutations to perform the pValue to evaluate clustering
#' @param nCluster, number of cluster in which the dataset is divided
#' @param nEpochs, number of Epochs for neural network training
#' @param projectName, might be different from the matrixname in order to perform different analysis on the same dataset
#' @param patiencePercentage, number of Epochs percentage of not training before to stop.
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param seed, important value to reproduce the same results with same input
#' @param lr, learning rate, the speed of learning. Higher value may increase the speed of convergence but may also be not very precise
#' @param loss, loss of function to use, for other loss of function check the keras loss of functions. 
#' @param clusterMethod, clustering methods: "GRIPH","SIMLR","SEURAT","SHARP"
#' @param pcaDimensions, number of dimensions to use for Seurat Pca reduction. 
#' @param Sp, minimun number of percentage of cells that has to be in common between two permutation to be the same cluster.
#' @return folders the complete autoencoder analysis.
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @examples
#' \dontrun{
#'  clusteringWADB_Wrapper(group=c("sudo"),scratch.folder="/home/user/scratch",file="/home/user/autoencoderClustering_v4/u/setA.csv",separator=",",nCluster=5,permutation=80,nEpochs=1000,patiencePercentage=5,seed=1111,projectName="yuppy",clusterMethod=c( "SIMLR"),lr=0.001)

#'}
#' @export


clusteringWADB_Wrapper=function(group=c("sudo","docker"),scratch.folder,file,separator,nCluster,permutation,nEpochs,patiencePercentage=5,seed=1111,projectName,lr=0.01,loss="mean_squared_error",clusterMethod=c( "GRIPH","SIMLR","SEURAT","SHARP"),pcaDimensions=5,Sp=0.8){
fileTemp=file

data.folder=dirname(file)
positions=length(strsplit(basename(file),"\\.")[[1]])
matrixNameC=strsplit(basename(file),"\\.")[[1]]
matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
format=strsplit(basename(basename(file)),"\\.")[[1]][positions]

autoencoderDB(group=group, scratch.folder=scratch.folder, file=file,separator=separator, permutation=permutation, nEpochs=nEpochs,patiencePercentage=patiencePercentage,seed=seed,projectName=projectName,lr=lr,loss=loss)

file=paste(dirname(file),"/Results/",projectName,"/",basename(file),sep="")
autocluster4clustering(group=group, scratch.folder=scratch.folder, file=file, separator=separator, nCluster=nCluster, projectName=projectName, clusterMethod=clusterMethod,seed=seed,pcaDimensions=pcaDimensions)

projectName2=paste(projectName,"_",clusterMethod,sep="")
file2=paste(dirname(fileTemp),"/Results/",projectName,"_",clusterMethod,"/",basename(fileTemp),sep="")

print(projectName2)
autoencoderAnalysis4Clustering(group=group, scratch.folder=scratch.folder, file=file2,separator=separator,projectName=projectName2,seed=seed,Sp=Sp)














} 
