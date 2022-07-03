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
#' @param bias, bias method to use : "mirna" , "TF", "CUSTOM", kinasi,immunoSignature,ALL 
#' @param cl, Clustering.output file. Can be the output of every clustering algorithm from rCASC or can be customized with first column cells names, second column cluster they belong.All path needs to be provided.  
#' @param bN, name of the custom bias file. This file need header, in the first column has to be the source and in the second column the gene symbol.All path needs to be provided, 
#' @param seed, important value to reproduce the same results with same input
#' @param lr, learning rate, the speed of learning. Higher value may increase the speed of convergence but may also be not very precise
#' @param beta_1, look at keras optimizer parameters
#' @param beta_2, look at keras optimizer parameters 
#' @param epsilon, look at keras optimizer parameters
#' @param decay, look at keras optimizer parameters
#' @param loss, loss of function to use, for other loss of function check the keras loss of functions. 
#' @param clusterMethod, clustering methods: "GRIPH","SIMLR","SEURAT","SHARP"
#' @param pcaDimensions, number of dimensions to use for Seurat Pca reduction. 
#' @param permAtTime, number of permutation in parallel
#' @param largeScale, boolean for SIMLR analysis, TRUE if rows are less then columns or if the computational time are huge
#' @param Sp, minimun number of percentage of cells that has to be in common between two permutation to be the same cluster.
#' @param threads, integer refering to the max number of process run in parallel default 1 max the number of clusters under analysis, i.e. nCluster
#' @param X, from 0 to 1 argument for XL-mHG default 0.15, for more info see cometsc help.
#' @param K, the number of gene combinations to be considered., possible values 2, 3, 4, default 2. WARNING increasing the number of combinations makes the matrices very big
#' @param counts, if set to True it will graph the log(expression+1). To be used if unlogged data are provided
#' @param skipvis, set to True to skip visualizations
#' @param variational, TRUE or FALSE if you want to use a variational autoencoder or the standard autoencoder
#' @param regularization, this parameter balances between reconstruction loss and enforcing a normal distribution in the latent space
#' @return folders the complete autoencoder analysis.
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @examples
#' \dontrun{
#'  wrapperAutoencoder(group="sudo",scratch.folder=scratch.folder,file="/home/lucastormreig/test/setA.csv",separator=",",nCluster=5,bias="mirna",permutation=10,nEpochs=10,cl="/home/lucastormreig/test/setA_clustering.output.csv",projectName="mirna",clusterMethod="GRIPH")
#'}
#' @export


wrapperAutoencoder=function(group=c("sudo","docker"),scratch.folder,file,separator,nCluster,bias=c("mirna","TF", "CUSTOM","kinasi","immunoSignature","ALL"),permutation,nEpochs,patiencePercentage=5,cl,seed=1111,projectName,bN="NULL",lr=0.01,beta_1=0.9,beta_2=0.999,epsilon=0.00000001,decay=0.0,loss="mean_squared_error",clusterMethod=c( "GRIPH","SIMLR","SEURAT","SHARP"),pcaDimensions=5,permAtTime=3,largeScale=FALSE,Sp=0.8,threads=1,  X=0.15, K=2, counts=c("False"), skipvis=c("False"),regularization=10,variational=FALSE){
fileTemp=file

data.folder=dirname(file)
positions=length(strsplit(basename(file),"\\.")[[1]])
matrixNameC=strsplit(basename(file),"\\.")[[1]]
matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
format=strsplit(basename(basename(file)),"\\.")[[1]][positions]



autoencoder(group=group, scratch.folder=scratch.folder, file=file,separator=separator, nCluster=nCluster, bias=bias, permutation=permutation, nEpochs=nEpochs,patiencePercentage=patiencePercentage, cl=cl,seed=seed,projectName=projectName,bN=bN,lr=lr,beta_1=beta_1,beta_2=beta_2,epsilon=epsilon,decay=decay,loss=loss,regularization=regularization,variational=variational)

file=paste(dirname(file),"/Results/",projectName,"/",basename(file),sep="")
autoencoderClustering(group=group, scratch.folder=scratch.folder, file=file,separator=separator, nCluster=nCluster,clusterMethod=clusterMethod,seed=seed,pcaDimensions=pcaDimensions,permAtTime=permAtTime,largeScale=largeScale)
projectName=paste(projectName,"_",clusterMethod,sep="")
file=paste(dirname(fileTemp),"/Results/",projectName,"/",basename(fileTemp),sep="")
autoencoderAnalysis(group=group, scratch.folder, file,separator, nCluster,seed=seed,Sp=Sp)
autoFeature(group=group, scratch.folder, file,separator, nCluster)
file=paste(dirname(file),"/",nCluster,"/freqMatrix.",format,sep="")
cometsc2(group=group, file=file, scratch.folder=scratch.folder, threads=threads,  X=X, K=K, counts=counts, skipvis=skipvis, nCluster=nCluster, separator=separator,clustering.output=cl)


}
