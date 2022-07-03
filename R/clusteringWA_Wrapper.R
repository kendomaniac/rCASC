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
#' @param bias, bias method to use : "mirna" , "TF", "CUSTOM", kinasi,immunoSignature, cytoBands,ALL 
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
#' @param Sp, minimun number of percentage of cells that has to be in common between two permutation to be the same cluster.
#' @param version, version 1 implements static batchsize, version 2 implements adaptive batchsize 
#' @return folders the complete autoencoder analysis.
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @examples
#' \dontrun{
#'  clusteringWA_Wrapper(group=c("sudo"),scratch.folder="/home/user/autoencoderClustering/test/Scratch/",file="/home/user/autoencoderClustering/test/Data/setA.csv",separator=",",nCluster=3,bias=c("ALL"),permutation=10,nEpochs=10,patiencePercentage=5,seed=1111,projectName="TEST2",bN="NULL",lr=0.01,beta_1=0.9,beta_2=0.999,epsilon=0.00000001,decay=0.0,loss="mean_squared_error",clusterMethod=c( "GRIPH"),pcaDimensions=5,Sp=0.8, version=2)
#'}
#' @export


clusteringWA_Wrapper=function(group=c("sudo","docker"),scratch.folder,file,separator,nCluster,bias=c("mirna","TF", "CUSTOM","kinasi","immunoSignature", "cytoBands", "ALL"),permutation,nEpochs,patiencePercentage=5,seed=1111,projectName,bN="NULL",lr=0.01,beta_1=0.9,beta_2=0.999,epsilon=0.00000001,decay=0.0,loss="mean_squared_error",clusterMethod=c( "GRIPH","SIMLR","SEURAT","SHARP"),pcaDimensions=5,Sp=0.8, version=2){

  if(version == 1){ 
    fileTemp=file

    data.folder=dirname(file)
    positions=length(strsplit(basename(file),"\\.")[[1]])
    matrixNameC=strsplit(basename(file),"\\.")[[1]]
    matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
    format=strsplit(basename(basename(file)),"\\.")[[1]][positions]

    autoencoder4clustering(group=group, scratch.folder=scratch.folder, file=file,separator=separator, bias=bias, permutation=permutation, nEpochs=nEpochs,patiencePercentage=patiencePercentage,seed=seed,projectName=projectName,bN=bN,lr=lr,beta_1=beta_1,beta_2=beta_2,epsilon=epsilon,decay=decay,loss=loss,regularization=regularization, version=1)

    file=paste(dirname(file),"/Results/",projectName,"/",basename(file),sep="")
    autocluster4clustering(group=group, scratch.folder=scratch.folder, file=file, separator=separator, nCluster=nCluster, projectName=projectName, clusterMethod=clusterMethod,seed=seed,pcaDimensions=pcaDimensions)

    projectName2=paste(projectName,"_",clusterMethod,sep="")
    file2=paste(dirname(fileTemp),"/Results/",projectName,"_",clusterMethod,"/",basename(fileTemp),sep="")

    print(projectName2)
    autoencoderAnalysis4Clustering(group=group, scratch.folder=scratch.folder, file=file2,separator=separator,projectName=projectName2,seed=seed,Sp=Sp)
  }else if(verison == 2){
    fileTemp=file
    
    data.folder=dirname(file)
    positions=length(strsplit(basename(file),"\\.")[[1]])
    matrixNameC=strsplit(basename(file),"\\.")[[1]]
    matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
    format=strsplit(basename(basename(file)),"\\.")[[1]][positions]
    
    autoencoder4clustering(group=group, scratch.folder=scratch.folder, file=file,separator=separator, bias=bias, permutation=permutation, nEpochs=nEpochs,patiencePercentage=patiencePercentage,seed=seed,projectName=projectName,bN=bN,lr=lr,beta_1=beta_1,beta_2=beta_2,epsilon=epsilon,decay=decay,loss=loss,regularization=regularization, version=2)
    
    file=paste(dirname(file),"/Results/",projectName,"/",basename(file),sep="")
    autocluster4clustering(group=group, scratch.folder=scratch.folder, file=file, separator=separator, nCluster=nCluster, projectName=projectName, clusterMethod=clusterMethod,seed=seed,pcaDimensions=pcaDimensions)
    
    projectName2=paste(projectName,"_",clusterMethod,sep="")
    file2=paste(dirname(fileTemp),"/Results/",projectName,"_",clusterMethod,"/",basename(fileTemp),sep="")
    
    print(projectName2)
    autoencoderAnalysis4Clustering(group=group, scratch.folder=scratch.folder, file=file2,separator=separator,projectName=projectName2,seed=seed,Sp=Sp)
    
  }

}
