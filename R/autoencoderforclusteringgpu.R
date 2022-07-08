#' @title Autoencoder4clustering
#' @description The present function compress data using autoencoder partially connected
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param permutation, number of permutations to perform the pValue to evaluate clustering
#' @param nEpochs, number of Epochs for neural network training
#' @param projectName, might be different from the matrixname in order to perform different analysis on the same dataset
#' @param patiencePercentage, number of Epochs percentage of not training before to stop.
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param bias, bias method to use : "mirna" , "TF", "CUSTOM", kinasi,immunoSignature, cytoBands,ALL 
#' @param bN, name of the custom bias file. This file need header, in the first column has to be the source and in the second column the gene symbol. All path needs to be provided. 
#' @param seed, important value to reproduce the same results with same input
#' @param lr, learning rate, the speed of learning. Higher value may increase the speed of convergence but may also be not very precise
#' @param beta_1, look at keras optimizer parameters
#' @param beta_2, look at keras optimizer parameters 
#' @param epsilon, look at keras optimizer parameters
#' @param decay, look at keras optimizer parameters
#' @param loss, loss of function to use, for other loss of function check the keras loss of functions. 
#' @param regularization, this parameter balances between reconstruction loss and enforcing a normal distribution in the latent space.
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @examples
#' \dontrun{
#autoencoder4clusteringGPU(group=c("sudo"), scratch.folder=scratch, file=file,separator=",", bias="ALL", permutation=10, nEpochs=100,patiencePercentage=5,seed=1111,projectName=projectName,bN="NULL",lr=0.01,beta_1=0.9,beta_2=0.999,epsilon=0.00000001,decay=0.0,loss="mean_squared_error",regularization=10)
#'}
#' @export
autoencoder4clusteringGPU <- function(group=c("sudo","docker"), scratch.folder, file,separator, bias, permutation, nEpochs,patiencePercentage=5,seed=1111,projectName,bN="NULL",lr=0.01,beta_1=0.9,beta_2=0.999,epsilon=0.00000001,decay=0.0,loss="mean_squared_error",regularization=10){

  data.folder=dirname(file)
positions=length(strsplit(basename(file),"\\.")[[1]])
matrixNameC=strsplit(basename(file),"\\.")[[1]]
matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
format=strsplit(basename(basename(file)),"\\.")[[1]][positions]


  #running time 1
  ptm <- proc.time()
  #setting the data.folder as working folder
  if (!file.exists(data.folder)){
    cat(paste("\nIt seems that the ",data.folder, " folder does not exist\n"))
    return(2)
  }

  #storing the position of the home folder
  home <- getwd()
  setwd(data.folder)
  #initialize status
  system("echo 0 > ExitStatusFile 2>&1")

  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    system("echo 10 > ExitStatusFile 2>&1")
    setwd(home)
    return(10)
  }



  #check  if scratch folder exist
  if (!file.exists(scratch.folder)){
    cat(paste("\nIt seems that the ",scratch.folder, " folder does not exist\n"))
    system("echo 3 > ExitStatusFile 2>&1")
    setwd(data.folder)
    return(3)
  }
  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folder=file.path(scratch.folder, tmp.folder)
  writeLines(scrat_tmp.folder,paste(data.folder,"/tempFolderID", sep=""))
  cat("\ncreating a folder in scratch folder\n")
  dir.create(file.path(scrat_tmp.folder))
  #preprocess matrix and copying files

if(separator=="\t"){
separator="tab"
}

system(paste("cp ",data.folder,"/",matrixName,".",format," ",scrat_tmp.folder,"/",sep=""))
  if(bias=="CUSTOM"){
system(paste("cp ",bN," ",scrat_tmp.folder,"/",sep=""))
}
  bN=basename(bN)
  #executing the docker job
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ", data.folder, ":/data --gpus all -d repbioinfo/autoencoderforclusteringgpu python3 /home/autoencoder.py ",matrixNameC,".",format," ",separator," ",bias," ",permutation," ",nEpochs," ",patiencePercentage," ",projectName," ",seed," ",bN," ",lr," ",beta_1," ",beta_2," ",epsilon," ",decay," ",loss,sep="")

resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){
    #system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
  }
  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(data.folder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
    con <- file("run.info", "r")
    tmp.run <- readLines(con)
    close(con)
    tmp.run[length(tmp.run)+1] <- paste("user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))


  #Copy result folder
  cat("Copying Result Folder")
  system(paste("cp -r ",scrat_tmp.folder,"/* ",data.folder,sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}  
