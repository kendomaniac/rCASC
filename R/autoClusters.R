#' @title Autoencoder clustering
#' @description This function Compress data using autoencoder partially connected
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included.Has to be the one in the projectName folder.Different so from the previous one. 
#' @param nCluster, number of cluster in which the dataset is divided
#' @param projectName, might be different from the matrixname in order to perform different analysis on the same dataset
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param pcaDimensions, number of dimensions to use for Seurat Pca reduction. 
#' @param seed, important value to reproduce the same results with same input
#' @param clusterMethod, clustering methods: "GRIPH","SIMLR","SEURAT","SHARP"
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return 
#' @examples
#' \dontrun{
#'  autoencoderClustering(group="docker", scratch.folder="/home/user/Riccardo/Riccardo/1_inDocker_2/scratch", file="/home/user/Riccardo/Riccardo/1_inDocker_2/data/Results/testDocker/setA.csv",separator=",", nCluster=5,clusterMethod=c("SEURAT"),seed=1111,"testDocker",13)

#'}
#' @export
autoencoderClustering <- function(group=c("sudo","docker"), scratch.folder, file,separator, nCluster,clusterMethod=c("GRIPH","SIMLR","SEURAT","SHARP"),seed=1111,projectName,pcaDimensions){

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

system(paste("cp -r ",data.folder," ",scrat_tmp.folder,"/",sep=""))

  #executing the docker job
  if(clusterMethod=="GRIPH"){
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch:Z -v ", data.folder, ":/data -d docker.io/repbioinfo/permutationclustering Rscript /home/mainGRIPH.R ",projectName," ",matrixName," ",separator," ",nCluster," ",seed,sep="")
}

  if(clusterMethod=="SIMLR"){
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch:Z -v ", data.folder, ":/data -d docker.io/repbioinfo/permutationclustering Rscript /home/mainSIMLR.R ",projectName," ",matrixName," ",separator," ",nCluster," ",seed,sep="")
}

  if(clusterMethod=="SHARP"){
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch:Z -v ", data.folder, ":/data -d docker.io/repbioinfo/permutationclusteringsharp Rscript /home/mainSHARP.R ",projectName," ",matrixName," ",separator," ",nCluster," ",seed,sep="")
}

 if(clusterMethod=="SEURAT"){
    params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,":/scratch:Z -v ", data.folder, ":/data -d docker.io/repbioinfo/seuratpermutation Rscript /home/mainSeurat.R ",projectName," ",matrixName," ",separator," ",nCluster," ",pcaDimensions," ",seed,sep="")
}
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
  system(paste("rm ",scrat_tmp.folder,"/*.pdf ",sep=""))
  system(paste("cp -r ",scrat_tmp.folder,"/* ",data.folder,"/../",sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}
