#' @title Seurat Integration Permutation
#' @description This function executes a ubuntu docker that performs seurat integration to identify associated clusters in two independent experiment. The analysis is repeated mutiple times removing 10% of the initial cells, to investigate how stable is the partitioning done with seuratIntegration function. This function requires that cells were clustered with any of the clustering tools included in rCASC.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file1, a character string indicating the path of the first matrix 
#' @param file2, a character string indicating the path of the second matrix 
#' @param separator1, separator used in count file, e.g. '\\t', ','
#' @param separator2, separator used in count file, e.g. '\\t', ','
#' @param cl1, path of clustering.output for file1
#' @param cl2, path of clustering.output for file2
#' @param permutation, number of permutation for statistic
#' @param seed, integer file necessary for reproducibility
#' @param outputFolder, path to the output folder
#' @param K, resolution for seurat analysis

#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#' @return A folder called ISC, the input data, the intermediate results and a comma separated file final_score.csv summarising the frequency by which the seuratIntegration function return a an integration cluster including at least 50% of the input clusters form X and Y datasets.
#' @examples
#' \dontrun{
#' seuratIntegrationPermutation(group="docker", scratch.folder="/home/user/scratch", file1="/home/user/dockerFile/Seurat_join_DAPUSHARE/function/example/set1.csv",file2="/home/user/dockerFile/Seurat_join_DAPUSHARE/function/example/setA.csv", separator1=",",separator2=",",cl1=, cl29,permutation=100, seed=111) 
#'}
#' @export
seuratIntegrationPermutation <- function(group=c("sudo","docker"), scratch.folder, file1, file2, separator1, separator2,cl1,cl2,permutation,seed,outputFolder,K=0.8){
dir.create(outputFolder)
  data.folder1=dirname(file1)
positions1=length(strsplit(basename(file1),"\\.")[[1]])
matrixNameC1=strsplit(basename(file1),"\\.")[[1]]
matrixName1=paste(matrixNameC1[seq(1,positions1-1)],collapse="")
format1=strsplit(basename(basename(file1)),"\\.")[[1]][positions1]


  data.folder2=dirname(file2)
positions2=length(strsplit(basename(file2),"\\.")[[1]])
matrixNameC2=strsplit(basename(file2),"\\.")[[1]]
matrixName2=paste(matrixNameC2[seq(1,positions2-1)],collapse="")
format2=strsplit(basename(basename(file2)),"\\.")[[1]][positions2]
data.folder=data.folder1
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

if(separator1=="\t"){
separator1="tab"
}
if(separator2=="\t"){
separator2="tab"
}
system(paste("cp ",data.folder1,"/",matrixName1,".",format1," ",scrat_tmp.folder,"/X_",matrixName1,".",format1,sep=""))
system(paste("cp ",data.folder2,"/",matrixName2,".",format2," ",scrat_tmp.folder,"/Y_",matrixName2,".",format2,sep=""))
system(paste("cp ",cl1," ",scrat_tmp.folder,"/X_",basename(cl1),sep=""))
system(paste("cp ",cl2," ",scrat_tmp.folder,"/Y_",basename(cl2),sep=""))
cl1=paste(strsplit(basename(cl1),"[.]")[[1]][seq(1,length(strsplit(basename(cl1),"[.]")[[1]])-1)],collapse=".")
cl2=paste(strsplit(basename(cl2),"[.]")[[1]][seq(1,length(strsplit(basename(cl2),"[.]")[[1]])-1)],collapse=".")
  #executing the docker job
    params <- paste("--cidfile ",data.folder1,"/dockerID -v ",scrat_tmp.folder,":/scratch -v ", data.folder1, ":/data -d docker.io/repbioinfo/seuratintegrationpermutation Rscript /home/pre_processing.R X_",matrixName1," ",format1," ",separator1," Y_",matrixName2," ",format2," ",separator2," X_",cl1," Y_",cl2," ",permutation," ",seed," ",K,sep="")

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
  system(paste("docker logs ", substr(container.id,1,12), " >& ",data.folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))


  #Copy result folder
  cat("Copying Result Folder")
  system(paste("cp -r ",scrat_tmp.folder,"/* ",outputFolder,sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}   
